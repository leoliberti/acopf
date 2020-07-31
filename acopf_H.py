#!/usr/bin/env python3

## ACOPF formulations
## acopf_H.py Python file
## Leo Liberti
## This Python script implements the complex SDP relaxation described in
##   ../../papers/acopf/formulations/4or/acopf_survey.pdf
##    it can integrate shunts, transformers or both
##   Formulation in complex numbers
##     decision variables: generated power, voltage
## NOTE: this formulation uses injected current bounds (either rel or restr)
## 191125 work started

import sys
import os.path
from amplpy import AMPL
import cvxpy as cp
import cvxopt
import time
import math
import cmath
import types
import numpy as np
import re
from pprint import pprint

######### CONFIGURABLE PARAMETERS ############

myZero = 1e-8
myLarge = 1e7
myEps = 1e-3
myCEps = myEps*(1+1j)

# fields in a bus record
busRecDesc = { 'busType':0, 'SDR':1, 'SDC':2, 'VL':3, 'VU':4, 'Vm':5, 'Va':6, 'shR':7, 'shC':8 }

# fields in a line record
lineRecDesc = { 'status':0, 'SU':1, 'r':2, 'x':3, 'bb':4, 'tau':5, 'theta':6, 'pdLB':7, 'pdUB':8 }

# fields in a generator record
genRecDesc = { 'SLR':0, 'SLC':1, 'SUR':2, 'SUC':3 }

# indexing starts from 0 (offset=1) or 1 (offset=0)
offset = 1

# offset for indexing of Schur complement matrix
s = 1

# using regexp to extract floats <https://stackoverflow.com/questions/4703390/how-to-extract-a-floating-number-from-a-string>
float_regexp = r"[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?"
floatrx = re.compile(float_regexp, re.VERBOSE)

# trace coefficient for +tr(X) in minimizing objective
#traceCoeff = 0.01
traceCoeff = 0.0

################ FUNCTIONS ###################

# read AMPL data file into a list of triplets (pnames,sets,vals)
#   for syntax "param : sets : pnames := vals"
def readDat(fd):
    p = []
    pset = None
    pname = None
    pval = None
    for line in fd:
        # ignore empty lines or comments
        if len(line) < 1 or line[0] == "#":
            continue
        try:
            line = line.decode('UTF-8').strip()
        except:
            pass
        # line elements
        l = [le.strip() for le in line.split()]
        lenl = len(l)
        if lenl > 0:
            if l[0] == "param":
                # AMPL parameter declaration
                ptype = 'P'
                lparm = " ".join(l[1:])
                ll = [lp.strip() for lp in lparm.split(':')]
                pset = []
                pname = []
                pval = []
                if len(ll) == 2:
                    # syntax [param pname :=]
                    pname.append(ll[0])
                elif len(ll) == 3:
                    # syntax [param : pname1 ... pnameN := ]
                    thenames = [lp.strip() for lp in ll[1].split()]
                    for pn in thenames:
                        pname.append(pn)
                elif len(ll) == 4:
                    # syntax [param : set1 ... setM : pname1 ... pnameN :=]
                    thesets = [sp.strip() for sp in ll[1].split()]
                    for ps in thesets:
                        pset.append(ps)
                    thenames = [lp.strip() for lp in ll[2].split()]
                    for pn in thenames:
                        pname.append(pn)
                if ll[-1][-1] == ';':
                    # end of param declaration on the line
                    pval.append([float(lp) for lp in ll[-1][0:-1].split() if lp != "="])
                    p.append((ptype,pname,pset,pval))
            elif l[0] == "set":
                # AMPL set declaration
                #   each set assumed defined on 1 line, no comments,
                #   final ';' separated from previous set elt by a space
                # "S[i] := j ;" => 'S' in pname, i in pset, j in pval
                # "S := j ;" => 'S' in pname, j in pval
                ptype = 'S'
                pset = []
                pname = []
                pval = []
                lparm = " ".join(l[1:])
                ll = [lp.strip() for lp in lparm.split(':')]
                longsetname = ll[0]
                setname = longsetname
                setindex = -1
                if '[' in longsetname:
                    # set is indexed
                    setname = longsetname.split('[')[0]
                    setindex = int(longsetname.split('[')[1].split(']')[0])
                pname = [setname]
                if setindex >= 0:
                    pset = [setindex]
                pval = [int(lp) for lp in ll[1].split()[1:-1]]
                p.append((ptype,pname,pset,pval))
            elif l[0] == ";":
                # end of current param declaration on a line by itself
                p.append((ptype,pname,pset,pval))
            else:
                # numerical data
                if line[-1] == ';':
                    # end of current param declaration
                    pval.append([float(lp.strip()) for lp in line[0:-1].split()])
                else:
                    pval.append([float(lp) for lp in l])
    return p

def printCvxPyExpr(e, varname, varprefix="var"):
    n = len(varname)
    s = str(e)
    for j in range(n):
        s = s.replace(varprefix + str(j), varname[j])
    fn = [float(f) for f in floatrx.findall(s)]
    for f in fn:
        sf = str(f)
        srf = str(round(f,3))
        if len(srf) >= len(sf):
            srf = np.format_float_scientific(f, unique=False, precision=3)
        if len(srf) < len(sf):
            s = s.replace(sf, srf)
    print(s)

# print out the flat form of the problem
def flatform(prob, varname, verbose=False):
    #verbose = True
    signs = dict()
    id2con = dict()
    data, chain, inverse_data = prob.get_problem_data(cp.SCS)
    print(len(chain.__dict__['reductions']), "reduction chain (including solver), ", end='')
    print(len(inverse_data), "corresponding inverse_data structures")
    #print("last inverse_data structure in chain (corresponding to solver):")
    #print(inverse_data[4])
    #print("last-but-first inverse_data structure in chain:")
    #pprint(inverse_data[3].__dict__)
    print("first inverse_data structure in chain:")
    #pprint(inverse_data[0].__dict__)
    print(" id_map:", inverse_data[0].id_map)
    print(" var_offsets:", inverse_data[0].var_offsets)
    print(" x_length:", inverse_data[0].x_length)
    print(" var_shapes:", inverse_data[0].var_shapes)
    for j in inverse_data[0].id2var.keys():
        vn = varname[j]
        print(" id2var[{0:d}] := {1}".format(j,vn))
    objn = str(prob.objective)
    for j in inverse_data[0].id2var.keys():
        vn = str(inverse_data[0].id2var[j])
        objn = objn.replace(vn,varname[j])
    print(" obj := {0}".format(objn))
    for h,i in enumerate(inverse_data[0].id2cons.keys()):
        id2con[h] = i
        cl = str(inverse_data[0].id2cons[i].__class__).split('.')        
        signs[i] = cl[3].split('\'')[0]
        cn = str(inverse_data[0].id2cons[i])
        for j in inverse_data[0].id2var.keys():
            vn = str(inverse_data[0].id2var[j])
            cn = cn.replace(vn,varname[j])
        fn = [float(f) for f in floatrx.findall(cn)]
        for f in fn:
            sf = str(f)
            srf = str(round(f,3))
            if len(srf) >= len(sf):
                srf = np.format_float_scientific(f, unique=False, precision=3)
            if len(srf) < len(sf):
                cn = cn.replace(sf, srf)
        print(" constr[{0:d}] := {1}".format(h,cn))
        
    if verbose:
        # this refers to the last reformulation in the chain
        A = data['A']
        b = data['b']
        c = data['c']
        cones = data['dims']
        print("cones:", cones)
        m,n = A.shape
        print("A is {0:d} x {1:d}".format(m,n))
        for i in range(m):
            firstterm = True
            numzeroes = 0
            for j in range(n):
                if abs(A[(i,j)]) > myZero:
                    if firstterm:
                        print("{0:d}: ".format(i+1), end='')
                        firstterm = False
                    else:
                        print(" + ", end='')
                    print("({0:.4f})x{1:d}".format(A[(i,j)], j+1), end='')
                else:
                    numzeroes += 1
            if numzeroes == n:
                print("{0:d}: ".format(i+1), end='')
                print("0.0", end='')
            sn = signs[id2con[i]]
            if sn == 'Inequality':
                print(" <= {0:.4f}".format(b[i]))
            else:
                print(" == {0:.4f}".format(b[i]))

################### MAIN #####################

busRecDescInv = dict()
for i in busRecDesc:
    busRecDescInv[busRecDesc[i]] = i
lineRecDescInv = dict()
for i in lineRecDesc:
    lineRecDescInv[lineRecDesc[i]] = i
genRecDescInv = dict()
for i in genRecDesc:
    genRecDescInv[genRecDesc[i]] = i

t0 = time.time()

## read command line
if len(sys.argv) < 2:
    exit('syntax is: ./acopf_H.py casefile.dat')

## read instance
with open(sys.argv[1], "r") as fd:
    p = readDat(fd)

## parse instance data
buses = {busRecDescInv[i]:[] for i in range(len(busRecDesc))}
bset = []
lines = {lineRecDescInv[i]:[] for i in range(len(lineRecDesc))}
lset = []
gens = {genRecDescInv[i]:[] for i in range(len(genRecDesc))}
gset = []
G = dict()
C = dict()
maxParBranches = 1
Kcard = 2
for (pt,pn,ps,pv) in p:
    if pt == 'P':
        # parameters
        if len(pn) == 1:
            # maxParBranches, Kcard
            if pn[0] == 'maxParBranches':
                maxParBranches = int(pv[0][0])
            elif pn[0] == 'Kcard':
                Kcard = int(pv[0][0])
            elif pn[0] == 'C':
                # generator costs
                for vlist in pv:
                    vl = len(vlist)
                    thebus = int(vlist[0]-offset)
                    theidx = int(vlist[1]-offset)
                    thedeg = int(vlist[2])
                    thecoeff = vlist[3]
                    C[(thebus,theidx,thedeg)] = thecoeff
        elif len(pn) == len(busRecDesc) and pn[0] == busRecDescInv[0]:
            # buses
            for vlist in pv:
                for i,v in enumerate(vlist):
                    if i == 0:
                        bset.append(int(v)-offset)
                    else:
                        buses[busRecDescInv[i-1]].append(v)
        elif len(pn) == len(lineRecDesc) and pn[0] == lineRecDescInv[0]:
            # lines
            for vlist in pv:
                for i,v in enumerate(vlist):
                    if i == 0:
                        lfrom = int(v)-offset
                    elif i == 1:
                        lto = int(v)-offset
                    elif i == 2:
                        lparidx = int(v)-offset
                        lset.append((lfrom,lto,lparidx))
                    else:
                        lines[lineRecDescInv[i-3]].append(v)
        elif len(pn) == len(genRecDesc) and pn[0] == genRecDescInv[0]:
            # generators
            for vlist in pv:
                for i,v in enumerate(vlist):
                    if i == 0:
                        thebus = int(v)-offset
                    elif i == 1:
                        theidx = int(v)-offset
                        gset.append((thebus,theidx))
                    else:
                        gens[genRecDescInv[i-2]].append(v)
                        
    elif pt == 'S':
        # sets
        if pn[0] == 'G':
            Gidx = ps[0]-offset
            G[Gidx] = [vv-offset for vv in pv]

            
# bus quantities
busType = {b:buses['busType'][i] for i,b in enumerate(bset)}
SDR = {b:buses['SDR'][i] for i,b in enumerate(bset)}
SDC = {b:buses['SDC'][i] for i,b in enumerate(bset)}
VL = {b:buses['VL'][i] for i,b in enumerate(bset)}
VU = {b:buses['VU'][i] for i,b in enumerate(bset)}
shR = {b:buses['shR'][i] for i,b in enumerate(bset)}
shC = {b:buses['shC'][i] for i,b in enumerate(bset)}

# line quantities
status = {bah:lines['status'][i] for i,bah in enumerate(lset)}
SU = {bah:lines['SU'][i] for i,bah in enumerate(lset)}
r = {bah:lines['r'][i] for i,bah in enumerate(lset)}
x = {bah:lines['x'][i] for i,bah in enumerate(lset)}
bb = {bah:lines['bb'][i] for i,bah in enumerate(lset)}
tau = {bah:lines['tau'][i] for i,bah in enumerate(lset)}
theta = {bah:lines['theta'][i] for i,bah in enumerate(lset)}
pdLB = {bah:lines['pdLB'][i] for i,bah in enumerate(lset)}
pdUB = {bah:lines['pdUB'][i] for i,bah in enumerate(lset)}

lset1 = [(a,b,h) for (b,a,h) in lset]

# generator quantities
SLR = {bg:gens['SLR'][i] for i,bg in enumerate(gset)}
SLC = {bg:gens['SLC'][i] for i,bg in enumerate(gset)}
SUR = {bg:gens['SUR'][i] for i,bg in enumerate(gset)}
SUC = {bg:gens['SUC'][i] for i,bg in enumerate(gset)}

### MP formulation

# sizes
n = len(bset)  # number of buses (B)
m = len(lset)  # number of lines (L0)
gn = len(gset)  # number of generators

## parameters

# the Y matrix for arcs in L0
Yff = {bah:(1/(r[bah]+1j*x[bah]) + 1j*bb[bah]/2)/tau[bah]**2 for bah in lset}
Yft = {bah:-1/((r[bah]+1j*x[bah])*tau[bah]*cmath.exp(-1j*theta[bah])) for bah in lset}
Ytf = {bah:-1/((r[bah]+1j*x[bah])*tau[bah]*cmath.exp(1j*theta[bah])) for bah in lset}
Ytt = {bah: 1/(r[bah]+1j*x[bah]) + 1j*bb[bah]/2 for bah in lset}

# complex power demand
SD = {b:SDR[b] + 1j*SDC[b] for b in bset}

# line shunt
A = {b:shR[b] + 1j*shC[b] for b in bset}

# in case costs are empty
for (b,g) in gset:
    if (b,g,1) not in C:
        C[(b,g,1)] = 0.0
    if (b,g,2) not in C:
        C[(b,g,2)] = 0.0

## decision variables
varName = dict()

# voltage
varName[0] = "V"
V = cp.Variable(n, complex=True, name=varName[0])

# generated power
varName[1] = "S"
S = cp.Variable(gn, complex=True, name=varName[1])

# Schur(X,VV*)
varName[2] = "X"
X = cp.Variable((n+1,n+1), hermitian=True, name=varName[2])

## variable indexing
bg2g = {bg:g for g,bg in enumerate(gset)}

## objective: minimize generation cost
obj = cp.Minimize(sum(np.real(C[b,g,2])*cp.real(S[d])**2 + np.real(C[b,g,1])*cp.real(S[d]) for d,(b,g) in enumerate(gset))) # + traceCoeff*cp.real(cp.trace(X)))

## constraints
constr = []

# # bounds on voltage
# voltageLR = [ cp.real(V[b]) >= -VU[b] for b in bset]
# constr.extend(voltageLR)
# voltageLC = [ cp.imag(V[b]) >= -VU[b] for b in bset]
# constr.extend(voltageLC)
# voltageUR = [ cp.real(V[b]) <= VU[b] for b in bset]
# constr.extend(voltageUR)
# voltageUC = [ cp.imag(V[b]) <= VU[b] for b in bset]
# constr.extend(voltageUC)

# bounds on generated power
genpowerLR = [ cp.real(S[g]) >= SLR[bg] for g,bg in enumerate(gset) ]
constr.extend(genpowerLR)
genpowerLC = [ cp.imag(S[g]) >= SLC[bg] for g,bg in enumerate(gset) ]
constr.extend(genpowerLC)
genpowerUR = [ cp.real(S[g]) <= SUR[bg] for g,bg in enumerate(gset) ]
constr.extend(genpowerUR)
genpowerUC = [ cp.imag(S[g]) <= SUC[bg] for g,bg in enumerate(gset) ]
constr.extend(genpowerUC)

# # bounds on line current
# #  (divide SU by VL to get relaxation, by VU to get restriction)
linecurrU0 = [ np.real(abs(Yff[(b,a,h)])**2) * cp.real(X[b+s,b+s]) + np.real(abs(Yft[(b,a,h)])**2) * cp.real(X[a+s,a+s]) + 2*cp.real(Yff[(b,a,h)]*np.conj(Yft[(b,a,h)])*X[b+s,a+s]) <= (SU[b,a,h]/VL[b])**2 for (b,a,h) in lset if (SU[b,a,h]>0 and SU[b,a,h]<myLarge) ]
# for c in linecurrU0:
#     e = c.expr
#     printCvxPyExpr(e, varName)
constr.extend(linecurrU0)
linecurrU1 = [ np.real(abs(Ytf[(a,b,h)])**2) * cp.real(X[a+s,a+s]) + np.real(abs(Ytt[(a,b,h)])**2) * cp.real(X[b+s,b+s]) + 2*cp.real(Ytf[(a,b,h)]*np.conj(Ytt[(a,b,h)])*X[a+s,b+s]) <= (SU[a,b,h]/VL[b])**2 for (b,a,h) in lset1 if (SU[a,b,h]>0 and SU[a,b,h]<myLarge) ]
constr.extend(linecurrU1)

# bounds on phase difference
for (b,a,h) in lset:
    if h == 1-offset:
        tanLB = math.tan(pdLB[(b,a,1-offset)])
        tanUB = math.tan(pdUB[(b,a,1-offset)])
        if tanLB > -myLarge:
            phasediffL = [ math.tan(pdLB[(b,a,1-offset)]) * cp.real(X[b+s,a+s]) <= cp.imag(X[b+s,a+s]) ]
            constr.extend(phasediffL)
        if tanUB < myLarge:
            phasediffU = [ math.tan(pdUB[(b,a,1-offset)]) * cp.real(X[b+s,a+s]) >= cp.imag(X[b+s,a+s]) ]
            constr.extend(phasediffU)
phasediffCon = [ cp.real(X[b+s,a+s]) >= 0 for (b,a,h) in lset if h == 1-offset ]
constr.extend(phasediffCon)

# voltage magnitude bounds
voltL = [ VL[b]**2 <= cp.real(X[b+s,b+s]) for b in bset ]
constr.extend(voltL)
voltU = [ cp.real(X[b+s,b+s]) <= VU[b]**2 for b in bset ]
constr.extend(voltU)

# reference bus
rb = -1
for b in bset:
    if busType[b] == 3:
        rb = b
if rb > -1:
    refbus = [ cp.real(V[rb]) >= 0, cp.imag(V[rb]) == 0 ]
    constr.extend(refbus)

# power balance equations
sumL0 = {b:sum(np.conj(Yff[(d,a,h)])*X[d+s,d+s] + np.conj(Yft[(d,a,h)])*X[d+s,a+s] for (d,a,h) in lset if d==b) for b in bset}
sumL1 = {b:sum(np.conj(Ytf[(a,d,h)])*X[d+s,a+s] + np.conj(Ytt[(a,d,h)])*X[d+s,d+s] for (d,a,h) in lset1 if d==b) for b in bset}
sumS = {b:sum(S[bg2g[(d,g)]] for (d,g) in gset if d==b) for b in bset}
power = [ sumL0[b] + sumL1[b] + SD[b] == -cp.conj(A[b])*X[b+s,b+s] + sumS[b] for b in bset ]
constr.extend(power)

# relations between V and X in the Schur complement
VXtop = [ X[0,b+s] == cp.conj(V[b]) for b in bset ]
constr.extend(VXtop)
VXleft = [ X[b+s,0] == V[b] for b in bset ]
constr.extend(VXleft)
## LEO191211: adding constraints below yields lower obj (wrong behaviour for a relaxation, something's wrong)
#VXdiagR = [ cp.real(X[b+s,b+s]) >= cp.real(V[b])**2 + cp.imag(V[b])**2 for b in bset ]
#constr.extend(VXdiagR)
VXdiagC = [ cp.imag(X[b+s,b+s]) == 0.0 for b in bset ]
constr.extend(VXdiagC)

# bottom right corner of Schur complement is identity
VXtopleft = [ X[0,0] == 1 ]
constr.extend(VXtopleft)

# conic constraint "X is Hermitian"
constr.extend([X >> 0])

##################################################
## debugging: fixing variable to optimal values
##   instance being tested: instances/oneline-simple.m
# constr.extend([cp.real(S[0]) <= np.real(1.9998+0j + myCEps), cp.real(S[0]) >= np.real(1.9998+0j - myCEps)])
# constr.extend([cp.imag(S[0]) <= np.imag(1.9998+0j + myCEps), cp.imag(S[0]) >= np.imag(1.9998+0j - myCEps)])
# constr.extend([cp.real(V[0]) <= np.real(2.0+0j + myCEps), cp.real(V[0]) >= np.real(2.0+0j - myCEps)])
# constr.extend([cp.imag(V[0]) <= np.imag(2.0+0j + myCEps), cp.imag(V[0]) >= np.imag(2.0+0j - myCEps)])
# constr.extend([cp.real(V[1]) <= cp.real(1.1+0j + myCEps), cp.real(V[1]) >= cp.real(1.1+0j - myCEps)])
# constr.extend([cp.imag(V[1]) <= cp.imag(1.1+0j + myCEps), cp.imag(V[1]) >= cp.imag(1.1+0j - myCEps)])
# constr.extend([cp.real(X[1,1]) <= 4 + myEps, cp.real(X[1,1]) >= 4 - myEps])
# constr.extend([cp.real(X[2,2]) <= 1.0002 + myEps, cp.real(X[2,2]) >= 1.0002 - myEps])
# constr.extend([cp.real(X[1,2]) <= cp.real(2.0002+0j + myCEps), cp.real(X[1,2]) >= cp.real(2.0002+0j - myCEps)])
# constr.extend([cp.imag(X[1,2]) <= cp.imag(2.0002+0j + myCEps), cp.imag(X[1,2]) >= cp.imag(2.0002+0j - myCEps)])
# constr.extend([cp.real(X[2,1]) <= cp.real(2.0002+0j + myCEps), cp.real(X[2,1]) >= cp.real(2.0002+0j - myCEps)])
# constr.extend([cp.imag(X[2,1]) <= cp.imag(2.0002+0j + myCEps), cp.imag(X[2,1]) >= cp.imag(2.0002+0j - myCEps)])
##################################################


prob = cp.Problem(obj, constr)

## print out the instance formulation ("flat form")
#flatform(prob, varName)

## solve the problem
prob.solve(solver=cp.SCS, verbose=True)
#prob.solve(solver=cp.ECOS, verbose=True)
#prob.solve(solver=cp.MOSEK, verbose=True)
#prob.solve(solver=cp.CVXOPT, verbose=True)

for b in bset:
    if abs(V.value[b]) > myZero:
        print("V[{0:d}] = {1:.3f}".format(b, V.value[b]), end=', ')
    else:
        print("|V[{0:d}]| < {1}".format(b, myZero), end=', ')
    if abs(X.value[b+s,b+s]) > myZero:
        print("X[{0:d},{0:d}] = {1:.3f}".format(b+s, np.real(X.value[b+s,b+s])))
    else:
        print("|X[{0:d},{0:d}]| < {1}".format(b+s, myZero))
for g,(b,d) in enumerate(gset):
    if abs(S.value[g]) > myZero:
        print("S[{0:d}] = {1}".format(g, S.value[g]))
    else:
        print("|S[{0:d}]| < {1}".format(g, myZero))
        
trX = np.real(np.trace(X.value))
rkX = np.linalg.matrix_rank(X.value)
print("trace(X) = {0:.4f}, rank(X) = {1:d}".format(trX,rkX))

objfunval = prob.value # - traceCoeff*trX
print("optimal obj. fun. value =", objfunval)

with open("acopf_H.out", "w") as f:
    print("param : csp ssp := ", file=f)
    for b in bset:
        print(" {0:d} {1:d}  {2} {3}".format(b+s, b+s, np.real(X.value[b+s,b+s]), np.imag(X.value[b+s,b+s])), file=f)
    for (b,a,h) in lset:
        if h == 0:
            print(" {0:d} {1:d}  {2} {3}".format(b+s, a+s, np.real(X.value[b+s,a+s]), np.imag(X.value[b+s,a+s])), file=f)
    print(";", file=f)

    
