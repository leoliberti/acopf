#!/usr/bin/env python

## mpc2dat.py
## Leo Liberti
## read a matpower case .m file and output the equivalent .dat file for acopf.mod
## 180907 work started
## 190930 major overhaul for new formulation after Dan's suggestions
## 191013 finally working for .dat output (matrix output still to be corrected)

################ imports ################

import sys
import os
import types
import math
import cmath
import numpy as np

################ user-configurable ###############

myInf = 1e30
myZero = 1e-10
verbose = False
#verbose = True

################ Classes ################

class mpcBus:
    def __init__(self, col):
        self.ID = col[0]
        self.Type = int(col[1])
        self.Pd = float(col[2])
        self.Qd = float(col[3])
        self.Gs = float(col[4])
        self.Bs = float(col[5])
        self.area = int(col[6])
        self.Vm = float(col[7]) # starting point for polar formulation on input
        self.Va = (math.pi) * float(col[8]) / 180 # startpt for polar on input
        self.baseKV = float(col[9])
        self.zone = int(col[10])
        self.Vmax = float(col[11])
        self.Vmin = float(col[12])
    def fieldnames(self):
        out = "ID,Type,Pd,Qd,Gs,Bs,area,Vm,Va,baseKV,zone,Vmax,Vmin"
        return out
    # printing
    def __str__(self):
        out = "bus " + str(self.ID) + ","
        out += str(self.Type) + "," + str(self.Pd) + "," + str(self.Qd) + "," + str(self.Gs) + ","
        out += str(self.Bs) + "," + str(self.area) + "," + str(self.Vm) + "," + str(self.Va) + ","
        out += str(self.baseKV) + "," + str(self.zone) + "," + str(self.Vmax) + "," + str(self.Vmin)
        return out    

class mpcGen:
    # gen and gencost dicts are indexed by same keys
    def __init__(self, col, count):
        self.bus = col[0]        # bus where generator is installed
        self.counter = count     # counter of generator at bus
        self.Pg = float(col[1])  # Pg is a decision var
        self.Qg = float(col[2])  # Qg is a decision var
        self.Qmax = float(col[3])
        self.Qmin = float(col[4])
        self.Vg = float(col[5])
        self.mBase = float(col[6])
        self.status = int(col[7])
        self.Pmax = float(col[8])
        self.Pmin = float(col[9])
        try:
            self.Pc1 = float(col[10])
            self.Pc2 = float(col[11])
            self.Qc1min = float(col[12])
            self.Qc1max = float(col[13])
            self.Qc2min = float(col[14])
            self.Qc2max = float(col[15])
            self.ramp_agc = float(col[16])
            self.ramp_10 = float(col[17])
            self.ramp_30 = float(col[18])
            self.ramp_q = float(col[19])
            self.apf = float(col[20])
        except:
            pass
    def fieldnames(self):
        out = "ID,Pg,Qg,Qmax,Qmin,Vg,mBase,status,Pmax,Pmin,Pc1,Pc2,Qc1min,Qc1max,Qc2min,Qc2max,ramp_agc,ramp_10,ramp_30,ramp_q,apf"
        return out
    # printing
    def __str__(self):
        out = "gen " + str(self.bus) + "," + str(self.counter) + ":"
        out += str(self.Pg) + "," + str(self.Qg) + "," + str(self.Qmax) + "," + str(self.Qmin) + ","
        out += str(self.Vg) + "," + str(self.mBase) + "," + str(self.status) + "," + str(self.Pmax) + ","
        out += str(self.Pmin) + ","
        try:
            out += str(self.Pc1) + "," + str(self.Pc2) + "," + str(self.Qc1min) + ","
            out += str(self.Qc1max) + "," + str(self.Qc2min) + "," + str(self.Qc2max) + "," + str(self.ramp_agc) + ","
            out += str(self.ramp_10) + "," + str(self.ramp_30) + "," + str(self.ramp_q) + "," + str(self.apf)
        except:
            pass
        return out
    
class mpcBranch:
    def __init__(self, col, parid):
        self.fbus = col[0]
        self.tbus = col[1]
        # ID of parallel edge
        self.parallelID = parid
        self.r = float(col[2])
        self.x = float(col[3])
        self.b = float(col[4])
        self.rateA = float(col[5])
        self.rateB = float(col[6])
        self.rateC = float(col[7])
        self.ratio = float(col[8])
        self.angle = float(col[9])
        self.status = int(col[10])
        self.angmin = float(col[11])
        self.angmax = float(col[12])
    def fieldnames(self):
        out = "fbus,tbus,r,x,b,rateA,rateB,rateC,ratio,angle,status,angmin,angmax"
        return out
    # printing
    def __str__(self):
        out = "branch " + str(self.fbus) + "," + str(self.tbus) + "," + str(self.parallelID) + ":"
        out += str(self.r) + "," + str(self.x) + "," + str(self.b) + "," + str(self.rateA) + ","
        out += str(self.rateB) + "," + str(self.rateC) + "," + str(self.ratio) + "," + str(self.angle) + ","
        out += str(self.status) + "," + str(self.angmin) + "," + str(self.angmax)
        return out

class mpcGenCost:
    # gen and gencost dicts are indexed by same keys
    def __init__(self, col):
        self.Type = int(col[0])
        self.startup = float(col[1])
        self.shutdown = float(col[2])
        self.n = int(col[3])
        self.data = [float(c) for c in col[4:]]
    def fieldnames(self):
        out = "type,startup,shutdown,n"
        for i in range(self.n):
            out += ",data" + str(i)
        return out
    # printing
    def __str__(self):
        out = "gencost " + str(self.Type) + "," + str(self.startup) + "," + str(self.shutdown) + "," + str(self.n)
        for i in range(self.n):
            out += "," + str(self.data[i])
        return out
    
class mpcCase:
    def __init__(self, mpcfilename):
        self.filename = mpcfilename
        f = open(self.filename, "r")
        lines = f.readlines()
        status = "outer"
        self.bus = dict()
        self.branch = dict()
        self.parbranch = dict()
        self.maxparbranches = 0
        self.gen = dict()
        self.gencost = dict()
        buscounter = 0
        branchcounter = 0
        gencounter = 0
        gencostcounter = 0
        for l in lines:
            line = l.strip(' \n').replace('\t', ' ')
            if len(line)>0 and line[0] != '%':
                if status == "outer":
                    if line.startswith("mpc.version"):
                        self.version = int(line.split('=')[1].strip(' ;\r').strip('\''))
                    elif line.startswith("mpc.baseMVA"):
                        self.baseMVA = float(line.split('=')[1].strip(' ;\r'))
                    elif line.startswith("mpc.bus"):
                        status = "bus"
                    elif line.startswith("mpc.gen") and not line.startswith("mpc.gencost"):
                        status = "gen"
                    elif line.startswith("mpc.branch"):
                        status = "branch"
                    elif line.startswith("mpc.gencost"):
                        status = "gencost"
                else:
                    if "]" in line:
                        status = "outer"
                        rows = line.split(']')[0]
                    elif ";" in line:
                        rows = line.split(';')
                    else:
                        rows = [line]
                    for r in rows:
                        if len(r) > 1:
                            cols = r.split()
                            if cols[0] == "%":
                                continue
                            if status == "bus" and len(cols) >= 12:
                                self.bus[buscounter] = mpcBus(cols)
                                if verbose:
                                    print self.bus[buscounter]
                                buscounter += 1
                            elif status == "gen" and len(cols) >= 10:
                                # gen and gencost are indexed by same keys
                                count = len([self.gen[g].bus for g in self.gen if self.gen[g].bus == cols[0]])
                                self.gen[gencounter] = mpcGen(cols, count)
                                if verbose:
                                    print self.gen[gencounter]
                                gencounter += 1 # gencost idx (starts from 1)
                            elif status == "branch" and len(cols) >= 13:
                                ID = (cols[0],cols[1]) # ID is link adjacencies
                                if ID in self.parbranch:
                                    self.parbranch[ID] += 1
                                else:
                                    self.parbranch[ID] = 1
                                self.branch[branchcounter] = mpcBranch(cols, self.parbranch[ID])
                                if self.parbranch[ID] > self.maxparbranches:
                                    self.maxparbranches = self.parbranch[ID]
                                if verbose:
                                    print self.branch[branchcounter]
                                branchcounter += 1
                            elif status == "gencost" and len(cols) >= 4:
                                # gen and gencost are indexed by same keys
                                self.gencost[gencostcounter] = mpcGenCost(cols)
                                if verbose:
                                    print self.gencost[gencostcounter]
                                gencostcounter += 1 
        f.close()

    def fieldnames(self):
        out = "CASEFN:filename,version,baseMVA"
        return out

    # printing
    def __str__(self,fnflag=True):
        if fnflag:
            out = self.fieldnames() + "\n"
        else:
            out = ""
        out += self.filename + ":" + str(self.version) + "," + str(self.baseMVA)
        if fnflag:
            out += "\n" + "BUSFN:" + self.bus.values()[0].fieldnames()
        for b in self.bus.values():
            out += "\n" + str(b)
        if fnflag:
            out += "\n" + "GENFN:" + self.gen.values()[0].fieldnames()
        for g in self.gen.values():
            out += "\n" + str(g)
        if fnflag:
            out += "\n" + "BRANCHFN:" + self.branch.values()[0].fieldnames()
        for l in self.branch.values():
            out += "\n" + str(l)
        if fnflag:
            out += "\n" + "GENCOSTFN:" + self.gencost.values()[0].fieldnames()
        for gc in self.gencost.values():
            out += "\n" + str(gc)
        return out

    def writeAMPLdat(self,filename):
        f = open(filename, "w")
        print >> f, "# AMPL .dat file translation of", self.filename
        print >> f, "param maxParBranches :=", self.maxparbranches, ";"
        ## set L0
        #print >> f, "\nset L0 :=",
        #for l in self.branch:
        #    FROM = self.branch[l].fbus
        #    TO = self.branch[l].tbus
        #    PARID = self.branch[l].parallelID
        #    print >> f, "({0:d},{1:d},{2:d})".format(FROM,TO,PARID),
        #print >> f, ";\n"
        # print data about buses
        print >> f, "\nparam : B : busType SDR SDC VL VU Vm Va shR shC :="
        for b in self.bus:
            # power demand at bus (any bus type)
            ID = self.bus[b].ID
            SDR = self.bus[b].Pd/self.baseMVA # from matpower's opf_setup.m
            SDC = self.bus[b].Qd/self.baseMVA # from matpower's opf_setup.m
            shR = self.bus[b].Gs/self.baseMVA # from matpower's makeYbus.m
            shC = self.bus[b].Bs/self.baseMVA # from matpower's makeYbus.m
            Vm = self.bus[b].Vm
            Va = self.bus[b].Va
            VL = self.bus[b].Vmin # see lib/opf_vlim_fcn.m in matpower
            VU = self.bus[b].Vmax # see lib/opf_vlim_fcn.m in matpower
            print >> f, " ", ID, " ", self.bus[b].Type, SDR, SDC, VL, VU, Vm, Va, shR, shC
        # print data about generators
        for g in self.gen:
            if self.gen[g].counter == 0:
                print >> f, ";"
                print >> f, "set G[{0:s}] := 1".format(self.gen[g].bus),
            else:
                print >> f, self.gen[g].counter + 1,
        print >> f, ";\nparam : SLR SLC SUR SUC :="
        for g in self.gen:
            # status=0 => inactive (from Escobar's email 191003)
            status = self.gen[g].status
            # matpower's opf_setup.m
            SLR = status * self.gen[g].Pmin / self.baseMVA 
            SLC = status * self.gen[g].Qmin / self.baseMVA 
            SUR = status * self.gen[g].Pmax / self.baseMVA 
            SUC = status * self.gen[g].Qmax / self.baseMVA
            print >> f, " ", self.gen[g].bus, self.gen[g].counter + 1, " ", SLR, SLC, SUR, SUC;
        print >> f, ";"
        # print data about links
        print >> f, "\nparam : L0 : status SU r x bb tau nu pdLB pdUB :="
        for l in self.branch:
            ## WARNING: 'l' here is the letter 'ell' not the digit 'one'
            FROM = self.branch[l].fbus
            TO = self.branch[l].tbus   ## WARNING: this is 'ell at index one'
            PARID = self.branch[l].parallelID
            status = self.branch[l].status
            if self.branch[l].rateA > 0:
                # semantics of current upper bound: only if >0 (qcqp_opf.m)
                # scaling by baseMVA taken from qcqp_opf.m, squaring in .mod,
                #   see (opf_branch_flow_fcn.m in matpower)
                SU = self.branch[l].rateA / self.baseMVA 
            else:
                SU = myInf
            r = self.branch[l].r
            x = self.branch[l].x
            # bb called "line charging susceptance" in MatPower:makeYbus.m
            bb = self.branch[l].b 
            tap = self.branch[l].ratio # a.k.a. "TAP"
            if abs(tap) < myZero:
                tap = 1.0;
            ang = self.branch[l].angle # a.k.a. "SHIFT"
            angmin = self.branch[l].angmin # phasediff across line ends LB
            angmax = self.branch[l].angmax # phasediff across line ends UB
            if angmin == 0 or angmin == -360:
                angmin = -90  # since we're going to apply tan() later
            if angmax == 0 or angmax == 360:
                angmax = 90   # since we're going to apply tan() later
            # transform to radians
            ang = (math.pi) * ang / 180
            angmin = (math.pi) * angmin / 180 
            angmax = (math.pi) * angmax / 180 
            print >> f, " ", FROM, TO, PARID, " ", status, SU, r, x, bb, tap, ang, angmin, angmax
        print >> f, ";"
        if len(self.gencost) > 0:
            # print coefficients of cost function
            #   this assumes that gencost data all have the same length
            print >> f, "\nparam Kcard :=", self.gencost[0].n, ";"
            print >> f, "\nparam C :="
            for g in self.gencost:
                K = range(self.gencost[g].n)
                for k in K:
                    ## degree after Escobar's acopf.m code
                    degree = self.gencost[g].n-k-1
                    coeffval = self.gen[g].status * self.gencost[g].data[k] * (self.baseMVA**degree)
                    if coeffval != 0:
                        print >> f, " ", int(self.gen[g].bus), int(self.gen[g].counter)+1, degree, " ", coeffval
            print >> f, ";"
        f.close()

    def writeMatrices(self):
        ## notation as per Zimmermann, with the exception of transformers from Bienstock
        nB = len(self.bus)
        busorder = self.bus.keys()
        nL = len(self.branch)
        lineorder = self.branch.keys()
        # compute bus shunt elements
        YshL = []
        for b in self.bus:
            yshb = self.bus[b].Gs + 1j*self.bus[b].Bs
            YshL.append(yshb) # weird bug <github.com/numpy/numpy/issues/4379> forces use of list first
        Ysh = np.array(YshL) / self.baseMVA # scaling from makeYbus.m
        # compute Cf,Ct (directed incidences line/bus)
        Cf = np.zeros((nL,nB))
        Ct = np.zeros((nL,nB))
        for lidx,line in enumerate(self.branch):
            b = line[0]
            a = line[1]
            bidx = busorder.index(b)
            aidx = busorder.index(a)
            Cf[lidx,bidx] = 1
            Ct[lidx,aidx] = 1
        # compute Ybr matrices
        Ybr = dict()
        for line in self.branch:
            zs = self.branch[line].r + 1j*self.branch[line].x
            ys = 1/zs
            ylinesh = 0 + 1j*self.branch[line].b
            bc = ylinesh.imag
            tau = self.branch[line].ratio
            sigma = self.branch[line].angle
            N = tau*cmath.exp(1j*sigma)
            if abs(N) < myZero:
                tau = 1
            ytt = (ys + 1j*bc/2) / (tau**2)
            yft = -ys / (tau*cmath.exp(-1j*sigma))
            ytf = -ys / (tau*cmath.exp(1j*sigma))
            yff = ys + 1j*bc/2
            Ybr[line] = np.array([[yff,yft],[yft,yff]])
        # compute Y vectors
        YffL = []
        YftL = []
        YtfL = []
        YttL = []
        for line in self.branch:
            YffL.append(Ybr[line][0,0])
            YftL.append(Ybr[line][0,1])
            YtfL.append(Ybr[line][1,0])
            YttL.append(Ybr[line][1,1])
        Yff = np.array(YffL)
        Yft = np.array(YftL)
        Ytf = np.array(YtfL)
        Ytt = np.array(YttL)
        Yf = np.dot(np.diag(Yff),Cf) + np.dot(np.diag(Yft),Ct)
        Yt = np.dot(np.diag(Ytf),Cf) + np.dot(np.diag(Ytt),Ct)
        Ybus = np.dot(np.transpose(Cf),Yf) + np.dot(np.transpose(Ct),Yt) + np.diag(Ysh)
        # RTE's qcqp_gen.m has Ybus = baseMVA*Ybus, MatPower doesn't appear to have this
            
        
################## main #####################

if len(sys.argv) < 2:
    print "command line:", sys.argv[0], "casefile.m"
    sys.exit(1)

    
caseFile = sys.argv[1]
print "mpc2dat: translating", caseFile

filename = os.path.basename(caseFile)
dirname = os.path.dirname(caseFile)
basename = ''.join([t + '.' for t in filename.split('.')[0:-1]])[0:-1]

mpc = mpcCase(caseFile)
#print mpc
if len(dirname) > 0:
    datname = dirname + "/" + basename + ".dat"
else:
    datname = basename + ".dat"

# write the .dat file
mpc.writeAMPLdat(datname)
print "  written output to", datname

#mpc.writeMatrices()

