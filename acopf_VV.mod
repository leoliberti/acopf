## ACOPF formulation
## acopf-V.mod AMPL file
## Leo Liberti
## This AMPL model implements a formulation described in
##   ../../papers/acopf/formulations/4or/acopf_survey.pdf 
##    it can integrate shunts, transformers or both
##   Formulation in real numbers (rather than complex),
##     decision variables: generated power, voltage
## NOTE: this formulation uses injected current bounds (either rel or restr)
## 191031  derived from acopf_SIV.mod through replacement I=YV, S=VI^*
##         => S = V(YV)^* = V Y^* V^* = Y^* VV^*,
##         expressed in function of vars VbR,VbC,|Vb|^2,(VbVa^*)^r,(VbVa^*)^c
## 191112  note: since we don't have upper current bounds on lines,
##         this is either a relaxation (SU/VL) or a restriction (SU/VU)
##         usually, get same answer as acopf_SIV or very similar

let formulation_type := "cartesian";

## common parameters read from .run

## local parameters

# various products of Y components
param Yff2{(b,a,i) in L0} default YffR[b,a,i]^2 + YffC[b,a,i]^2;
param Yft2{(b,a,i) in L0} default YftR[b,a,i]^2 + YftC[b,a,i]^2;
param Ytf2{(b,a,i) in L0} default YtfR[b,a,i]^2 + YtfC[b,a,i]^2;
param Ytt2{(b,a,i) in L0} default YttR[b,a,i]^2 + YttC[b,a,i]^2;
param YffYftR{(b,a,i) in L0} default YffR[b,a,i]*YftR[b,a,i] + YffC[b,a,i]*YftC[b,a,i];
param YffYftC{(b,a,i) in L0} default YffR[b,a,i]*YftC[b,a,i] - YffC[b,a,i]*YftR[b,a,i];
param YtfYttR{(b,a,i) in L0} default YtfR[b,a,i]*YttR[b,a,i] + YtfC[b,a,i]*YttC[b,a,i];
param YtfYttC{(b,a,i) in L0} default YtfR[b,a,i]*YttC[b,a,i] - YtfC[b,a,i]*YttR[b,a,i];

## common decision variables read from .run

## local decision variables

# V2 = |V|^2
var V2{b in B} >= VL[b]^2, <= VU[b]^2;
# VbaR = VR[b]*VR[a]+VC[b]*VC[a]
var VbaR{b in B, a in B : (b,a,1) in L};
# VbaC = VC[b]*VR[a]-VR[b]*VC[a]
var VbaC{b in B, a in B : (b,a,1) in L};

## objective function

# generation cost (WARNING: to remove quadratic power terms set C[g,2]=0)
minimize gencost_VV: sum{b in B, g in G[b]} (C[b,g,2]*SgenR[b,g]^2 + C[b,g,1]*SgenR[b,g] + C[b,g,0]);

## constraints

## power flow, derived from replacing S_ba with V_b(Y_bb V_b - Y_ba V_a)^*

# #################################################################

## this version is computed by hand (last check 191113)
subject to powerflowR_VV{b in B}:  # real
   SDR[b] +
   sum{(b,a,i) in L0}
     (YffR[b,a,i]*V2[b] + YftR[b,a,i]*VbaR[a,b] - YftC[b,a,i]*VbaC[a,b]) +
   sum{(b,a,i) in L1}
     (YttR[a,b,i]*V2[b] + YtfR[a,b,i]*VbaR[b,a] + YtfC[a,b,i]*VbaC[b,a]) =
   -shR[b]*V2[b] + sum{g in G[b]} SgenR[b,g];
subject to powerflowC_VV{b in B}:  # imaginary
   SDC[b] +
   sum{(b,a,i) in L0}
     (-YffC[b,a,i]*V2[b] + YftR[b,a,i]*VbaC[a,b] - YftC[b,a,i]*VbaR[a,b]) +
   sum{(b,a,i) in L1}
     (-YttC[a,b,i]*V2[b] - YtfR[a,b,i]*VbaC[b,a] - YtfC[a,b,i]*VbaR[b,a]) =
    shC[b]*V2[b] + sum{g in G[b]} SgenC[b,g];

## power bound on lines b->a defined on |S| <= SU:
##   reformulate to current bounds:
##   SU=>|S|=|V||I|>=VL|I| --> |I|<=SU/VL for relaxation
##   |I|<=SU/VU for restriction
subject to powerbound1_VV{(b,a,i) in L0 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
 Yff2[b,a,i]*V2[b] + Yft2[b,a,i]*V2[a] +
   2*(YffYftR[b,a,i]*VbaR[b,a] -
      YffYftC[b,a,i]*VbaC[b,a]) <= (SU[b,a,i]/VL[b])^2;
subject to powerbound2_VV{(b,a,i) in L1 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
 Ytf2[a,b,i]*V2[a] + Ytt2[a,b,i]*V2[b] +
   2*(YtfYttR[a,b,i]*VbaR[a,b] -
      YtfYttC[a,b,i]*VbaC[a,b]) <= (SU[b,a,i]/VL[b])^2;
## above constr indexed on L1 is equivalent to below indexed on L0:
#subject to powerbound2_VV{(b,a,i) in L0 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
# Ytf2[b,a,i]*V2[b] + Ytt2[b,a,i]*V2[a] +
#   2*(YtfYttR[b,a,i]*VbaR[b,a] -
#      YtfYttC[b,a,i]*VbaC[b,a]) <= (SU[b,a,i]/VU[b])^2;

# definition of V2
subject to V2def_VV{b in B}: V2[b] = VR[b]^2 + VC[b]^2;
# definition of VbaR (re(VbVa*))
subject to VbaRdef_VV{b in B, a in B : (b,a,1) in L}:
  VbaR[b,a] = VR[b]*VR[a] + VC[b]*VC[a];
# definition of VbaC (im(VbVa*))
subject to VbaCdef_VV{b in B, a in B : (b,a,1) in L}:
  VbaC[b,a] = VC[b]*VR[a] - VR[b]*VC[a];

# bounds on phase difference
subject to phasediff1_VV{(b,a,i) in L0}: VbaC[b,a]<=tan(pdUB[b,a,i])*VbaR[b,a];
subject to phasediff2_VV{(b,a,i) in L0}: VbaC[b,a]>=tan(pdLB[b,a,i])*VbaR[b,a];
subject to phasediff3_VV{(b,a,1) in L0}: VbaR[b,a] >= 0;

# reference bus: there had better be just one reference -- check in .run
subject to reference1_VV{b in B : busType[b] == 3}: VC[b] = 0;
subject to reference2_VV{b in B : busType[b] == 3}: VR[b] >= 0;

