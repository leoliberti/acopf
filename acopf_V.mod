## ACOPF formulation
## acopf_V.mod AMPL file
## Leo Liberti
## This AMPL model implements a formulation described in
##   ../../papers/acopf/formulations/4or/acopf_survey.pdf 
##    it can integrate shunts, transformers or both
##   Formulation in real numbers (rather than complex),
##     decision variables: generated power, voltage
## NOTE: this formulation uses injected current bounds (either rel or restr)
## 191031  derived from acopf_SIV.mod through replacement I=YV, S=VI^*
##         => S = V(YV)^* = V Y^* V^* = Y^* VV^*
##         expressed in function of variables VbR,VbC only
## 191112  note: since we don't have upper current bounds on lines,
##         this is either a relaxation (SU/VL) or a restriction (SU/VU)
##         usually, get same answer as acopf_SIV or very similar
## 191203  fixed some bugs using SymPy output

let formulation_type := "cartesian";

## common parameters and variables read from .run

## objective function

# generation cost (WARNING: to remove quadratic power terms set C[g,2]=0)
minimize gencost_V: sum{b in B, g in G[b]} (C[b,g,2]*SgenR[b,g]^2 + C[b,g,1]*SgenR[b,g] + C[b,g,0]);

## constraints

## power flow, derived by SymPy from repl S_ba with V_b(Y_bb V_b - Y_ba V_a)^*
subject to powerflowR_V{b in B}:  # real
   SDR[b] +
   sum{(b,a,i) in L0}
     (VC[a]*VC[b]*YftR[b,a,i] - VC[a]*VR[b]*YftC[b,a,i] + VC[b]*VR[a]*YftC[b,a,i] + VR[a]*VR[b]*YftR[b,a,i] + YffR[b,a,i]*(VC[b]^2 + VR[b]^2)) +
   sum{(b,a,i) in L1}
     (VC[a]*VC[b]*YtfR[a,b,i] - VC[a]*VR[b]*YtfC[a,b,i] + VC[b]*VR[a]*YtfC[a,b,i] + VR[a]*VR[b]*YtfR[a,b,i] + YttR[a,b,i]*(VC[b]^2 + VR[b]^2)) =
   -shR[b]*(VR[b]^2 + VC[b]^2) + sum{g in G[b]} SgenR[b,g];
subject to powerflowC_V{b in B}:  # imaginary
   SDC[b] +
   sum{(b,a,i) in L0}
     (-VC[a]*VC[b]*YftC[b,a,i] - VC[a]*VR[b]*YftR[b,a,i] + VC[b]*VR[a]*YftR[b,a,i] - VR[a]*VR[b]*YftC[b,a,i] - YffC[b,a,i]*(VC[b]^2 + VR[b]^2)) +
   sum{(b,a,i) in L1}
     (-VC[a]*VC[b]*YtfC[a,b,i] - VC[a]*VR[b]*YtfR[a,b,i] + VC[b]*VR[a]*YtfR[a,b,i] - VR[a]*VR[b]*YtfC[a,b,i] - YttC[a,b,i]*(VC[b]^2 + VR[b]^2)) =
    shC[b]*(VR[b]^2 + VC[b]^2) + sum{g in G[b]} SgenC[b,g];

## power bound on lines b->a defined on |S| <= SU
##   reformulate to current bounds:
##     SU=>|S|=|V||I|>=VL|I| --> |I|<=SU/VL for relaxation
##     |I|<=SU/VU for restriction
subject to powerbound1_V{(b,a,i) in L0 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  ## by SymPy
  VC[a]^2*YftC[b,a,i]^2 + VC[a]^2*YftR[b,a,i]^2 + 2*VC[a]*VC[b]*YffC[b,a,i]*YftC[b,a,i] + 2*VC[a]*VC[b]*YffR[b,a,i]*YftR[b,a,i] + 2*VC[a]*VR[b]*YffC[b,a,i]*YftR[b,a,i] - 2*VC[a]*VR[b]*YffR[b,a,i]*YftC[b,a,i] + VC[b]^2*YffC[b,a,i]^2 + VC[b]^2*YffR[b,a,i]^2 - 2*VC[b]*VR[a]*YffC[b,a,i]*YftR[b,a,i] + 2*VC[b]*VR[a]*YffR[b,a,i]*YftC[b,a,i] + VR[a]^2*YftC[b,a,i]^2 + VR[a]^2*YftR[b,a,i]^2 + 2*VR[a]*VR[b]*YffC[b,a,i]*YftC[b,a,i] + 2*VR[a]*VR[b]*YffR[b,a,i]*YftR[b,a,i] + VR[b]^2*YffC[b,a,i]^2 + VR[b]^2*YffR[b,a,i]^2 <= (SU[b,a,i]/VL[b])^2;
  ## by Leo
  #   (YftR[b,a,i]^2 + YftC[b,a,i]^2)*VC[a]^2 +
  # 2*(YffR[b,a,i]*YftR[b,a,i] + YffC[b,a,i]*YftC[b,a,i])*VC[a]*VC[b] -
  # 2*(YffR[b,a,i]*YftC[b,a,i] - YffC[b,a,i]*YftR[b,a,i])*VC[a]*VR[b] +
  #   (YffR[b,a,i]^2 + YffC[b,a,i]^2)*VC[b]^2 +
  # 2*(YffR[b,a,i]*YftC[b,a,i] - YffC[b,a,i]*YftR[b,a,i])*VR[a]*VC[b] +
  #   (YffR[b,a,i]^2 + YffC[b,a,i]^2)*VR[b]^2 +
  #   (YftR[b,a,i]^2 + YftC[b,a,i]^2)*VR[a]^2 +
  # 2*(YffR[b,a,i]*YftR[b,a,i] + YffC[b,a,i]*YftC[b,a,i])*VR[a]*VR[b]
  #   <= (SU[b,a,i]/VL[b])^2;

subject to powerbound2_V{(b,a,i) in L1 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  ## by SymPy
  VC[a]^2*YtfC[a,b,i]^2 + VC[a]^2*YtfR[a,b,i]^2 + 2*VC[a]*VC[b]*YtfC[a,b,i]*YttC[a,b,i] + 2*VC[a]*VC[b]*YtfR[a,b,i]*YttR[a,b,i] - 2*VC[a]*VR[b]*YtfC[a,b,i]*YttR[a,b,i] + 2*VC[a]*VR[b]*YtfR[a,b,i]*YttC[a,b,i] + VC[b]^2*YttC[a,b,i]^2 + VC[b]^2*YttR[a,b,i]^2 + 2*VC[b]*VR[a]*YtfC[a,b,i]*YttR[a,b,i] - 2*VC[b]*VR[a]*YtfR[a,b,i]*YttC[a,b,i] + VR[a]^2*YtfC[a,b,i]^2 + VR[a]^2*YtfR[a,b,i]^2 + 2*VR[a]*VR[b]*YtfC[a,b,i]*YttC[a,b,i] + 2*VR[a]*VR[b]*YtfR[a,b,i]*YttR[a,b,i] + VR[b]^2*YttC[a,b,i]^2 + VR[b]^2*YttR[a,b,i]^2	<= (SU[b,a,i]/VL[b])^2;
  ## by Leo (for b,a,i in L0)
  #   (YtfC[b,a,i]^2 + YtfR[b,a,i]^2)*VC[a]^2 +
  # 2*(YtfC[b,a,i]*YttC[b,a,i] + YtfR[b,a,i]*YttR[b,a,i])*VC[a]*VC[b] +
  # 2*(YtfR[b,a,i]*YttC[b,a,i] - YtfC[b,a,i]*YttR[b,a,i])*VC[a]*VR[b] +
  #   (YttR[b,a,i]^2 + YttC[b,a,i]^2)*VC[b]^2 -
  # 2*(YtfR[b,a,i]*YttC[b,a,i] - YtfC[b,a,i]*YttR[b,a,i])*VR[a]*VC[b] +
  #   (YtfR[b,a,i]^2 + YtfC[b,a,i]^2)*VR[a]^2 +
  #   (YttR[b,a,i]^2 + YttC[b,a,i]^2)*VR[b]^2 +
  # 2*(YtfR[b,a,i]*YttR[b,a,i] + YtfC[b,a,i]*YttC[b,a,i])*VR[a]*VR[b]
  #   <= (SU[b,a,i]/VL[b])^2;

# bounds on phase difference
subject to phasediff1_V{(b,a,i) in L0}: 
  tan(pdUB[b,a,i])*(VR[b]*VR[a]+VC[b]*VC[a]) + (VR[b]*VC[a]-VC[b]*VR[a]) >= 0;
subject to phasediff2_V{(b,a,i) in L0}:
  tan(pdLB[b,a,i])*(VR[b]*VR[a]+VC[b]*VC[a]) + (VR[b]*VC[a]-VC[b]*VR[a]) <= 0;
subject to phasediff3_V{(b,a,1) in L0}: VR[b]*VR[a]+VC[b]*VC[a] >= 0;

# reference bus: there had better be just one reference -- check in .run
subject to reference1_V{b in B : busType[b] == 3}: VC[b] = 0;
subject to reference2_V{b in B : busType[b] == 3}: VR[b] >= 0;

# bounds on voltage magnitude
subject to vmagbound1_V{b in B}: VR[b]^2 + VC[b]^2 <= VU[b]^2;
subject to vmagbound2_V{b in B}: VR[b]^2 + VC[b]^2 >= VL[b]^2;

