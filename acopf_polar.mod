## ACOPF formulation
## acopf_V_polar.mod AMPL file
## Leo Liberti
## This AMPL model implements a formulation described in
##   ../../papers/acopf/formulations/4or/acopf_survey.pdf 
##    it can integrate shunts, transformers or both
##   Formulation in real numbers (rather than complex),
##     decision variables: generated power, voltage
## NOTE: this formulation uses injected power bounds, so it is "exact"
## 191031  derived from acopf_SIV.mod through replacement
##         I=YV with polar coordinates, then S=VI^*
##         expressed in function of polar variables v_b, theta_b
## 191203  polar formulation derived by SymPy

let formulation_type := "polar";

## common parameters and variables read from .run file

## objective function

# generation cost (WARNING: to remove quadratic power terms set C[g,2]=0)
minimize gencost_pol: sum{b in B, g in G[b]} (C[b,g,2]*SgenR[b,g]^2 + C[b,g,1]*SgenR[b,g] + C[b,g,0]);

### constraints

## power flow balance
subject to powerflowR_pol{b in B}:  # real
   SDR[b] +
   sum{(b,a,i) in L0}
     (YffR[b,a,i]*v[b]^2 + YftC[b,a,i]*v[b]*v[a]*sin(theta[b] - theta[a]) + YftR[b,a,i]*v[b]*v[a]*cos(theta[b] - theta[a])) +
   sum{(b,a,i) in L1}
     (YtfC[a,b,i]*v[b]*v[a]*sin(theta[b] - theta[a]) + YtfR[a,b,i]*v[b]*v[a]*cos(theta[b] - theta[a]) + YttR[a,b,i]*v[b]^2) =
   -shR[b]*v[b]^2 + sum{g in G[b]} SgenR[b,g];
subject to powerflowC_pol{b in B}:  # imaginary
   SDC[b] +
   sum{(b,a,i) in L0}
     (-YffC[b,a,i]*v[b]^2 - YftC[b,a,i]*v[b]*v[a]*cos(theta[b] - theta[a]) + YftR[b,a,i]*v[b]*v[a]*sin(theta[b] - theta[a])) +
   sum{(b,a,i) in L1}
     (-YtfC[a,b,i]*v[b]*v[a]*cos(theta[b] - theta[a]) + YtfR[a,b,i]*v[b]*v[a]*sin(theta[b] - theta[a]) - YttC[a,b,i]*v[b]^2) =
    shC[b]*v[b]^2 + sum{g in G[b]} SgenC[b,g];

## power bound on lines b->a defined on |S|^2 <= SU^2: (computed using sympy)
subject to powerbound1_pol{(b,a,i) in L0 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  YffC[b,a,i]^2*v[b]^4 + 2*YffC[b,a,i]*YftC[b,a,i]*v[a]*v[b]^3*cos(theta[b] - theta[a]) - 2*YffC[b,a,i]*YftR[b,a,i]*v[a]*v[b]^3*sin(theta[b] - theta[a]) + YffR[b,a,i]^2*v[b]^4 + 2*YffR[b,a,i]*YftC[b,a,i]*v[a]*v[b]^3*sin(theta[b] - theta[a]) + 2*YffR[b,a,i]*YftR[b,a,i]*v[a]*v[b]^3*cos(theta[b] - theta[a]) + YftC[b,a,i]^2*v[a]^2*v[b]^2 + YftR[b,a,i]^2*v[a]^2*v[b]^2 <= SU[b,a,i]^2;
subject to powerbound2_pol{(b,a,i) in L1 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  YtfC[a,b,i]^2*v[a]^2*v[b]^2 + 2*YtfC[a,b,i]*YttC[a,b,i]*v[a]*v[b]^3*cos(theta[b] - theta[a]) + 2*YtfC[a,b,i]*YttR[a,b,i]*v[a]*v[b]^3*sin(theta[b] - theta[a]) + YtfR[a,b,i]^2*v[a]^2*v[b]^2 - 2*YtfR[a,b,i]*YttC[a,b,i]*v[a]*v[b]^3*sin(theta[b] - theta[a]) + 2*YtfR[a,b,i]*YttR[a,b,i]*v[a]*v[b]^3*cos(theta[b] - theta[a]) + YttC[a,b,i]^2*v[b]^4 + YttR[a,b,i]^2*v[b]^4 <= SU[b,a,i]^2;

# bounds on phase difference
subject to phasediff1_pol{(b,a,i) in L0}: pdLB[b,a,i] <= theta[b] - theta[a];
subject to phasediff2_pol{(b,a,i) in L0}: theta[b] - theta[a] <= pdUB[b,a,i];

## reference bus: there had better be just one reference -- check in .run
#subject to reference_pol{b in B : busType[b] == 3}: theta[b] = 0;

## see acopf_polar-decl.mod
##problem acopf_polar: v, theta, SgenR, SgenC, gencost_pol, powerflowR_pol, powerflowC_pol, powerbound1_pol, powerbound2_pol, phasediff1_pol, phasediff2_pol, reference_pol;
