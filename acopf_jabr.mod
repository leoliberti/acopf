## ACOPF formulation
## acopf_jabr.mod AMPL file
## Leo Liberti
## This AMPL model implements a formulation described in
##   ../../papers/acopf/formulations/4or/acopf_survey.pdf 
##    it can integrate shunts, transformers or both
##   Formulation in real numbers (rather than complex),
##     decision variables: generated power, voltage
## 191031  derived from acopf_SIV.mod through replacement
##         I=YV with polar coordinates, then S=VI^*
##         expressed in function of polar variables v_b, theta_b
## 191204  Jabr SOCP relaxation derived from polar formulation acopf_V_polar

let formulation_type := "polar";

## common parameters read from .run

## local parameters

# cos/sin representation variable index set
set R := {b in B, a in B : (b,a,1) in L or b==a};

# starting point (used for verifying SDP relaxation solution)
param csp{R} default 0;
param ssp{R} default 0;

## common decision variables read from .run

## local decision variables

# c[b,a] represents v_b*v_a*cos(theta_b - theta_a)
var c{(b,a) in R} >= if b==a then VL[b]^2 else 0, <= if b==a then VU[b]^2 else VU[b]*VU[a];
# s[b,a] represents v_b*v_a*sin(theta_b - theta_a)
var s{(b,a) in R} >= if b==a then 0 else -VU[b]*VU[a], <= if b==a then 0 else VU[b]*VU[a];

## objective function

# generation cost (WARNING: to remove quadratic power terms set C[g,2]=0)
minimize gencost_jabr: sum{b in B, g in G[b]}
  (C[b,g,2]*SgenR[b,g]^2 + C[b,g,1]*SgenR[b,g] + C[b,g,0]);

### constraints

## power flow balance
subject to powerflowR_jabr{b in B}:  # real
   SDR[b] +
   sum{(b,a,i) in L0}
     (YffR[b,a,i]*c[b,b] + YftC[b,a,i]*s[b,a] + YftR[b,a,i]*c[b,a]) +
   sum{(b,a,i) in L1}
     (YtfC[a,b,i]*s[b,a] + YtfR[a,b,i]*c[b,a] + YttR[a,b,i]*c[b,b]) =
   -shR[b]*c[b,b] + sum{g in G[b]} SgenR[b,g];
subject to powerflowC_jabr{b in B}:  # imaginary
   SDC[b] +
   sum{(b,a,i) in L0}
     (-YffC[b,a,i]*c[b,b] + YftR[b,a,i]*s[b,a] - YftC[b,a,i]*c[b,a]) +
   sum{(b,a,i) in L1}
     (-YttC[a,b,i]*c[b,b] + YtfR[a,b,i]*s[b,a] - YtfC[a,b,i]*c[b,a]) =
    shC[b]*c[b,b] + sum{g in G[b]} SgenC[b,g];

## power bound on lines b->a defined on |S|^2 <= SU^2: (computed using sympy)
subject to powerbound1_jabr{(b,a,i) in L0 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  YffC[b,a,i]^2*c[b,b]^2 + 2*YffC[b,a,i]*YftC[b,a,i]*c[b,b]*c[b,a] - 2*YffC[b,a,i]*YftR[b,a,i]*c[b,b]*s[b,a] + YffR[b,a,i]^2*c[b,b]^2 + 2*YffR[b,a,i]*YftC[b,a,i]*c[b,b]*s[b,a] + 2*YffR[b,a,i]*YftR[b,a,i]*c[b,b]*c[b,a] + YftC[b,a,i]^2*c[b,b]*c[a,a] + YftR[b,a,i]^2*c[b,b]*c[a,a] <= SU[b,a,i]^2;

subject to powerbound2_jabr{(b,a,i) in L1 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  YtfC[a,b,i]^2*c[b,b]*c[a,a] + 2*YtfC[a,b,i]*YttC[a,b,i]*c[b,b]*c[b,a] + 2*YtfC[a,b,i]*YttR[a,b,i]*c[b,b]*s[b,a] + YtfR[a,b,i]^2*c[b,b]*c[a,a] - 2*YtfR[a,b,i]*YttC[a,b,i]*c[b,b]*s[b,a] + 2*YtfR[a,b,i]*YttR[a,b,i]*c[b,b]*c[b,a] + YttC[a,b,i]^2*c[b,b]^2 + YttR[a,b,i]^2*c[b,b]^2 <= SU[b,a,i]^2;

# bounds on phase difference
subject to phasediff1_jabr{(b,a,i) in L0}: tan(pdLB[b,a,i])*c[b,a] <= s[b,a];
subject to phasediff2_jabr{(b,a,i) in L0}: s[b,a] <= tan(pdUB[b,a,i])*c[b,a];

# reference bus: there had better be just one reference -- check in .run
subject to reference_jabr{b in B : busType[b] == 3}: theta[b] = 0;

## original defining constraints
#subject to cdeforg_jabr{(b,a) in R}: c[b,a] = v[b]*v[a]*cos(theta[b]-theta[a]);
#subject to sdeforg_jabr{(b,a) in R}: s[b,a] = v[b]*v[a]*sin(theta[b]-theta[a]);

# ## derived defining constraints
# subject to defcon1_jabr{(b,a,1) in L0}: c[b,a]^2 + s[b,a]^2 = c[b,b]*c[a,a];
# subject to defcon2_jabr{b in B}: v[b]^2 = c[b,b];
# subject to defcon3_jabr{b in B}: s[b,b] = 0;

## Jabr relaxation

# symmetry
subject to symm1_jabr{(b,a,1) in L0}: c[b,a] =  c[a,b];
subject to symm2_jabr{(b,a,1) in L0}: s[b,a] = -s[a,b];

# relaxation
subject to relax_jabr1{(b,a,1) in L0}: c[b,a]^2 + s[b,a]^2 <= c[b,b]*c[a,a];
subject to relax_jabr2{b in B}: v[b]^2 <= c[b,b];
subject to relax_jabr3{b in B}: s[b,b] = 0;

