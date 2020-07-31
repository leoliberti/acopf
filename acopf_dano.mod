## ACOPF formulation
## acopf_dano.mod AMPL file
## Leo Liberti
## This AMPL model implements a formulation described in
##   ../../papers/acopf/formulations/4or/acopf_survey.pdf 
##    it can integrate shunts, transformers or both
##   Formulation in real numbers (rather than complex),
##     decision variables: generated power, voltage
## 191031  derived from acopf_SIV.mod through replacement
##         I=YV with polar coordinates, then S=VI^*
##         expressed in function of polar variables v_b, theta_b
## 191206  Dan's formulation: Jabr SOCP relaxation mixed with SIV formulation

let formulation_type := "cartesian";

## common parameters read from .run

## local parameters

# cos/sin representation variable index set
set R := {b in B, a in B : (b,a,1) in L or b==a};

## common decision variables read from .run

## local decision variables

# c[b,a] represents v_b*v_a*cos(theta_b - theta_a)
var c{(b,a) in R} >= if b==a then VL[b]^2 else -VU[b]*VU[a], <= if b==a then VU[b]^2 else VU[b]*VU[a];
# s[b,a] represents v_b*v_a*sin(theta_b - theta_a)
var s{(b,a) in R} >= if b==a then 0 else -VU[b]*VU[a], <= if b==a then 0 else VU[b]*VU[a];

## objective function

# generation cost (WARNING: to remove quadratic power terms set C[g,2]=0)
minimize gencost_dano: sum{b in B, g in G[b]}
  (C[b,g,2]*SgenR[b,g]^2 + C[b,g,1]*SgenR[b,g] + C[b,g,0]);

### constraints

## power flow balance
subject to powerflowR_dano{b in B}:  # real
   SDR[b] +
   sum{(b,a,i) in L0}
     (YffR[b,a,i]*c[b,b] + YftC[b,a,i]*s[b,a] + YftR[b,a,i]*c[b,a]) +
   sum{(b,a,i) in L1}
     (YtfC[a,b,i]*s[b,a] + YtfR[a,b,i]*c[b,a] + YttR[a,b,i]*c[b,b]) =
   -shR[b]*c[b,b] + sum{g in G[b]} SgenR[b,g];
subject to powerflowC_dano{b in B}:  # imaginary
   SDC[b] +
   sum{(b,a,i) in L0}
     (-YffC[b,a,i]*c[b,b] - YftC[b,a,i]*c[b,a] + YftR[b,a,i]*s[b,a]) +
   sum{(b,a,i) in L1}
     (-YtfC[a,b,i]*c[b,a] + YtfR[a,b,i]*s[b,a] - YttC[a,b,i]*c[b,b]) =
    shC[b]*c[b,b] + sum{g in G[b]} SgenC[b,g];

## power bound on lines b->a defined on |S|^2 <= SU^2: (computed using sympy)
subject to powerbound1_dano{(b,a,i) in L0 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  YffC[b,a,i]^2*c[b,b]^2 + 2*YffC[b,a,i]*YftC[b,a,i]*c[b,b]*c[b,a] - 2*YffC[b,a,i]*YftR[b,a,i]*c[b,b]*s[b,a] + YffR[b,a,i]^2*c[b,b]^2 + 2*YffR[b,a,i]*YftC[b,a,i]*c[b,b]*s[b,a] + 2*YffR[b,a,i]*YftR[b,a,i]*c[b,b]*c[b,a] + YftC[b,a,i]^2*c[b,b]*c[a,a] + YftR[b,a,i]^2*c[b,b]*c[a,a] <= SU[b,a,i]^2;

subject to powerbound2_dano{(b,a,i) in L1 : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  YtfC[a,b,i]^2*c[b,b]*c[a,a] + 2*YtfC[a,b,i]*YttC[a,b,i]*c[b,b]*c[b,a] + 2*YtfC[a,b,i]*YttR[a,b,i]*c[b,b]*s[b,a] + YtfR[a,b,i]^2*c[b,b]*c[a,a] - 2*YtfR[a,b,i]*YttC[a,b,i]*c[b,b]*s[b,a] + 2*YtfR[a,b,i]*YttR[a,b,i]*c[b,b]*c[b,a] + YttC[a,b,i]^2*c[b,b]^2 + YttR[a,b,i]^2*c[b,b]^2 <= SU[b,a,i]^2;

# bounds on phase difference
subject to phasediff1_dano{(b,a,i) in L0}: tan(pdLB[b,a,i])*c[b,a] <= s[b,a];
subject to phasediff2_dano{(b,a,i) in L0}: s[b,a] <= tan(pdUB[b,a,i])*c[b,a];

## reference bus: there had better be just one reference -- check in .run
subject to reference_dano{b in B : busType[b] == 3}: VC[b] = 0;

## Jabr relaxation

# symmetry
subject to symm1_dano{(b,a,1) in L0}: c[b,a] =  c[a,b];
subject to symm2_dano{(b,a,1) in L0}: s[b,a] = -s[a,b];

# relaxation
subject to relax_dano1{(b,a,1) in L0}: c[b,a]^2 + s[b,a]^2 = c[b,b]*c[a,a];
subject to relax_dano2{b in B}: v[b]^2 <= c[b,b];
subject to relax_dano3{b in B}: s[b,b] = 0;

## links with V variables
subject to csVrel1_dano{b in B}: c[b,b] = VR[b]^2 + VC[b]^2;
subject to csVrel2_dano{(b,a,1) in L0}: c[b,a] = VR[b]*VR[a] + VC[b]*VC[a];
subject to csVrel3_dano{(b,a,1) in L0}: s[b,a] = VC[b]*VR[a] - VR[b]*VC[a];

## links with v,theta variables
#subject to csvrel1_dano{b in B}: c[b,b] = v[b]^2;
#subject to csvrel2_dano{(b,a,1) in L0}:c[b,a]=v[b]*v[a]*cos(theta[b]-theta[a]);
#subject to csvrel3_dano{(b,a,1) in L0}:s[b,a]=v[b]*v[a]*sin(theta[b]-theta[a]);

