## ACOPF formulation
## acopf.mod AMPL file
## Leo Liberti
## This AMPL model implements a formulation described in
##   ../../papers/acopf/formulations/4or/acopf_survey.pdf 
##    it can integrate shunts, transformers or both
##   Formulation in real numbers (rather than complex), decision variables:
##     power,current,voltage; current,power completely determined by voltage
## NOTE: this formulation uses injected power bounds, so it is "exact"
## 180907 work started (IASI@CNR) -- old version with A,B coefficients
##        Dan Bienstock found serious mistake in expression for Ohm's law
## 190925 major overhaul after consulting with Dan Bienstock @Columbia:
##        see "flow" as "influence", consider half-edges; also, more than one
##        generator can be at a single node
## 191025 as in MatPower cases, we have power magnitude bounds on lines
##        but not current bounds

## matpower: extract VR, VI from Vm, Vang in results = runopf(case5, mpopt)
## V = [results.bus(:, 8).*cos(results.bus(:,9)*(pi/180)), results.bus(:, 8).*sin(results.bus(:,9)*(pi/180))]

let formulation_type := "cartesian";

## common parameters and variables read in run file

## local decision variables

# V2 = |V|^2
var V2{b in B} >= VL[b]^2, <= VU[b]^2;

# current (real, imaginary)
var IR{L};
var IC{L};

# power injected on line at bus (real, imaginary)
var SR{L};
var SC{L};

## objective function

# generation cost (WARNING: to remove quadratic power terms set C[g,2]=0)
minimize gencost: sum{b in B, g in G[b]} (C[b,g,2]*SgenR[b,g]^2 + C[b,g,1]*SgenR[b,g] + C[b,g,0]);

# # reactive power loss
# minimize reactive_power_loss: sum{b in B, g in G[b]} (C[b,g,1]*SgenC[b,g] + C[b,g,0]);

## constraints

# power flow (real, imaginary)
subject to powerflowR{b in B}:
   SDR[b] + sum{(b,a,i) in L} SR[b,a,i] = -shR[b]*V2[b] + sum{g in G[b]} SgenR[b,g];
subject to powerflowC{b in B}:
   SDC[b] + sum{(b,a,i) in L} SC[b,a,i] =  shC[b]*V2[b] + sum{g in G[b]} SgenC[b,g];

# definition of power in function of current and voltage (real, imaginary)
subject to powerinjR{(b,a,i) in L}: SR[b,a,i] = VR[b]*IR[b,a,i] + VC[b]*IC[b,a,i];
subject to powerinjC{(b,a,i) in L}: SC[b,a,i] = VC[b]*IR[b,a,i] - VR[b]*IC[b,a,i];

# Ohm's law ((b,a):real,imaginary; (a,b):real,imaginary)
subject to ohm1R{(b,a,i) in L0}:
  IR[b,a,i] = YffR[b,a,i]*VR[b]-YffC[b,a,i]*VC[b] + YftR[b,a,i]*VR[a]-YftC[b,a,i]*VC[a];
subject to ohm1C{(b,a,i) in L0}:
  IC[b,a,i] = YffR[b,a,i]*VC[b]+YffC[b,a,i]*VR[b] + YftR[b,a,i]*VC[a]+YftC[b,a,i]*VR[a];
subject to ohm2R{(b,a,i) in L0}:
  IR[a,b,i] = YtfR[b,a,i]*VR[b]-YtfC[b,a,i]*VC[b] + YttR[b,a,i]*VR[a]-YttC[b,a,i]*VC[a];
subject to ohm2C{(b,a,i) in L0}:
  IC[a,b,i] = YtfR[b,a,i]*VC[b]+YtfC[b,a,i]*VR[b] + YttR[b,a,i]*VC[a]+YttC[b,a,i]*VR[a];

# power bound on lines b->a defined on I
subject to powerbound{(b,a,i) in L : SU[b,a,i]>0 and SU[b,a,i]<Inf}:
  SR[b,a,i]^2 + SC[b,a,i]^2 <= SU[b,a,i]^2;

# definition of V2
subject to V2def{b in B}: V2[b] = VR[b]^2 + VC[b]^2;

# bounds on phase difference
subject to phasediff1{(b,a,i) in L0}:
  VC[b]*VR[a] - VR[b]*VC[a] <= tan(pdUB[b,a,i]) * (VR[b]*VR[a] + VC[b]*VC[a]);
subject to phasediff2{(b,a,i) in L0}:
  VC[b]*VR[a] - VR[b]*VC[a] >= tan(pdLB[b,a,i]) * (VR[b]*VR[a] + VC[b]*VC[a]);
subject to phasediff3{(b,a,1) in L0}: VR[b]*VR[a] + VC[b]*VC[a] >= 0;

# reference bus: there had better be just one reference -- check in .run
subject to reference1{b in B : busType[b] == 3}: VC[b] = 0;
subject to reference2{b in B : busType[b] == 3}: VR[b] >= 0;


## see acopf_SIV-decl.mod
##problem acopf_SIV: VR, VC, V2, SgenR, SgenC, IR, IC, SR, SC, gencost, powerflowR, powerflowC, powerinjR, powerinjC, ohm1R, ohm1C, ohm2R, ohm2C, powerbound, V2def, phasediff1, phasediff2, phasediff3, reference1, reference2;
