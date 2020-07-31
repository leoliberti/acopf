
#### user-configurable ####
param Inf := 1e30;
param Eps := 1e-6;
param Pi := 4*atan(1);
param myZero := 1e-9;

## sets

# max number of parallel branches in instance
param maxParBranches integer, >0, default 1;
set PB := 1..maxParBranches;

# buses
set B;

## lines (directed: edges correspond to antisymmetric arcs)
## 191026: can't deal with parallel arcs (eg case57)
# the given edges (b,a): the transformer is on b
set L0 within {B,B,PB};  
# set of all antiparallel arcs
set L default L0 union {(a,b,i) in {B,B,PB} : (b,a,i) in L0}; 
set L1 default L diff L0;

# set of generators at node
set G{b in B} default {};

## parameters

# bus type (2=generator, 3=reference)
param busType{B} integer;

# cost coefficients (minimization of power generation costs)
param Kcard integer, >= 0, default 2;
set K := 0..Kcard;
# initialization: only linear terms in P (quadratic in V)
param C{b in B, g in G[b], k in K} default if k == 1 then 1 else 0;

# real power demand at buses (there can be nodes with negative demands)
param SDR{B} default 0; 
# reactive power demand at buses (appears in MATPOWER documentation)
param SDC{B} default 0;

# power bounds at generators (real, imaginary)
param SLR{b in B, G[b]} default -Inf;
param SLC{b in B, G[b]} default -Inf;
param SUR{b in B, g in G[b]} >= SLR[b,g], default Inf;
param SUC{b in B, g in G[b]} >= SLC[b,g], default Inf;

# upper power magnitude bounds on links - symmetric
param SU{L} >= 0, default Inf; 

# voltage magnitude bounds at buses
param VL{B} default 0; # can't have negative moduli
param VU{b in B} >= VL[b], default Inf;
param Vm{B} default 1; # starting point for polar formulation
param Va{B} default 0; # starting point for polar formulation

# shunt parameters at buses
param shR{B} default 0;  # MatPower's Gs
param shC{B} default 0;  # MatPower's Bs

# status of a branch
param status{L} default 1;

# Y matrix (Ohm's law in AC) data, only defined on given arcs in L0
param r{L0} default 0;
param x{L0} default 0;
param bb{L0} default 0;
param tau{L0} default 1;
param nu{L0} default 0;  # translated to radians by mpc2dat.py

# phase difference bounds (only across lines in L0)
#   translated to radians by mpc2dat.py ##191026: not using these
param pdLB{L0} default -Pi;
param pdUB{L0} default Pi;

# the 2x2 complex Y matrix appearing in Ohm's law for a line (b,a) in L0
#   (computation in .run file)
param YffR{L0} default 0;
param YffC{L0} default 0;
param YftR{L0} default 0;
param YftC{L0} default 0;
param YtfR{L0} default 0;
param YtfC{L0} default 0;
param YttR{L0} default 0;
param YttC{L0} default 0;

# used in cartesian representation
param Vre{B} default 0;
param Vim{B} default 0;

# used in polar representation
param Vmag{b in B} default 0;
param Vang{b in B} default 0;

