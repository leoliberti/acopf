## decision variables for all models in ACOPF

# voltage (amplitude, angle)
var v{b in B} <= VU[b], >= max(0,VL[b]);  # amplitude
var theta{B} <= 2*atan(1), >= -2*atan(1); # phase

# voltage (real, imaginary)
var VR{b in B} <= VU[b], >= -VU[b];  # real
var VC{b in B} <= VU[b], >= -VU[b];  # imaginary

# power generation (real, imaginary)
var SgenR{b in B, g in G[b]} >= SLR[b,g], <= SUR[b,g];
var SgenC{b in B, g in G[b]} >= SLC[b,g], <= SUC[b,g];

