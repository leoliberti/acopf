# acopf
Code for ACOPF mathematical programming formulations

Minimal instructions:
solve an instance (see in instances/ and case/)
- edit the solver choice part in acopf.run
- run ./acopf.sh instance_file.dat model_file.mod
- example: ./acopf.sh instances/case5.dat acopf_SIV.mod
  (you need to have the AMPL version of the Baron solver installed)
  
syntax is ./acopf.sh instance.dat [model.mod [sp]]
  sp in {rnd, dat}
    rnd = random, dat = read Va,Vm from dat file
  default model is acopf_SIV.mod
  default starting point (sp) is rnd
