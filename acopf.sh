#!/bin/bash

MOD=acopf_SIV.mod
SP=rnd # random starting point

if [ "$1" == "" ]; then
    echo "syntax is $0 instance.dat [model.mod [sp]]"
    echo "  sp in {rnd, dat}"
    echo "    rnd = random, dat = read Va,Vm from dat file"
    echo "  default model is $MOD"
    echo "  default starting point (sp) is $SP"
    exit 1
fi

if [ "$2" != "" ]; then
    MOD=$2
fi

if [ "$3" != "" ]; then
    SP=$3
fi

modtype=`basename $MOD .mod | cut -d '_' -f 2`
if [[ "$modtype" =~ ^(SIV|V|VV)$ ]]; then
    # cartesian
    echo "$0: using cartesian model $MOD"
else
    # polar
    echo "$0: using polar model $MOD"
fi

ln -sf $1 acopf.dat
ln -sf $MOD acopf.mod

startpt=$SP time ampl acopf.run

rm -f acopf.dat acopf.mod 
