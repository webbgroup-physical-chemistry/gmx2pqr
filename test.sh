#1 /bin/bash

clear
rm gmx2pqr examples/state0.pqr examples/field_projection.xvg examples/state*.pqr
make -j 8
if [ ! -f gmx2pqr ] ; then
    exit
    fi
if [ ! -z $1 ] ; then
    ./gmx2pqr -h
else
    gro=examples/gro.gro
    CD=`grep "CNC" $gro | grep CD | awk '{print $3}'`
    NE=`grep "CNC" $gro | grep NE | awk '{print $3}'`
    SG=`grep "CNC" $gro | grep SG | awk '{print $3}'`
    CB=`grep "CNC" $gro | grep CB | awk '{print $3}'`
    HB1=`grep "CNC" $gro | grep HB1 | awk '{print $3}'`
    HB2=`grep "CNC" $gro | grep HB2 | awk '{print $3}'`
time ./gmx2pqr -s examples/tpr.tpr -f examples/xtc.xtc -o examples/state.pqr -of examples/field_projection.xvg -a1 $CD -a2 $NE -exclude "$SG $CB $HB1 $HB2" -site "$CB $SG" -select "(not resname SOL and not name Na)"  -bw examples/bw.dat -b 5301 -e 6000 -nodopqr
#cat examples/field_projection.xvg
#grep "CNC\|DUM" examples/state0.pqr
    fi

