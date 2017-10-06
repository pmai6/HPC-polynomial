#!/bin/bash
make all

read step
for i in `seq 1 24`;
do
        echo $i
        mpirun -np $i ./poly-eval sample-constants-long.txt sample-values.txt > nproc-$i-run1-$step.dat
        mpirun -np $i ./poly-eval sample-constants-long.txt sample-values.txt > nproc-$i-run2-$step.dat
        mpirun -np $i ./poly-eval sample-constants-long.txt sample-values.txt > nproc-$i-run3-$step.dat
        mpirun -np $i ./poly-eval sample-constants-long.txt sample-values.txt > nproc-$i-run4-$step.dat
        mpirun -np $i ./poly-eval sample-constants-long.txt sample-values.txt > nproc-$i-run5-$step.dat
done

rm -r $step
mkdir $step
mv *.dat $step/
