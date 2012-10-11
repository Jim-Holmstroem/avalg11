#!/bin/bash

for ((i=4;i<48;i=i+4))
do
    for ((j=4;j<=i;j=j+4))
    do
	     make run_test SIZE="$i"x"$j"
    done
done

exit 0
