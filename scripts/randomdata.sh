#!/bin/bash

i=$1

while [ $i -le 500000 ]
do
  echo $i 
  (eval "python3 10000x10000RandomPoints.py $i && triangle -zn 10000x10000_$i.node && ./generatetrivertex 10000x10000_$i.1")
  i=$(($i + 1000000))
done
