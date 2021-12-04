#!/bin/bash

#i = $1
#file = $2

read -p "Enter start number " i
#read -p "Enter increment " increment
while [ $i -le 1000000000000 ]
do
  echo $i 
  (eval "python3 10000x10000RandomPoints.py $i && triangle -zn 10000x10000_$i.node && ./generatetrivertex 10000x10000_$i.1")
  #(eval "rm $i.node")
  i=$(($i*10))
done
