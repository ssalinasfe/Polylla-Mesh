#!/bin/bash

i=10

echo "pnumber, tnumber, num_regs, edges, max_edges, min_edges, avg_edges, time"
while :
do
  (eval "rbox $i D2 z B1e3 > data.dat  && ./uwu data.dat 0 dat")
  i=$(($i * 10))
done