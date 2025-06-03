#!/usr/bin/env bash

for n in $(seq 1 576); do
# for n in $(seq 76 76); do
    input_file=$(printf "output_state_%04d.nc\n" $n)
    output_file=$(printf "data/output_state_%04d.nc\n" $n)
    echo "nccmp --force -d -g -t 1e-8 $input_file $output_file"
    nccmp -d -g -t 1e-8 $input_file $output_file
done

#for n in $(seq 1 576); do
#    input_file=$(printf "input_state_%04d.nc\n" $n)
#    output_file=$(printf "data/input_state_%04d.nc\n" $n)
#    echo "nccmp --force -d -g -t 1e-9 $input_file $output_file"
#    nccmp -d -g -t 1e-13 $input_file $output_file
#done
