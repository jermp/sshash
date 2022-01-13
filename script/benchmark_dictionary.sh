#!/bin/bash

echo "input file =" $1
echo "k =" $2
echo "m =" $3

./build $1 $2 $3 -o out.index
./build $1 $2 $3 --canonical-parsing -o out.canon.index

./bench out.index
./bench out.canon.index
