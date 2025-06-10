#!/bin/bash

echo "input file =" $1
echo "k =" $2
echo "m =" $3

./sshash build -i $1 -k $2 -m $3 -o index.sshash
./sshash build -i $1 -k $2 -m $3 --canonical -o index.canon.sshash

./sshash bench -i index.sshash
./sshash bench -i index.canon.sshash
