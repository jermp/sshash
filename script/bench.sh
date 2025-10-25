#!/bin/bash

echo "output log file =" $1
echo "input dir for indexes =" $2

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=Off
make -j

mkdir results-$1
mkdir results-$1/k31 results-$1/k63
mkdir $2

prefix=/mnt/hd2/pibiri/DNA

for i in {1..3}; do
    ./sshash bench -i $2/cod.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/kestrel.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/human.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/axolotl.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/hprc.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/ec.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/se.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/ncbi-virus.k31.sshash >> results-$1/k31/regular-bench.log
done

for i in {1..3}; do
    ./sshash bench -i $2/cod.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/kestrel.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/human.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/axolotl.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/hprc.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/ec.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/se.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/ncbi-virus.k31.canon.sshash >> results-$1/k31/canon-bench.log
done

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=On
make -j

for i in {1..3}; do
    ./sshash bench -i $2/cod.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/kestrel.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/human.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/axolotl.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/hprc.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/ec.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/se.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/ncbi-virus.k63.sshash >> results-$1/k63/regular-bench.log
done

for i in {1..3}; do
    ./sshash bench -i $2/cod.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/kestrel.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/human.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/axolotl.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/hprc.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/ec.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/se.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i $2/ncbi-virus.k63.canon.sshash >> results-$1/k63/canon-bench.log
done

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=Off
make -j
