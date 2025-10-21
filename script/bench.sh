#!/bin/bash

echo "output log file =" $1

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=Off
make -j

mkdir results-$1
mkdir results-$1/k31 results-$1/k63


for i in {1..3}; do
    ./sshash bench -i cod.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i kestrel.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i human.k31.sshash >> results-$1/k31/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i hprc.k31.sshash >> results-$1/k31/regular-bench.log
done

for i in {1..3}; do
    ./sshash bench -i cod.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i kestrel.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i human.k31.canon.sshash >> results-$1/k31/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i hprc.k31.canon.sshash >> results-$1/k31/canon-bench.log
done

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=On
make -j

for i in {1..3}; do
    ./sshash bench -i cod.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i kestrel.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i human.k63.sshash >> results-$1/k63/regular-bench.log
done
for i in {1..3}; do
    ./sshash bench -i hprc.k63.sshash >> results-$1/k63/regular-bench.log
done

for i in {1..3}; do
    ./sshash bench -i cod.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i kestrel.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i human.k63.canon.sshash >> results-$1/k63/canon-bench.log
done
for i in {1..3}; do
    ./sshash bench -i hprc.k63.canon.sshash >> results-$1/k63/canon-bench.log
done

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=Off
make -j
