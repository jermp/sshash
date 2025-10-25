#!/bin/bash

echo "output log file =" $1
echo "output dir for indexes =" $2

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=Off
make -j

mkdir results-$1
mkdir results-$1/k31 results-$1/k63
mkdir $2

prefix=/mnt/hd2/pibiri/DNA

rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/cod.k31.eulertigs.fa.gz -k 31 -m 20 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/cod.k31.sshash >> results-$1/k31/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/kestrel.k31.eulertigs.fa.gz -k 31 -m 20 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/kestrel.k31.sshash >> results-$1/k31/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/human.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/human.k31.sshash >> results-$1/k31/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/axolotl.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/axolotl.k31.sshash >> results-$1/k31/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/hprc.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/hprc.k31.sshash >> results-$1/k31/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/ec.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/ec.k31.sshash >> results-$1/k31/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/se.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/se.k31.sshash >> results-$1/k31/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/ncbi-virus.k31.eulertigs.fa.gz -k 31 -m 19 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/ncbi-virus.k31.sshash >> results-$1/k31/regular-build.log

rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/cod.k31.eulertigs.fa.gz -k 31 -m 20 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/cod.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/kestrel.k31.eulertigs.fa.gz -k 31 -m 20 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/kestrel.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/human.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/human.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/axolotl.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/axolotl.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/hprc.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/hprc.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/ec.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/ec.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/se.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/se.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/ncbi-virus.k31.eulertigs.fa.gz -k 31 -m 19 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/ncbi-virus.k31.canon.sshash >> results-$1/k31/canon-build.log

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=On
make -j

rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/cod.k63.eulertigs.fa.gz -k 63 -m 24 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/cod.k63.sshash >> results-$1/k63/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/kestrel.k63.eulertigs.fa.gz -k 63 -m 24 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/kestrel.k63.sshash >> results-$1/k63/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/human.k63.eulertigs.fa.gz -k 63 -m 25 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/human.k63.sshash >> results-$1/k63/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/axolotl.k63.eulertigs.fa.gz -k 63 -m 25 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/axolotl.k63.sshash >> results-$1/k63/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/hprc.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/hprc.k63.sshash >> results-$1/k63/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/ec.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/ec.k63.sshash >> results-$1/k63/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/se.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/se.k63.sshash >> results-$1/k63/regular-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/ncbi-virus.k63.eulertigs.fa.gz -k 63 -m 23 -g 16 -t 64 --verbose -d $prefix/tmp_dir -o $2/ncbi-virus.k63.sshash >> results-$1/k63/regular-build.log

rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/cod.k63.eulertigs.fa.gz -k 63 -m 24 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/cod.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/kestrel.k63.eulertigs.fa.gz -k 63 -m 24 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/kestrel.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/human.k63.eulertigs.fa.gz -k 63 -m 25 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/human.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/axolotl.k63.eulertigs.fa.gz -k 63 -m 25 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/axolotl.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/hprc.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/hprc.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/ec.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/ec.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/se.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/se.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf $prefix/tmp_dir/*; ./sshash build -i $prefix/eulertigs/ncbi-virus.k63.eulertigs.fa.gz -k 63 -m 23 -g 16 -t 64 --canonical --verbose -d $prefix/tmp_dir -o $2/ncbi-virus.k63.canon.sshash >> results-$1/k63/canon-build.log

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=Off
make -j
