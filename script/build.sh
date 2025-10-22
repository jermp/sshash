#!/bin/bash

echo "output log file =" $1

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=Off
make -j

mkdir results-$1
mkdir results-$1/k31 results-$1/k63

rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/cod.k31.eulertigs.fa.gz -k 31 -m 20 -g 16 -t 64 --verbose -d tmp_dir -o cod.k31.sshash >> results-$1/k31/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/kestrel.k31.eulertigs.fa.gz -k 31 -m 20 -g 16 -t 64 --verbose -d tmp_dir -o kestrel.k31.sshash >> results-$1/k31/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/human.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d tmp_dir -o human.k31.sshash >> results-$1/k31/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/hprc.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d tmp_dir -o hprc.k31.sshash >> results-$1/k31/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/ec.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d tmp_dir -o ec.k31.sshash >> results-$1/k31/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/se.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --verbose -d tmp_dir -o se.k31.sshash >> results-$1/k31/regular-build.log

rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/cod.k31.eulertigs.fa.gz -k 31 -m 20 -g 16 -t 64 --canonical --verbose -d tmp_dir -o cod.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/kestrel.k31.eulertigs.fa.gz -k 31 -m 20 -g 16 -t 64 --canonical --verbose -d tmp_dir -o kestrel.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/human.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d tmp_dir -o human.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/hprc.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d tmp_dir -o hprc.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/ec.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d tmp_dir -o ec.k31.canon.sshash >> results-$1/k31/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/se.k31.eulertigs.fa.gz -k 31 -m 21 -g 16 -t 64 --canonical --verbose -d tmp_dir -o se.k31.canon.sshash >> results-$1/k31/canon-build.log

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=On
make -j

rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/cod.k63.eulertigs.fa.gz -k 63 -m 24 -g 16 -t 64 --verbose -d tmp_dir -o cod.k63.sshash >> results-$1/k63/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/kestrel.k63.eulertigs.fa.gz -k 63 -m 24 -g 16 -t 64 --verbose -d tmp_dir -o kestrel.k63.sshash >> results-$1/k63/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/human.k63.eulertigs.fa.gz -k 63 -m 25 -g 16 -t 64 --verbose -d tmp_dir -o human.k63.sshash >> results-$1/k63/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/hprc.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --verbose -d tmp_dir -o hprc.k63.sshash >> results-$1/k63/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/ec.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --verbose -d tmp_dir -o ec.k63.sshash >> results-$1/k63/regular-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/se.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --verbose -d tmp_dir -o se.k63.sshash >> results-$1/k63/regular-build.log

rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/cod.k63.eulertigs.fa.gz -k 63 -m 24 -g 16 -t 64 --canonical --verbose -d tmp_dir -o cod.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/kestrel.k63.eulertigs.fa.gz -k 63 -m 24 -g 16 -t 64 --canonical --verbose -d tmp_dir -o kestrel.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/human.k63.eulertigs.fa.gz -k 63 -m 25 -g 16 -t 64 --canonical --verbose -d tmp_dir -o human.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/hprc.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --canonical --verbose -d tmp_dir -o hprc.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/ec.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --canonical --verbose -d tmp_dir -o ec.k63.canon.sshash >> results-$1/k63/canon-build.log
rm -rf tmp_dir/; ./sshash build -i /mnt/hd2/pibiri/DNA/eulertigs/se.k63.eulertigs.fa.gz -k 63 -m 31 -g 16 -t 64 --canonical --verbose -d tmp_dir -o se.k63.canon.sshash >> results-$1/k63/canon-build.log

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DSSHASH_USE_ARCH_NATIVE=On -DSSHASH_USE_SANITIZERS=Off -DSSHASH_USE_MAX_KMER_LENGTH_63=Off
make -j
