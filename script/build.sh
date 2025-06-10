#!/bin/bash

### regular indexes

./sshash build -i /data2/DNA/lphash_datasets/celegans.k31.unitigs.fa.ust.fa.gz -k 31 -m 16 --verbose -o celegans.k31.sshash -d tmp_dir >> sshash.regular.build_log
./sshash build -i /data2/DNA/lphash_datasets/celegans.k47.unitigs.fa.ust.fa.gz -k 47 -m 20 --verbose -o celegans.k47.sshash -d tmp_dir >> sshash.regular.build_log
./sshash build -i /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz -k 63 -m 20 --verbose -o celegans.k63.sshash -d tmp_dir >> sshash.regular.build_log

./sshash build -i /data2/DNA/lphash_datasets/cod.k31.unitigs.fa.ust.fa.gz -k 31 -m 20 --verbose -o cod.k31.sshash -d tmp_dir >> sshash.regular.build_log
./sshash build -i /data2/DNA/lphash_datasets/cod.k47.unitigs.fa.ust.fa.gz -k 47 -m 22 --verbose -o cod.k47.sshash -d tmp_dir >> sshash.regular.build_log
./sshash build -i /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz -k 63 -m 24 --verbose -o cod.k63.sshash -d tmp_dir >> sshash.regular.build_log

./sshash build -i /data2/DNA/lphash_datasets/kestrel.k31.unitigs.fa.ust.fa.gz -k 31 -m 20 --verbose -o kestrel.k31.sshash -d tmp_dir >> sshash.regular.build_log
./sshash build -i /data2/DNA/lphash_datasets/kestrel.k47.unitigs.fa.ust.fa.gz -k 47 -m 22 --verbose -o kestrel.k47.sshash -d tmp_dir >> sshash.regular.build_log
./sshash build -i /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz -k 63 -m 24 --verbose -o kestrel.k63.sshash -d tmp_dir >> sshash.regular.build_log

./sshash build -i /data2/DNA/lphash_datasets/human.k31.unitigs.fa.ust.fa.gz -k 31 -m 21 --verbose -o human.k31.sshash -d tmp_dir >> sshash.regular.build_log
./sshash build -i /data2/DNA/lphash_datasets/human.k47.unitigs.fa.ust.fa.gz -k 47 -m 23 --verbose -o human.k47.sshash -d tmp_dir >> sshash.regular.build_log
./sshash build -i /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz -k 63 -m 25 --verbose -o human.k63.sshash -d tmp_dir >> sshash.regular.build_log

### canonical indexes

./sshash build -i /data2/DNA/lphash_datasets/celegans.k31.unitigs.fa.ust.fa.gz -k 31 -m 15 --canonical --verbose -o celegans.k31.canon.sshash -d tmp_dir >> sshash.canon.build_log
./sshash build -i /data2/DNA/lphash_datasets/celegans.k47.unitigs.fa.ust.fa.gz -k 47 -m 19 --canonical --verbose -o celegans.k47.canon.sshash -d tmp_dir >> sshash.canon.build_log
./sshash build -i /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz -k 63 -m 19 --canonical --verbose -o celegans.k63.canon.sshash -d tmp_dir >> sshash.canon.build_log

./sshash build -i /data2/DNA/lphash_datasets/cod.k31.unitigs.fa.ust.fa.gz -k 31 -m 19 --canonical --verbose -o cod.k31.canon.sshash -d tmp_dir >> sshash.canon.build_log
./sshash build -i /data2/DNA/lphash_datasets/cod.k47.unitigs.fa.ust.fa.gz -k 47 -m 21 --canonical --verbose -o cod.k47.canon.sshash -d tmp_dir >> sshash.canon.build_log
./sshash build -i /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz -k 63 -m 23 --canonical --verbose -o cod.k63.canon.sshash -d tmp_dir >> sshash.canon.build_log

./sshash build -i /data2/DNA/lphash_datasets/kestrel.k31.unitigs.fa.ust.fa.gz -k 31 -m 19 --canonical --verbose -o kestrel.k31.canon.sshash -d tmp_dir >> sshash.canon.build_log
./sshash build -i /data2/DNA/lphash_datasets/kestrel.k47.unitigs.fa.ust.fa.gz -k 47 -m 21 --canonical --verbose -o kestrel.k47.canon.sshash -d tmp_dir >> sshash.canon.build_log
./sshash build -i /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz -k 63 -m 23 --canonical --verbose -o kestrel.k63.canon.sshash -d tmp_dir >> sshash.canon.build_log

./sshash build -i /data2/DNA/lphash_datasets/human.k31.unitigs.fa.ust.fa.gz -k 31 -m 20 --canonical --verbose -o human.k31.canon.sshash -d tmp_dir >> sshash.canon.build_log
./sshash build -i /data2/DNA/lphash_datasets/human.k47.unitigs.fa.ust.fa.gz -k 47 -m 22 --canonical --verbose -o human.k47.canon.sshash -d tmp_dir >> sshash.canon.build_log
./sshash build -i /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz -k 63 -m 24 --canonical --verbose -o human.k63.canon.sshash -d tmp_dir >> sshash.canon.build_log