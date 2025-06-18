#!/bin/bash


echo "output log file =" $1

### regular indexes

./sshash build -i ~/sshash_datasets/cod.k31.unitigs.fa.ust.fa.gz -k 31 -m 20 -t 8 -g 16 --verbose -o cod.k31.sshash -d tmp_dir >> $1.regular.build_log
# ./sshash build -i ~/sshash_datasets/cod.k63.unitigs.fa.ust.fa.gz -k 63 -m 24 -t 8 -g 16 --verbose -o cod.k63.sshash -d tmp_dir >> $1.regular.build_log

./sshash build -i ~/sshash_datasets/kestrel.k31.unitigs.fa.ust.fa.gz -k 31 -m 20 -t 8 -g 16 --verbose -o kestrel.k31.sshash -d tmp_dir >> $1.regular.build_log
# ./sshash build -i ~/sshash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz -k 63 -m 24 -t 8 -g 16 --verbose -o kestrel.k63.sshash -d tmp_dir >> $1.regular.build_log

./sshash build -i ~/sshash_datasets/human.k31.unitigs.fa.ust.fa.gz -k 31 -m 21 -t 8 -g 16 --verbose -o human.k31.sshash -d tmp_dir >> $1.regular.build_log
# ./sshash build -i ~/sshash_datasets/human.k63.unitigs.fa.ust.fa.gz -k 63 -m 25 -t 8 -g 16 --verbose -o human.k63.sshash -d tmp_dir >> $1.regular.build_log

### canonical indexes

./sshash build -i ~/sshash_datasets/cod.k31.unitigs.fa.ust.fa.gz -k 31 -m 19 -t 8 -g 16 --canonical --verbose -o cod.k31.canon.sshash -d tmp_dir >> $1.canon.build_log
# ./sshash build -i ~/sshash_datasets/cod.k63.unitigs.fa.ust.fa.gz -k 63 -m 23 -t 8 -g 16 --canonical --verbose -o cod.k63.canon.sshash -d tmp_dir >> $1.canon.build_log

./sshash build -i ~/sshash_datasets/kestrel.k31.unitigs.fa.ust.fa.gz -k 31 -m 19 -t 8 -g 16 --canonical --verbose -o kestrel.k31.canon.sshash -d tmp_dir >> $1.canon.build_log
# ./sshash build -i ~/sshash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz -k 63 -m 23 -t 8 -g 16 --canonical --verbose -o kestrel.k63.canon.sshash -d tmp_dir >> $1.canon.build_log

./sshash build -i ~/sshash_datasets/human.k31.unitigs.fa.ust.fa.gz -k 31 -m 20 -t 8 -g 16 --canonical --verbose -o human.k31.canon.sshash -d tmp_dir >> $1.canon.build_log
# ./sshash build -i ~/sshash_datasets/human.k63.unitigs.fa.ust.fa.gz -k 63 -m 24 -t 8 -g 16 --canonical --verbose -o human.k63.canon.sshash -d tmp_dir >> $1.canon.build_log