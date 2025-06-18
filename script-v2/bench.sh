#!/bin/bash

echo "output log file =" $1

### regular indexes

./sshash bench -i cod.k31.sshash >> $1.regular.bench_log
# ./sshash bench -i cod.k63.sshash >> $1.regular.bench_log

./sshash bench -i kestrel.k31.sshash >> $1.regular.bench_log
# ./sshash bench -i kestrel.k63.sshash >> $1.regular.bench_log

./sshash bench -i human.k31.sshash >> $1.regular.bench_log
# ./sshash bench -i human.k63.sshash >> $1.regular.bench_log

### canonical indexes

./sshash bench -i cod.k31.canon.sshash >> $1.canon.bench_log
# ./sshash bench -i cod.k63.canon.sshash >> $1.canon.bench_log

./sshash bench -i kestrel.k31.canon.sshash >> $1.canon.bench_log
# ./sshash bench -i kestrel.k63.canon.sshash >> $1.canon.bench_log

./sshash bench -i human.k31.canon.sshash >> $1.canon.bench_log
# ./sshash bench -i human.k63.canon.sshash >> $1.canon.bench_log