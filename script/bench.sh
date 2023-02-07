#!/bin/bash

### regular indexes

./sshash bench -i celegans.k31.sshash >> sshash.benchmarking_log
./sshash bench -i celegans.k47.sshash >> sshash.benchmarking_log
./sshash bench -i celegans.k63.sshash >> sshash.benchmarking_log

./sshash bench -i cod.k31.sshash >> sshash.benchmarking_log
./sshash bench -i cod.k47.sshash >> sshash.benchmarking_log
./sshash bench -i cod.k63.sshash >> sshash.benchmarking_log

./sshash bench -i kestrel.k31.sshash >> sshash.benchmarking_log
./sshash bench -i kestrel.k47.sshash >> sshash.benchmarking_log
./sshash bench -i kestrel.k63.sshash >> sshash.benchmarking_log

./sshash bench -i human.k31.sshash >> sshash.benchmarking_log
./sshash bench -i human.k47.sshash >> sshash.benchmarking_log
./sshash bench -i human.k63.sshash >> sshash.benchmarking_log

### canonical indexes

./sshash bench -i celegans.k31.canon.sshash >> sshash.benchmarking_log
./sshash bench -i celegans.k47.canon.sshash >> sshash.benchmarking_log
./sshash bench -i celegans.k63.canon.sshash >> sshash.benchmarking_log

./sshash bench -i cod.k31.canon.sshash >> sshash.benchmarking_log
./sshash bench -i cod.k47.canon.sshash >> sshash.benchmarking_log
./sshash bench -i cod.k63.canon.sshash >> sshash.benchmarking_log

./sshash bench -i kestrel.k31.canon.sshash >> sshash.benchmarking_log
./sshash bench -i kestrel.k47.canon.sshash >> sshash.benchmarking_log
./sshash bench -i kestrel.k63.canon.sshash >> sshash.benchmarking_log

./sshash bench -i human.k31.canon.sshash >> sshash.benchmarking_log
./sshash bench -i human.k47.canon.sshash >> sshash.benchmarking_log
./sshash bench -i human.k63.canon.sshash >> sshash.benchmarking_log