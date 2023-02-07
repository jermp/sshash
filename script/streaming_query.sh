#!/bin/bash

### regular indexes

./sshash query -i celegans.k31.sshash -q /data2/DNA/queries/SRR5901135.fastq.gz >> sshash.streaming_query_log
./sshash query -i celegans.k47.sshash -q /data2/DNA/queries/SRR5901135.fastq.gz >> sshash.streaming_query_log
./sshash query -i celegans.k63.sshash -q /data2/DNA/queries/SRR5901135.fastq.gz >> sshash.streaming_query_log

./sshash query -i cod.k31.sshash -q /data2/DNA/queries/SRR12858649.fastq.gz >> sshash.streaming_query_log
./sshash query -i cod.k47.sshash -q /data2/DNA/queries/SRR12858649.fastq.gz >> sshash.streaming_query_log
./sshash query -i cod.k63.sshash -q /data2/DNA/queries/SRR12858649.fastq.gz >> sshash.streaming_query_log

./sshash query -i kestrel.k31.sshash -q /data2/DNA/queries/SRR11449743.fastq.gz >> sshash.streaming_query_log
./sshash query -i kestrel.k47.sshash -q /data2/DNA/queries/SRR11449743.fastq.gz >> sshash.streaming_query_log
./sshash query -i kestrel.k63.sshash -q /data2/DNA/queries/SRR11449743.fastq.gz >> sshash.streaming_query_log

./sshash query -i human.k31.sshash -q /data2/DNA/queries/SRR5833294.fastq.gz >> sshash.streaming_query_log
./sshash query -i human.k47.sshash -q /data2/DNA/queries/SRR5833294.fastq.gz >> sshash.streaming_query_log
./sshash query -i human.k63.sshash -q /data2/DNA/queries/SRR5833294.fastq.gz >> sshash.streaming_query_log

### canonical indexes

./sshash query -i celegans.k31.canon.sshash -q /data2/DNA/queries/SRR5901135.fastq.gz >> sshash.streaming_query_log
./sshash query -i celegans.k47.canon.sshash -q /data2/DNA/queries/SRR5901135.fastq.gz >> sshash.streaming_query_log
./sshash query -i celegans.k63.canon.sshash -q /data2/DNA/queries/SRR5901135.fastq.gz >> sshash.streaming_query_log

./sshash query -i cod.k31.canon.sshash -q /data2/DNA/queries/SRR12858649.fastq.gz >> sshash.streaming_query_log
./sshash query -i cod.k47.canon.sshash -q /data2/DNA/queries/SRR12858649.fastq.gz >> sshash.streaming_query_log
./sshash query -i cod.k63.canon.sshash -q /data2/DNA/queries/SRR12858649.fastq.gz >> sshash.streaming_query_log

./sshash query -i kestrel.k31.canon.sshash -q /data2/DNA/queries/SRR11449743.fastq.gz >> sshash.streaming_query_log
./sshash query -i kestrel.k47.canon.sshash -q /data2/DNA/queries/SRR11449743.fastq.gz >> sshash.streaming_query_log
./sshash query -i kestrel.k63.canon.sshash -q /data2/DNA/queries/SRR11449743.fastq.gz >> sshash.streaming_query_log

./sshash query -i human.k31.canon.sshash -q /data2/DNA/queries/SRR5833294.fastq.gz >> sshash.streaming_query_log
./sshash query -i human.k47.canon.sshash -q /data2/DNA/queries/SRR5833294.fastq.gz >> sshash.streaming_query_log
./sshash query -i human.k63.canon.sshash -q /data2/DNA/queries/SRR5833294.fastq.gz >> sshash.streaming_query_log
