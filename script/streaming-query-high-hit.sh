#!/bin/bash

echo "output log file =" $1

### regular indexes

./sshash query -i cod.k31.sshash -q ~/sshash_queries/SRR12858649.fastq.gz >> $1.regular.high-hit.streaming_query_log
# ./sshash query -i cod.k63.sshash -q ~/sshash_queries/SRR12858649.fastq.gz >> $1.regular.high-hit.streaming_query_log

./sshash query -i kestrel.k31.sshash -q ~/sshash_queries/SRR11449743_1.fastq.gz >> $1.regular.high-hit.streaming_query_log
# ./sshash query -i kestrel.k63.sshash -q ~/sshash_queries/SRR11449743_1.fastq.gz >> $1.regular.high-hit.streaming_query_log

./sshash query -i human.k31.sshash -q ~/sshash_queries/SRR5833294.fastq.gz >> $1.regular.high-hit.streaming_query_log
# ./sshash query -i human.k63.sshash -q ~/sshash_queries/SRR5833294.fastq.gz >> $1.regular.high-hit.streaming_query_log

### canonical indexes

./sshash query -i cod.k31.canon.sshash -q ~/sshash_queries/SRR12858649.fastq.gz >> $1.canon.high-hit.streaming_query_log
# ./sshash query -i cod.k63.canon.sshash -q ~/sshash_queries/SRR12858649.fastq.gz >> $1.canon.high-hit.streaming_query_log

./sshash query -i kestrel.k31.canon.sshash -q ~/sshash_queries/SRR11449743_1.fastq.gz >> $1.canon.high-hit.streaming_query_log
# ./sshash query -i kestrel.k63.canon.sshash -q ~/sshash_queries/SRR11449743_1.fastq.gz >> $1.canon.high-hit.streaming_query_log

./sshash query -i human.k31.canon.sshash -q ~/sshash_queries/SRR5833294.fastq.gz >> $1.canon.high-hit.streaming_query_log
# ./sshash query -i human.k63.canon.sshash -q ~/sshash_queries/SRR5833294.fastq.gz >> $1.canon.high-hit.streaming_query_log
