#!/bin/bash

mkdir DNA_datasets

wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -O DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
~/bcalm/build/bcalm -in ~/DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 64
~/UST/ust -k 31 -i ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa
gzip Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa.ust.fa
mv Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa.ust.fa.gz DNA_datasets/

wget http://ftp.ensembl.org/pub/current_fasta/gadus_morhua/dna/Gadus_morhua.gadMor3.0.dna.toplevel.fa.gz -O DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa.gz
~/bcalm/build/bcalm -in ~/DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 64
~/UST/ust -k 31 -i ~/Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs.fa
gzip Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs.fa.ust.fa
mv Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs.fa.ust.fa.gz DNA_datasets/

wget http://ftp.ensembl.org/pub/current_fasta/falco_tinnunculus/dna/Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.gz -O DNA_datasets/Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.gz
~/bcalm/build/bcalm -in ~/DNA_datasets/Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 64
~/UST/ust -k 31 -i ~/Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs.fa
gzip Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs.fa.ust.fa
mv Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs.fa.ust.fa.gz DNA_datasets/

wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -O DNA_datasets/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
~/bcalm/build/bcalm -in ~/DNA_datasets/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 64
~/UST/ust -k 31 -i ~/Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs.fa
gzip Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs.fa.ust.fa
mv Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs.fa.ust.fa.gz DNA_datasets/

wget https://zenodo.org/record/995689/files/bacterial.genome.fixed.fa -O DNA_datasets/bacterial.genome.fixed.fa
~/bcalm/build/bcalm -in ~/DNA_datasets/bacterial.genome.fixed.fa -kmer-size 31 -abundance-min 1 -nb-cores 64
~/UST/ust -k 31 -i ~/bacterial.genome.fixed.fa
gzip bacterial.genome.fixed.fa.unitigs.fa.ust.fa
mv bacterial.genome.fixed.fa.unitigs.fa.ust.fa.gz DNA_datasets/
