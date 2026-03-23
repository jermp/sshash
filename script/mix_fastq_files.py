#!/usr/bin/env python3

import gzip
import random
import argparse

def fastq_reader(file_path):
    """
    Generator that yields one complete FASTQ read (4 lines) at a time
    from a gzipped file.
    """
    with gzip.open(file_path, 'rt') as f:
        while True:
            line1 = f.readline()
            if not line1:
                break # EOF
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline()
            
            # Yield the entire 4-line read as a single string
            yield line1 + line2 + line3 + line4

def mix_fastq(file1, file2, output_file):
    print(f"Mixing reads from:\n  - {file1}\n  - {file2}")
    print(f"Writing mixed reads to: {output_file} ...")
    
    # Initialize generators for both files
    iter1 = fastq_reader(file1)
    iter2 = fastq_reader(file2)
    
    # Keep track of active iterators
    active_iters = [iter1, iter2]
    
    written_count = 0
    
    with gzip.open(output_file, 'wt') as out:
        while active_iters:
            # Pick a random index from the currently active iterators
            idx = random.randrange(len(active_iters))
            
            try:
                # Fetch the next read from the randomly chosen file
                read = next(active_iters[idx])
                out.write(read)
                written_count += 1
                
                if written_count == 3_000_000:
                    break

                # Optional: print progress every 1M reads
                if written_count % 1_000_000 == 0:
                    print(f"Processed {written_count:,} reads...")
                    
            except StopIteration:
                # If the chosen file is exhausted, remove it from the active pool
                active_iters.pop(idx)
                
    print(f"Done! Successfully mixed {written_count:,} total reads into {output_file}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mix two gzipped FASTQ files with uniform probability.")
    parser.add_argument("-1", "--file1", required=True, help="Path to the first compressed fastq file (.fastq.gz)")
    parser.add_argument("-2", "--file2", required=True, help="Path to the second compressed fastq file (.fastq.gz)")
    parser.add_argument("-o", "--output", required=True, help="Path to the output mixed fastq file (.fastq.gz)")
    
    args = parser.parse_args()
    
    mix_fastq(args.file1, args.file2, args.output)
