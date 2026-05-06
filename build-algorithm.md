# SSHash build algorithm

This note describes how `sshash build` constructs a dictionary while keeping
peak resident memory bounded by the user-supplied `-g` (in GiB).

The design has two ideas, applied uniformly:

1. **Spill, don't accumulate.** Every intermediate that grows with the input
   size is written to a tmp file under `-d` (tmp dir) rather than held in a
   `std::vector` / bit-vector in RAM. Producers append through a small write
   buffer; consumers re-read through a small read buffer
   (`buffered_record_stream<T>`).
2. **Cap working buffers at a fraction of `-g`.** Buffers that live
   only inside one step are sized as `ram_limit_in_GiB · GiB / N` (with `N`
   typically 2 or 8). The constants are picked so that even when several
   buffers are alive at the same time across overlapping steps, their sum
   stays under the user budget while heap fragmentation across step
   transitions is absorbed.

The build never materializes the final index in RAM. Instead, step 8
streams it directly to the user-supplied output file (or a tmp file, deleted
on exit, when the user did not pass `-o`). To run `--check`, `tools/build.cpp`
mmaps the saved file and runs the correctness queries against the mmap'd
dictionary.

---

## Pipeline overview

The orchestration is in `include/builder/dictionary_builder.hpp`
(`run_steps_1_through_7` + `build`). Per-step details are in
`src/builder/{encode_strings,compute_minimizer_tuples,build_sparse_and_skew_index}.cpp`.

| Step | What it produces | Where it lives between steps |
|------|------------------|------------------------------|
| 1    | Encoded `strings` bit-vector + `strings_offsets` | tmp files (`disk_backed_strings`, `disk_backed_offsets_builder`) |
| 1.1  | Compressed weights (only if `--weighted`) | `weights::builder` (in-RAM, bounded by run-length structure) |
| 2    | Per-thread sorted runs of `minimizer_tuple` | tmp files, one per flushed buffer |
| 3    | Single sorted run of all `minimizer_tuple`s | tmp file (k-way external merge) |
| 4    | Minimizers MPHF F | tmp file (spilled at end of step 5) |
| 5    | Minimizer values replaced by F(minimizer); buffers re-flushed in F-order | new sorted runs, tmp files |
| 6    | Single sorted run keyed by F(minimizer) | tmp file |
| 7.1  | Sparse-index components (`control_codewords`, `mid_load_buckets`) | tmp files |
| 7.2  | Skew-index components (`heavy_load_buckets`, per-partition MPHFs and `positions`) | tmp files |
| 8    | Final on-disk index file | streamed to output, tmp files removed |

After step 8 the dictionary object `d` is **not** query-ready: the spilled
components were copied into the output file but never read back into `d`.
`finalize_stats` reports `index_size_in_bytes` via `std::filesystem::file_size`
on the saved path.

---

## Step 1 — encode strings (`encode_strings.cpp`)

Iterates the input FASTA, producing the 2-bit-packed `strings` bit-vector
and the `strings_offsets` array (one offset per sequence + a sentinel).
Both go through disk-backed builders:

- **`disk_backed_strings`**: appends 2-bit characters into a small in-RAM
  word buffer; flushes the buffer to a tmp file when full.
- **`disk_backed_offsets_builder<Offsets>`**: appends one `uint64_t` offset
  per sequence into a small write buffer; flushes to a tmp file.

In-RAM footprint of step 1 is `O(buffer)` regardless of input size.

## Step 1.1 — weights (optional)

Only runs with `--weighted`. The weights builder uses run-length encoding:
its in-RAM size is proportional to the number of distinct weights, not to
the number of k-mers.

## Step 2 — compute minimizer tuples (`compute_minimizer_tuples.cpp`)

Each thread streams its assigned slice of the input via the disk-backed
strings/offsets readers and emits `minimizer_tuple` records into a private
in-RAM buffer:

```cpp
buffer_size = (ram_limit · GiB) / (2 · sizeof(minimizer_tuple) · num_threads)
```

When the buffer fills, the thread sorts it in parallel and flushes a sorted
run to a tmp file (`minimizers_tuples::sort_and_flush`). The factor of 2
in the denominator leaves headroom for `std::sort`'s allocations and
inter-thread contention; the per-thread split makes the total in-RAM tuple
buffer ≈ `ram_limit / 2`.

## Step 3 — k-way external merge (`minimizers_tuples::merge`)

The N tmp files from step 2 are merged into a single sorted run via a
**winner-tree-based external-merge iterator** (`file_merging_iterator<T>`).
Each input file is read through a `buffered_record_stream<minimizer_tuple>`
with `default_buffer_records = 4096` records, so the total in-RAM merge
state is `N · 4096 · sizeof(minimizer_tuple)` ≈ tens of MB even for very
many runs. The output is written through a small `std::ofstream` buffer.

When N == 1 the merge degenerates to a rename + a single streaming scan to
collect bucket statistics; same RAM bound.

## Step 4 — build minimizers MPHF

Builds an external-memory partitioned PHF over distinct minimizers, using
pthash's `build_in_external_memory`. The minimizers are streamed from the
sorted run via `streaming_minimizers_iterator` (one buffered ifstream),
and pthash spills its own working hashes under `tmp_dirname` capped by
`mphf_build_config.ram = ram_limit / 2`.

## Step 5 — replace minimizer values with F(minimizer)

The merged file is re-read in fixed-size blocks; each block is hashed in
parallel and re-flushed as a new sorted run. Two RAM caps are combined:

```cpp
RAM_available    = ram_limit · GiB − sizeof(F) − offsets_builder.num_bytes()
buffer_unbounded = RAM_available / (3 · sizeof(minimizer_tuple))   // 3× = read+sort scratch+write
buffer_cap       = (ram_limit · GiB / 8) / sizeof(minimizer_tuple)
buffer_size      = min(buffer_unbounded, buffer_cap)
```

The `/ 8` cap exists because step 5 leaves heap pages dirtied that linger
into later steps' allocations; capping at one-eighth of the budget keeps
the cumulative RSS under `ram_limit` when steps 6/7 start allocating.

After step 5, the minimizers MPHF F is **spilled to disk** and the in-RAM
copy is freed: subsequent steps only ever use F(minimizer) values, not F
itself.

## Step 6 — re-merge in F-order

Same machinery as step 3, applied to the new sorted runs from step 5.

## Step 7.1 — sparse index (`build_sparse_and_skew_index.cpp`)

Constructs `control_codewords` and `mid_load_buckets`. Both are produced as
on-disk `bits::compact_vector` files via `streaming_compact_vector_writer`,
so neither is ever materialized in RAM.

## Step 7.2 — skew index

The most RAM-sensitive step; it has three internal phases, all
disk-backed:

- **Phase B (k-mer extraction requests).** Heavy-bucket entries become
  `kmer_extraction_request` records. They are external-sorted by
  `starting_pos` so that k-mer extraction reduces to a single forward
  scan over `strings`. The request buffer is capped at
  `ram_limit / 8 / sizeof(kmer_extraction_request)`; flushed runs are
  merged with `file_merging_iterator`.
- **Per-partition kmer files.** While walking `strings` in request-sorted
  order, each extracted k-mer is written to its partition's tmp file via
  a buffered writer; this file is the input to the partition's MPHF.
- **Phase C (per-partition MPHF + `positions`).** For each skew partition:
  1. Build the partition MPHF with pthash external-memory (`ram = ram_limit / 2`,
     iterator: `skew_partition_kmer_iterator` over the partition's tmp file).
  2. Stream-read the partition file again, emit `(F(kmer), pos_in_bucket)`
     tuples; external-sort them in `ram_limit / 8`-sized buffers and merge.
  3. Pack the sorted tuples into the partition's `positions`
     compact_vector via `streaming_compact_vector_writer`.

  Only the freshly-built MPHF for the *current* partition lives in RAM
  during phase C; once spilled (`essentials::save`), it is freed before the
  next partition starts. `positions` is fully on-disk.

## Step 8 — stream-save (`include/builder/streaming_save.hpp`)

The dictionary `d` is walked by `essentials::saver`, but every spilled
component is intercepted via an **address+type-keyed substitution map**
(`typed_address_sub`): when the saver visits a registered (address, type)
pair, it appends the bytes of the corresponding tmp file straight into the
output stream instead of reading from `d`. The strings bit-vector goes
through the same mechanism via `disk_backed_strings`.

Concretely, the registered substitutions are:

| Component                         | Source tmp file                        |
|-----------------------------------|----------------------------------------|
| `m_ssi.codewords.control_codewords` | step 7.1                                |
| `m_ssi.mid_load_buckets`            | step 7.1                                |
| `m_ssi.ski.heavy_load_buckets`      | step 7.2 phase B                        |
| `m_ssi.codewords.mphf`              | step 5 spill                            |
| `m_ssi.ski.positions[i]`            | step 7.2 phase C, per partition         |
| `m_ssi.ski.mphfs[i]`                | step 7.2 phase C, per partition         |
| `m_spss.strings`                    | step 1 (`disk_backed_strings`)          |

Because the substitutions are by `(address, type)` pair, a struct's address
coinciding with its first member's address does not cause confusion.

After step 8 returns, the tmp files are removed and `finalize_stats` reads
the saved file's size with `std::filesystem::file_size`.

---

## How the RAM cap is enforced — summary

The on-disk index size grows with `num_kmers`. The build's **resident**
memory does not, because every component that scales with input size is
either:

- **always on disk** (`strings`, `strings_offsets`, all sorted minimizer
  runs, the merged minimizers file, the sparse-index compact_vectors, the
  skew-index per-partition kmer/positions files, the codewords MPHF and
  per-partition MPHFs), or
- **bounded by a working buffer** sized as a fraction of `ram_limit`:

  | Buffer                                 | Cap                      |
  |----------------------------------------|--------------------------|
  | Step 2 per-thread minimizer buffer     | `ram_limit / 2 / num_threads` |
  | Step 5 hashing buffer                  | `min(ram/8, RAM_available/3)` |
  | Step 7.2 kmer-request external sort    | `ram_limit / 8`          |
  | Step 7.2 phase-C `position_tuple` sort | `ram_limit / 8`          |
  | pthash external-memory builds          | `ram_limit / 2` (its own `ram` field) |
  | Every disk-backed reader/writer        | `default_buffer_records ≈ 32 KiB` |
  | Every external merge front (per run)   | `4096 · sizeof(T)`       |

The fractions (`/2` for the dominant per-step buffer, `/8` for buffers
that span step boundaries) are chosen so that overlapping allocations and
heap fragmentation between steps stay under `ram_limit` in practice. There
is a hard floor of `min_ram_limit_in_GiB` (enforced in
`validate_and_normalize_build_config`) below which step 4's MPHF builder
no longer has enough room to make progress.

The result: peak RSS during the build is governed by `-g`, not by
the input size or by the on-disk index size, and the saved index is
identical (byte-for-byte) to one written by an in-RAM builder followed by
`essentials::save`.
