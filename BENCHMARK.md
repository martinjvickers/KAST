
# Benchmark hardware

44 core Intel(R) Xeon(R) CPU E5-2699A v4 @ 2.40GHz (2xsocket)
512GB RAM

# Software

* diamond version 0.9.24
* kast version: 0.0.29
* ncbi-blast+ 2.8.1

# Data

Microbial genome sequences and taxonomic information based on the Genometa 2012 data set

https://doi.org/10.20391/e6974906-f30f-4976-90fb-ea1679eedef0

## Data-prep

I needed to prep the nucl data because it contains characters that ALF doesn't processed, so I replaced those characters with an N, which is what KAST does automagically.

```
cat allgenomes_subset.list.fa | sed  '/^>/! s/H\|B\|D\|V\|W\|S\|Y\|R\|M\|K/N/g' > allgenomes_subset.list_modded.fa
```

# Results

## Pairwise

|   | KAST   |   |   |   |
|---|---|---|---|---|
|   |   |   |   |   |
|   |   |   |   |   |
|   |   |   |   |   |

## Query/Ref mode

|   | KAST  |   |   |   |
|---|---|---|---|---|
|   |   |   |   |   |
|   |   |   |   |   |
|   |   |   |   |   |

# Raw Data results

# Query/Search mode

Ref: allgenomes_april_2012_v6_one_per_species.fa.zip
Qry: allgenomes_subset.list.fa.zip

## BLAST+

#### making db

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./ncbi-blast-2.8.1+/bin/makeblastdb -in allgenomes_april_2012_v6_one_per_species.fa -dbtype nucl -out allgenomes_april_2012_v6_one_per_species


Building a new DB, current time: 07/22/2019 14:37:59
New DB name:   /users/hpctestuser/mvickers/allgenomes_april_2012_v6_one_per_species
New DB title:  allgenomes_april_2012_v6_one_per_species.fa
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 1052 sequences in 24.2921 seconds.
        Command being timed: "./ncbi-blast-2.8.1+/bin/makeblastdb -in allgenomes_april_2012_v6_one_per_species.fa -dbtype nucl -out allgenomes_april_2012_v6_one_per_species"
        User time (seconds): 23.08
        System time (seconds): 1.36
        Percent of CPU this job got: 73%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:33.21
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 299584
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 20286
        Voluntary context switches: 51075
        Involuntary context switches: 27203
        Swaps: 0
        File system inputs: 0
        File system outputs: 1738160
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4097
        Exit status: 0

```

### blast results

```
```

## Diamond

#### making db

```
Total time = 69.0017s
        Command being timed: "./diamond makedb --threads 44 --in allgenomes_april_2012_v6_one_per_species.fa -d allgenomes_april_2012_v6_one_per_species"
        User time (seconds): 911.97
        System time (seconds): 12.02
        Percent of CPU this job got: 1338%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:09.00
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 8169664
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 2684095
        Voluntary context switches: 11001
        Involuntary context switches: 6077
        Swaps: 0
        File system inputs: 0
        File system outputs: 6951416
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

Alignment step

```
        Command being timed: "./diamond blastx -d allgenomes_april_2012_v6_one_per_species -q allgenomes_subset.list.fa -o diamond_queryref.m8 --threads 44"
        User time (seconds): 2984.95
        System time (seconds): 48.13
        Percent of CPU this job got: 2167%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:19.93
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 34521568
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 8045357
        Voluntary context switches: 160869
        Involuntary context switches: 37658
        Swaps: 0
        File system inputs: 56
        File system outputs: 104
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```


## KAST

### d2 k=3

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 3 -t d2 -s dna -o results_kast_k3_d2.txt 
        Command being timed: "./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 3 -t d2 -s dna -o results_kast_k3_d2.txt"
        User time (seconds): 61.40
        System time (seconds): 8.86
        Percent of CPU this job got: 145%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:48.19
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 31144000
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3054120
        Voluntary context switches: 1868
        Involuntary context switches: 104
        Swaps: 0
        File system inputs: 0
        File system outputs: 1432
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

### d2 k=5

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 5 -t d2 -s dna -o results_kast_k5_d2.txt 
        Command being timed: "./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 5 -t d2 -s dna -o results_kast_k5_d2.txt"
        User time (seconds): 72.67
        System time (seconds): 9.00
        Percent of CPU this job got: 147%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:55.19
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 31250592
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3072186
        Voluntary context switches: 1831
        Involuntary context switches: 115
        Swaps: 0
        File system inputs: 0
        File system outputs: 1424
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0

```

### d2 k=7

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 7 -t d2 -s dna -o results_kast_k7_d2.txt
        Command being timed: "./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 7 -t d2 -s dna -o results_kast_k7_d2.txt"
        User time (seconds): 96.86
        System time (seconds): 9.43
        Percent of CPU this job got: 169%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:02.54
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 31396048
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3079209
        Voluntary context switches: 1936
        Involuntary context switches: 148
        Swaps: 0
        File system inputs: 0
        File system outputs: 1424
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0

```

### d2 k=9

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 9 -t d2 -s dna -o results_kast_k9_d2.txt
        Command being timed: "./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 9 -t d2 -s dna -o results_kast_k9_d2.txt"
        User time (seconds): 498.29
        System time (seconds): 11.78
        Percent of CPU this job got: 494%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:43.24
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 33455984
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3189834
        Voluntary context switches: 1863
        Involuntary context switches: 2645
        Swaps: 0
        File system inputs: 0
        File system outputs: 1424
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0

```

### d2 k=11

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 11 -t d2 -s dna -o results_kast_k11_d2.txt
        Command being timed: "./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 11 -t d2 -s dna -o results_kast_k11_d2.txt"
        User time (seconds): 10253.30
        System time (seconds): 32.17
        Percent of CPU this job got: 2828%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 6:03.67
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 104622160
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3813533
        Voluntary context switches: 1978
        Involuntary context switches: 40586
        Swaps: 0
        File system inputs: 0
        File system outputs: 1704
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0

```

### d2s k=3 m=0

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 3 -t d2s -s dna -o results_kast_k3_d2s.txt
        Command being timed: "./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 3 -t d2s -s dna -o results_kast_k3_d2s.txt"
        User time (seconds): 87.59
        System time (seconds): 9.15
        Percent of CPU this job got: 149%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:04.68
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 31317584
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3069038
        Voluntary context switches: 1878
        Involuntary context switches: 145
        Swaps: 0
        File system inputs: 0
        File system outputs: 1424
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

### d2s k=5 m=0

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 5 -t d2s -s dna -o results_kast_k5_d2s.txt 
        Command being timed: "./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 5 -t d2s -s dna -o results_kast_k5_d2s.txt"
        User time (seconds): 104.42
        System time (seconds): 9.15
        Percent of CPU this job got: 158%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:11.75
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 31489408
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3060907
        Voluntary context switches: 1933
        Involuntary context switches: 153
        Swaps: 0
        File system inputs: 0
        File system outputs: 1416
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

### d2s k=7 m=0

```
[hpctestuser@node08(rvchpc1) mvickers]$ /usr/bin/time -v ./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 7 -t d2s -s dna -o results_kast_k7_d2s.txt
        Command being timed: "./kast -r allgenomes_april_2012_v6_one_per_species.fa -q allgenomes_subset.list.fa -c 44 -n 10 -k 7 -t d2s -s dna -o results_kast_k7_d2s.txt"
        User time (seconds): 212.68
        System time (seconds): 9.73
        Percent of CPU this job got: 264%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:24.17
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 32118416
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3082862
        Voluntary context switches: 1963
        Involuntary context switches: 344
        Swaps: 0
        File system inputs: 0
        File system outputs: 1416
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```







# Pairwise analysis

## KAST

AminoAcid mode

### d2 k = 3 AminoAcid

```
ubuntu@vickers:~/kast_performance_test$ /usr/bin/time -v ./kast -p yeastSC.clean.fa -k 3 -t d2 -c 16 -o pw_kast_d2_k3.txt -s aa
 Command being timed: "./kast -p yeastSC.clean.fa -k 3 -t d2 -c 16 -o pw_kast_d2_k3.txt -s aa"
 User time (seconds): 942.56
 System time (seconds): 9.17
 Percent of CPU this job got: 1188%
 Elapsed (wall clock) time (h:mm:ss or m:ss): 1:20.06
 Average shared text size (kbytes): 0
 Average unshared data size (kbytes): 0
 Average stack size (kbytes): 0
 Average total size (kbytes): 0
 Maximum resident set size (kbytes): 904808
 Average resident set size (kbytes): 0
 Major (requiring I/O) page faults: 0
 Minor (reclaiming a frame) page faults: 225762
 Voluntary context switches: 315940
 Involuntary context switches: 4463
 Swaps: 0
 File system inputs: 8
 File system outputs: 775336
 Socket messages sent: 0
 Socket messages received: 0
 Signals delivered: 0
 Page size (bytes): 4096
 Exit status: 0
```

### d2 k = 3 ReducedAminoAcid


### d2 k = 5 ReducedAminoAcid

```

```
