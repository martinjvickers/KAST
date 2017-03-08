###AWS instance CPU
```
processor       : 15
vendor_id       : GenuineIntel
cpu family      : 6
model           : 79
model name      : Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz
stepping        : 1
microcode       : 0xb000014
cpu MHz         : 2698.816
cache size      : 46080 KB
physical id     : 0
siblings        : 16
core id         : 7
cpu cores       : 8
apicid          : 15
initial apicid  : 15
fpu             : yes
fpu_exception   : yes
cpuid level     : 13
wp              : yes
flags           : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ht syscall nx pdpe1gb rdtscp lm constant_tsc rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq ssse3 fma cx16 pcid sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand hypervisor lahf_lm abm 3dnowprefetch fsgsbase bmi1 hle avx2 smep bmi2 erms invpcid rtm rdseed adx xsaveopt
bugs            :
bogomips        : 4600.11
clflush size    : 64
cache_alignment : 64
address sizes   : 46 bits physical, 48 bits virtual
power management:
```


###Individual timings

```
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u"
	User time (seconds): 188.78
	System time (seconds): 182.80
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 6:11.98
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 137200
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135358561
	Voluntary context switches: 65
	Involuntary context switches: 533
	Swaps: 0
	File system inputs: 0
	File system outputs: 4152
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 16
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 16"
	User time (seconds): 355.04
	System time (seconds): 296.08
	Percent of CPU this job got: 722%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:30.15
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 693496
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135311284
	Voluntary context switches: 15708
	Involuntary context switches: 1175
	Swaps: 0
	File system inputs: 0
	File system outputs: 4072
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 8
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 8"
	User time (seconds): 307.18
	System time (seconds): 260.76
	Percent of CPU this job got: 686%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:22.67
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 383980
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135234829
	Voluntary context switches: 10296
	Involuntary context switches: 646
	Swaps: 0
	File system inputs: 0
	File system outputs: 4072
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 12
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 12"
	User time (seconds): 351.02
	System time (seconds): 293.94
	Percent of CPU this job got: 722%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:29.28
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 515320
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135297243
	Voluntary context switches: 13199
	Involuntary context switches: 994
	Swaps: 0
	File system inputs: 0
	File system outputs: 4072
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 4
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 4"
	User time (seconds): 236.13
	System time (seconds): 228.33
	Percent of CPU this job got: 375%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:03.70
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 243456
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135234624
	Voluntary context switches: 3534
	Involuntary context switches: 540
	Swaps: 0
	File system inputs: 0
	File system outputs: 4088
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 2
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 2"
	User time (seconds): 212.45
	System time (seconds): 205.55
	Percent of CPU this job got: 194%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:34.82
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 176284
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135296787
	Voluntary context switches: 2037
	Involuntary context switches: 531
	Swaps: 0
	File system inputs: 0
	File system outputs: 4104
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ uptime
 12:29:15 up 27 min,  2 users,  load average: 1.03, 1.85, 1.65
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 6
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 6"
	User time (seconds): 286.91
	System time (seconds): 258.74
	Percent of CPU this job got: 539%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:41.16
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 317880
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135277441
	Voluntary context switches: 6354
	Involuntary context switches: 603
	Swaps: 0
	File system inputs: 0
	File system outputs: 4080
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 10
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 10"
	User time (seconds): 333.58
	System time (seconds): 283.28
	Percent of CPU this job got: 698%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:28.28
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 450468
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135251608
	Voluntary context switches: 12467
	Involuntary context switches: 789
	Swaps: 0
	File system inputs: 0
	File system outputs: 4072
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 14
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 14"
	User time (seconds): 361.61
	System time (seconds): 304.76
	Percent of CPU this job got: 739%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:30.16
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 585800
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135238959
	Voluntary context switches: 14238
	Involuntary context switches: 1078
	Swaps: 0
	File system inputs: 0
	File system outputs: 4080
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u"
	User time (seconds): 188.78
	System time (seconds): 182.80
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 6:11.98
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 137200
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135358561
	Voluntary context switches: 65
	Involuntary context switches: 533
	Swaps: 0
	File system inputs: 0
	File system outputs: 4152
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 16
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 16"
	User time (seconds): 355.04
	System time (seconds): 296.08
	Percent of CPU this job got: 722%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:30.15
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 693496
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135311284
	Voluntary context switches: 15708
	Involuntary context switches: 1175
	Swaps: 0
	File system inputs: 0
	File system outputs: 4072
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 8
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 8"
	User time (seconds): 307.18
	System time (seconds): 260.76
	Percent of CPU this job got: 686%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:22.67
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 383980
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135234829
	Voluntary context switches: 10296
	Involuntary context switches: 646
	Swaps: 0
	File system inputs: 0
	File system outputs: 4072
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 12
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 12"
	User time (seconds): 351.02
	System time (seconds): 293.94
	Percent of CPU this job got: 722%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:29.28
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 515320
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135297243
	Voluntary context switches: 13199
	Involuntary context switches: 994
	Swaps: 0
	File system inputs: 0
	File system outputs: 4072
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 4
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 4"
	User time (seconds): 236.13
	System time (seconds): 228.33
	Percent of CPU this job got: 375%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:03.70
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 243456
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135234624
	Voluntary context switches: 3534
	Involuntary context switches: 540
	Swaps: 0
	File system inputs: 0
	File system outputs: 4088
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 2
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 2"
	User time (seconds): 212.45
	System time (seconds): 205.55
	Percent of CPU this job got: 194%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:34.82
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 176284
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135296787
	Voluntary context switches: 2037
	Involuntary context switches: 531
	Swaps: 0
	File system inputs: 0
	File system outputs: 4104
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ uptime
 12:29:15 up 27 min,  2 users,  load average: 1.03, 1.85, 1.65
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 6
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 6"
	User time (seconds): 286.91
	System time (seconds): 258.74
	Percent of CPU this job got: 539%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:41.16
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 317880
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135277441
	Voluntary context switches: 6354
	Involuntary context switches: 603
	Swaps: 0
	File system inputs: 0
	File system outputs: 4080
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 10
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 10"
	User time (seconds): 333.58
	System time (seconds): 283.28
	Percent of CPU this job got: 698%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:28.28
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 450468
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135251608
	Voluntary context switches: 12467
	Involuntary context switches: 789
	Swaps: 0
	File system inputs: 0
	File system outputs: 4072
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
[ec2-user@ip-172-31-17-61 ~]$ /usr/bin/time -v ./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 14
	Command being timed: "./alfsc -q SRR042642_10k.fastq -r yeast.fasta -o meh.txt -u -c 14"
	User time (seconds): 361.61
	System time (seconds): 304.76
	Percent of CPU this job got: 739%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:30.16
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 585800
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 135238959
	Voluntary context switches: 14238
	Involuntary context switches: 1078
	Swaps: 0
	File system inputs: 0
	File system outputs: 4080
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
