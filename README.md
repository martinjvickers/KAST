[![Build Status](https://travis-ci.org/martinjvickers/KAST.svg?branch=master)](https://travis-ci.org/martinjvickers/KAST)

# Summary #

Perform Alignment-free k-tuple frequency comparisons from sequences. This can be in the form of two input files (e.g. a reference and a query) or a single file for pairwise comparisons to be made.

## Documentation

There are several features of KAST documented in the github wiki. 

https://github.com/martinjvickers/KAST/wiki

## Installation

KAST can be installed in a variety of different ways depending on your environment. For detailed installation documentation, take a look at the wiki.

https://github.com/martinjvickers/KAST/wiki/Installation

## Quick install and start ###

If you have a modern Linux OS, you should be able to use the static binary.

Get static binary and extract tar.gz

```
wget https://github.com/martinjvickers/KAST/releases/download/1.0.0/KAST_v1.0.0.tar.gz
tar xvfz KAST_v1.0.0.tar.gz
```

Get some example data to run;

```
wget -c https://github.com/martinjvickers/KAST/blob/master/example_data/yeast.fasta?raw=true -O yeast.fasta
wget -c https://github.com/martinjvickers/KAST/blob/master/example_data/SRR042642_100k.fastq.gz?raw=true -O SRR042642_100k.fastq.gz
```

And now try the program (`-c` assuming 4 CPU cores)

```
/usr/bin/time -v ./kast -q SRR042642_100k.fastq.gz -r yeast.fasta -k 3 -c 4 -o output.txt
```

If you wish to do pairwise comparison the following command;

```
./kast -p example_data/yeast.fasta -o test.txt
```

### Notes on the binary release ###

Using standard CMAKE static build options a binary release has been created that works on 64bit x86 GNU/Linux machines, however it is not fully backwards compatible for old kernels since alfsc uses c++14 and SeqAn which make creating a backward compatible binary for very old kernels difficult. This does not mean that you will not be able to use alfsc, it just means that you will have to build alfsc from source or use a container. 

The current binary release has been built on Ubuntu LTS 20.04. If there is sufficient demand I can make older static binaries, otherwise one can use the installation on biocontainers.

### Who do I talk to? ###

* Martin Vickers (martinj.vickers@gmail.com)
