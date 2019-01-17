[![Build Status](https://travis-ci.org/martinjvickers/KAST.svg?branch=master)](https://travis-ci.org/martinjvickers/KAST)

### Quick summary ###
Perform Alignment-free k-tuple frequency comparisons from sequences. This can be in the form of two input files (e.g. a reference and a query) or a single file for pairwise comparisons to be made.

### Version ###
Not far from a complete release but still a bit to go. This is version 0.0.24.

### Manual and Usage ###

https://github.com/martinjvickers/KAST/wiki

### How do I get set up on Ubuntu 16.04? ###

```
sudo apt-get install git g++ build-essential cmake zlib1g-dev libbz2-dev libboost-all-dev
git clone https://github.com/seqan/seqan.git seqan
git clone https://github.com/martinjvickers/KAST.git
cd KAST
cmake ../KAST -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
make
```

### Quick run test ###

```
/usr/bin/time -v ./kast -q example_data/SRR042642_1.fastq.gz -r example_data/yeast.fasta -k 3 -c 4 -o output.txt
```

If you wish to do pairwise comparison the following command;

```
./kast -p example_data/yeast.fasta -o test.txt
```

### Binary Release ###

Using standard CMAKE static build options a binary release has been created that works on 64bit x86 GNU/Linux machines, however it is not fully backwards compatible for old kernels since alfsc uses c++14 and SeqAn which make creating a backward compatible binary for very old kernels difficult. This does not mean that you will not be able to use alfsc, it just means that you will have to build alfsc from source. 

The current beta binary release has been tested on CentOS6 and 7 (which implies that it should work on SL6 and 7), Ubuntu 12.04.5 and 16.04.1. 

### Who do I talk to? ###

* Martin Vickers (martinj.vickers@gmail.com)
