### Quick summary ###
This is a rewrite

### Version ###
Early Early Beta! 0.0.3 :-)

### How do I get set up on Ubuntu 16.04? ###

```
sudo apt-get install git g++ build-essential cmake zlib1g-dev libbz2-dev libboost-all-dev
git clone https://github.com/seqan/seqan.git seqan
git clone https://github.com/martinjvickers/alfsc.git
cd alfsc
cmake ../alfsc -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
make
```

### Quick run test ###

```
/usr/bin/time -v ./alfsc -q example_data/SRR042642_1.fastq.gz -r example_data/yeast.fasta -k 3 -c 4
```

### Binary Release ###

Using standard CMAKE static build options a binary release has been created that works on 64bit x86 GNU/Linux machines, however it is not fully backwards compatible for old kernels since alfsc uses c++14 and SeqAn which make creating a backward compatible binary for very old kernels difficult. This does not mean that you will not be able to use alfsc, it just means that you will have to build alfsc from source. 

The current beta binary release has been tested on CentOS6 and 7 (which implies that it should work on SL6 and 7), Ubuntu 12.04.5 and 16.04.1. 

### Who do I talk to? ###

* Martin Vickers (martinj.vickers@gmail.com)
