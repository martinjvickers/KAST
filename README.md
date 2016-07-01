### Quick summary ###
This is alfsc - Alignment-free Sequence Comparsion (until maybe we can think of a better name).

### Version ###
Early Early Beta! 0.0.1 :-)

### How do I get set up on Ubuntu 16.04? ###

* sudo apt-get install git g++ build-essential cmake zlib1g-dev libbz2-dev libboost-all-dev

* git clone https://github.com/seqan/seqan.git seqan

* git clone https://martinjvickers@bitbucket.org/martinjvickers/alfsc-stable.git

* mkdir -p alfsc-build/release
* cd alfsc-build/release
* cmake ../../alfsc-stable -DCMAKE_MODULE_PATH=~/seqan/util/cmake -DSEQAN_INCLUDE_PATH=~/seqan/include -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_BUILD_TYPE=Release
* make

### Who do I talk to? ###

* Martin Vickers (martinj.vickers@gmail.com)
