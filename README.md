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

### TODO ###

* Create a decent method to store the reference count/probability data in memory. 
** I want this to play nice with memory, e.g. either a --mem flag, or maybe it loads things into memory and if it is too much for the machine, reverts to a low memory state?

* Don't create all of the kmers, just store the ones that are actually there. This is very important for very large kmer sizes e.g. 4^9 isn't that unreasonable but is massive.
** In order to do this, I need to reimplement all the calculations and distances to be okay with missing data, e.g. reverting to zero if there is no count from say a hashtable
** The order also needs to be preserved for markov doesn't it?
** If I can manage this, it will reduce the amount of memory used.

* When no files are entered, it segfaults. I need so solve that.

* Reimplement Matrix and pairwise

* Add write to file functionality

* Write unit and app tests.

### Who do I talk to? ###

* Martin Vickers (martinj.vickers@gmail.com)
