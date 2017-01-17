/*
MIT License

Copyright (c) 2017 Martin James Vickers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include "boost/multi_array.hpp"
#include <boost/unordered_map.hpp>
#include <cassert>
#include <string>
#include <thread>
#include <mutex>
#include <unordered_map>
using namespace seqan;
using namespace std;

/*
For a given kmer when calculating d2s and d2star you'll store 
the counts of that kmer in the sequence as well as the markov
probability. markov_dat is this data structure.
*/
struct markov_dat
{
        long long int count;
        double prob;
};

struct count_obj
{
	unordered_map<string, long long int> kmer_counts;
	int total;
};

struct markov_obj
{
        unordered_map<string, markov_dat> markov_counts;
        int total_count;
	double sum_prob;
};

/*
User defined options struct
*/
struct ModifyStringOptions
{
        unsigned klen;
        int nohits;
        int markovOrder;
        CharString type;
	CharString output_format;
        bool noreverse;
        CharString queryFileName = NULL;
        CharString referenceFileName = NULL;
	CharString pairwiseFileName = NULL;
        int num_threads;
        bool debug;
        bool useram;
	CharString outputFileName = NULL;
};

