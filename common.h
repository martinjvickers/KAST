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

struct thingy
{
        long long int count;
        double prob;
};
