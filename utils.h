#ifndef UTILS_H
#define UTILS_H

/*
#include <seqan/seq_io.h>
#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>      
#include <seqan/store.h>
#include <string>
#include <thread>
#include <map>
#include <vector>
*/

#include "common.h"

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv);
AminoAcid getRevCompl(AminoAcid const & nucleotide);
Dna5String doRevCompl(Dna5String seq);
map<string, unsigned int> count(String<AminoAcid> sequence, int klen);
double gc_ratio(String<AminoAcid> sequence);
map<string, double> markov(int klen, String<AminoAcid> sequence, int markovOrder, map<string, bool> kmer_count_map);
map<string, bool> makecomplete(ModifyStringOptions options);
map<string, bool> makequick(ModifyStringOptions options, StringSet<Dna5String> referenceseqs);
//int makequick(ModifyStringOptions options, StringSet<Dna5String> referenceseqs, map<string, bool> quickkmers);
int ipow(int base, int exp);
CharString namecut(CharString seq, int val);
#endif
