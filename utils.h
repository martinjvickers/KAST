#ifndef UTILS_H
#define UTILS_H

#include "common.h"

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv);
AminoAcid getRevCompl(AminoAcid const & nucleotide);
String<AminoAcid> doRevCompl(String<AminoAcid> seq);
map<string, unsigned int> count(String<AminoAcid> sequence, int klen);
map<string, unsigned int> count(String<AminoAcid> sequence, int klen, 
                                vector<CharString> mask);
double gc_ratio(String<AminoAcid> sequence);
map<string, double> markov(int klen, String<AminoAcid> sequence, 
                           int markovOrder, map<string, bool> kmer_count_map);
map<string, bool> makecomplete(ModifyStringOptions options);
map<string, bool> makequick(ModifyStringOptions options, 
                            StringSet<String<AminoAcid>> referenceseqs);
int ipow(int base, int exp);
CharString namecut(CharString seq, int val);
int parseMask(ModifyStringOptions options, int &effectiveKlen);
int printPhylyp(ModifyStringOptions options,    
                vector<pair<CharString, map<string, unsigned int>>> pw_counts,
                vector< vector<double>> &array_threaded_internal);

#endif
