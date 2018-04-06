#ifndef UTILS_H
#define UTILS_H

#include "common.h"

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv);
AminoAcid getRevCompl(AminoAcid const & nucleotide);
String<AminoAcid> doRevCompl(String<AminoAcid> seq);

/*
   Generic function to return kmer counts
*/
template <typename TAlphabet>
map<string, unsigned int> count(String<TAlphabet> const & seq, int klen)
{
   map<string, unsigned int> map;
   for(unsigned i = 0; i <= length(seq) - klen; i++)
   {
      string kmer;
      assign(kmer,infix(seq, i, i+klen));
      size_t found = kmer.find("N");
      if(found > kmer.size())
      {
         unsigned int count = map[kmer];
         map[kmer] = count + 1;
      }
   }
   return map;
};

/*
   Generic function to return kmer counts but with a mask vector
*/
template <typename TAlphabet>
map<string, unsigned int> count(String<TAlphabet> sequence, int klen, 
                                vector<CharString> mask)
{
   map<string, unsigned int> map;
   for(int i = 0; i <= length(sequence)-klen; i++)
   {
      string kmer;
      assign(kmer,infix(sequence, i, i+klen));

      // so kmer has the wider bit we care about
      // now go through the mask and append kmers
      for(auto m : mask)
      {
         CharString masked_kmer;
         for(int loc = 0; loc < length(m); loc++)
         {
            if(m[loc] == '1')
               append(masked_kmer,kmer[loc]);
         }
                        
         unsigned int count = map[toCString(masked_kmer)];
         map[toCString(masked_kmer)] = count + 1;
      }
   }
   return map;
};

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
int printResult(ModifyStringOptions options, CharString &queryid,
                ofstream &outfile, String<AminoAcid> &queryseq,
                map<double, int> &results, StringSet<CharString> &referenceids,
                StringSet<String<AminoAcid>> &referenceseqs);
#endif
