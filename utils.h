#ifndef UTILS_H
#define UTILS_H

#include "common.h"

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <unordered_map>

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv);
AminoAcid getRevCompl(AminoAcid const & nucleotide);
String<AminoAcid> doRevCompl(String<AminoAcid> seq);

/*
   Generic function to return kmer counts
map<string, unsigned int> count(String<TAlphabet> const & seq, int klen)

   // I'm hoping to phase this out


*/
template <typename TAlphabet>
map<string, unsigned int> count(String<TAlphabet> seq, int klen)
{
   map<string, unsigned int> counts;
   if(length(seq) >= klen)
   {
      for(unsigned i = 0; i <= length(seq) - klen; i++)
      {
         string kmer;
         assign(kmer,infix(seq, i, i+klen));
         size_t found = kmer.find("N");
         if(found > kmer.size())
         {
            unsigned int count = counts[kmer];
            counts[kmer] = count + 1;
         }
      }
   }

   return counts;
};

/*
This is what I'm aiming to get across the board now.
*/
//template <typename AminoAcid>
map<String<AminoAcid>, unsigned int> count_test(String<AminoAcid> seq, int klen, bool noreverse);

//template <typename Dna5>
map<String<Dna5>, unsigned int> count_test(String<Dna5> seq, int klen, bool noreverse);

/*
   Generic function to return kmer counts but with a mask vector
*/
template <typename TAlphabet>
map<string, unsigned int> count(String<TAlphabet> sequence, int klen, 
                                vector<CharString> mask)
{
   map<string, unsigned int> map;

   if(length(sequence) >= klen)
   {
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
   }
   return map;
};

double gc_ratio(String<AminoAcid> sequence);
map<string, double> markov(int klen, String<AminoAcid> sequence, 
                           int markovOrder, map<string, bool> kmer_count_map);


// new template markov function
template <typename TAlphabet>
map<String<TAlphabet>, double> markov_test(unsigned int klen, String<TAlphabet> sequence,
                                           unsigned int markovOrder, 
                                           vector<String<TAlphabet>> kmer_count_map,
                                           bool noreverse)
{
   map<String<TAlphabet>, double> markovmap;
   double sum_prob = 0.0;

   map<String<TAlphabet>, unsigned int> markovcounts = count_test(sequence, markovOrder, noreverse);
   double tot = 0;

   // calculate sum of counts
   for(auto p: markovcounts)
   {
      tot = tot + p.second;
   }

   for(auto kmer : kmer_count_map)
   {
      double prob = 1.0;
      for(int l = 0; l <= (length(kmer))-markovOrder; l++)
      {
         int j = l + markovOrder;
         string inf;
         assign(inf,infix(kmer, l, j));
         prob = prob * ((double)markovcounts[inf] / (double)tot);
      }
      markovmap[kmer] = prob;
   }

   return markovmap;
};



map<string, bool> makecomplete(ModifyStringOptions options);


// New create all kmer function
template <typename TAlphabet>
vector<String<TAlphabet>> makecomplete(unsigned int klen, TAlphabet const &)
{
   vector<String<TAlphabet>> result;

   // create our alphabet
   String<TAlphabet> bases;
   typedef typename ValueSize<TAlphabet>::Type TSize;
   TSize alphSize = ValueSize<TAlphabet>::VALUE;

   // the minus one prevents the 'N' from being added,
   // therefore making Dna5 into Dna
   for(TSize i = 0; i < alphSize - 1; ++i)
      appendValue(bases, TAlphabet(i));

   // store temp data, populate with alphabet
   vector<String<TAlphabet>> temp;

   for(auto b : bases)
      appendValue(temp, b);

   for(unsigned int i = 1; i < klen; i++)
   {
      for(auto s : temp)
      {
         for(auto b : bases)
         {
            String<TAlphabet> meh = s;
            meh += b;
            appendValue(result, meh);
         }
      }

      temp.clear();
      temp = result;
      result.clear();
   }

   return temp;

};


map<string, bool> makequick(ModifyStringOptions options, 
                            StringSet<String<AminoAcid>> referenceseqs);
int ipow(int base, int exp);
CharString namecut(CharString seq, int val);
int parseMask(ModifyStringOptions options, int &effectiveKlen);


template <typename TAlphabet>
int printPhylyp(ModifyStringOptions options,
                vector<pair<CharString, map<String<TAlphabet>, unsigned int>>> &pw_counts,
                vector<vector<double>> &array_threaded_internal)
{
   ofstream outfile_new;

   try
   {
      outfile_new.open(toCString(options.outputFileName), std::ios_base::out);
   }
   catch (const ifstream::failure& e)
   {
      cout << "Error: could not open output file ";
      cout << toCString(options.outputFileName) << endl;
      return 1;
   }

   if(options.outputFileName == NULL)
   {
      cout << pw_counts.size() << endl;
   }
   else
   {
      outfile_new << pw_counts.size() << endl;
   }

   for(int i = 0; i < pw_counts.size(); i++)
   {
      StringSet<CharString> split;
      strSplit(split, pw_counts[i].first);
      int cutsize = 10;
      CharString qName = split[0];
      qName = namecut(qName, cutsize);

      if(options.outputFileName == NULL)
      {
         cout << qName << "\t";
         for(int j = 0; j < pw_counts.size(); j++)
            cout << array_threaded_internal[i][j] << "\t";
         cout << endl;
      }
      else
      {
         outfile_new << qName << "\t";
         for(int j = 0; j < pw_counts.size(); j++)
            outfile_new << array_threaded_internal[i][j] << "\t";
         outfile_new << endl;
      }
   }

   close(outfile_new);
   return 0;
};

#endif
