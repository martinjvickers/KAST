#ifndef UTILS_H
#define UTILS_H

#include "common.h"

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <unordered_map>
//#include <seqan/reduced_aminoacid.h>

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

//typedef SimpleType<unsigned char, ReducedAminoAcid_<Murphy10> > ReducedAminoAcidMurphy10;
map<String<ReducedAminoAcidMurphy10>, unsigned int> count_test(String<ReducedAminoAcidMurphy10> seq, int klen, bool noreverse);

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
map<String<TAlphabet>, double> markov_test(unsigned int klen, String<TAlphabet> &sequence,
                                           unsigned int markovOrder, 
                                           vector<String<TAlphabet>> &kmer_count_map,
                                           bool noreverse)
{
   map<String<TAlphabet>, double> markovmap;
   double sum_prob = 0.0;

   if(markovOrder == 0)
   {
      map<String<TAlphabet>, unsigned int> markovcounts = count_test(sequence, markovOrder+1, (markovOrder+1)<=1);
      double tot = 0;

      // sum of markov counts
      for(auto p: markovcounts)
      {
         tot = tot + p.second;
      }

      for(auto kmer : kmer_count_map)
      {
         double prob = 1.0;
         for(int l = 0; l <= (length(kmer))-markovOrder-1; l++)
         {
            int j = l + (markovOrder+1);
            string inf;
            assign(inf,infix(kmer, l, j));
            prob = prob * ((double)markovcounts[inf] / (double)tot);
         }
         markovmap[kmer] = prob;
      }
   }
   else
   {
      map<String<TAlphabet>, unsigned int> markovcounts_more = count_test(sequence, markovOrder+1, (markovOrder+1)<=1);
      map<String<TAlphabet>, unsigned int> markovcounts = count_test(sequence, markovOrder, markovOrder<=1);
      double tot = 0;
      for(auto p: markovcounts)
         tot = tot + p.second;

      for(auto kmer : kmer_count_map)
      {
         string w_start;
         assign(w_start, infix(kmer, 0, markovOrder));
         double p_temp = markovcounts[w_start] / tot;
         for(int l = 0; l <= (length(kmer))-markovOrder-1; l++)
         {
            int j = l + (markovOrder);
            string inf;
            assign(inf,infix(kmer, l, j));

            // now I need to add A,G,C,T to the end of inf
            string w1;
            assign(w1,infix(kmer, l, j+1));

            int nk1 = markovcounts_more[w1];
            int nk2 = markovcounts_more[inf+'A']+markovcounts_more[inf+'G']+markovcounts_more[inf+'C']+markovcounts_more[inf+'T'];
            if(nk2 != 0)
               p_temp = p_temp * ((double)markovcounts_more[w1] / (double)nk2);
            else
               p_temp = 0;
         }
         markovmap[kmer] = p_temp;
      }
   }

   return markovmap;
};


/*
   This is the markov calculation used by ALFSC-python

*/
template <typename TAlphabet>
map<String<TAlphabet>, double> markov_old(unsigned int klen, String<TAlphabet> &sequence,
                                           unsigned int markovOrder,
                                           vector<String<TAlphabet>> &kmer_count_map,
                                           bool noreverse)
{
   double sum_prob = 0.0;
   map<String<TAlphabet>, double> markovmap;
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
         prob = (double)prob * (double)((double)markovcounts[inf] / (double)tot);
      }
      markovmap[kmer] = prob;
   }
   return markovmap;
};

// new template markov function
template <typename TAlphabet>
void markov_test(unsigned int klen, String<TAlphabet> &sequence,
                 unsigned int markovOrder,
                 vector<String<TAlphabet>> &kmer_count_map,
                 bool noreverse, map<String<TAlphabet>, double> &markovmap)
{
//   map<String<TAlphabet>, double> markovmap;
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
         prob = (double)prob * (double)((double)markovcounts[inf] / (double)tot);
      }
      markovmap[kmer] = prob;
   }

   //return markovmap;
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
