#ifndef UTILS_H
#define UTILS_H

#include "common.h"

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <unordered_map>

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv);

int countKmersNew(String<unsigned> & kmerCounts, Dna5String const & sequence, 
                   unsigned const k);

template <class T>
int safe_increment(T& value)
{
   if(value < numeric_limits<T>::max())
   {
      value++;
      return 0;
   }
   else
   {
      cerr << "Error: Integer overflow detected. Exiting" << endl;
      return 1;
   }
};


//void markov(String<double> & markovCounts, String<unsigned> const & kmerCounts,
//            String<Dna5> const & sequence, unsigned const k, unsigned const markovOrder);

//void markov(String<double> & markovCounts, String<unsigned> const & kmerCounts,
 //           String<AminoAcid> const & sequence, unsigned const k, unsigned const markovOrder);

//void markov(String<double> & markovCounts, String<unsigned> const & kmerCounts,
//            String<ReducedAminoAcid> const & sequence, unsigned const k, unsigned const markovOrder);


// this is the markov method used within ALFSC-python
/*
template <typename TAlphabet>
void markov(String<double> & markovCounts, String<unsigned> const & kmerCounts,
            String<TAlphabet> const & sequence, unsigned const k, unsigned const markovOrder)
{
   // setup markovCounts
   Shape<TAlphabet> myShape;
   resize(myShape, k);
   int kmerNumber = _intPow((unsigned)ValueSize<TAlphabet>::VALUE, weight(myShape));

   cout << typeid(TAlphabet).name() << endl;

   seqan::clear(markovCounts);
   seqan::resize(markovCounts, kmerNumber, 0);

   // Now create the background model
   String<unsigned> markovbg;
   countKmersNew(markovbg, sequence, markovOrder);
   unsigned tot = 0;

   // sum the occurances
   for(unsigned i = 0; i < length(markovbg); i++)
      tot = tot + markovbg[i];

   for(unsigned i = 0; i < length(markovCounts); i++)
   {
      String<TAlphabet> inf;
      unhash(inf, i, k);
      String<unsigned> occurances;
      countKmersNew(occurances, inf, markovOrder);
      double prob = 1.0;
      for(unsigned i = 0; i < length(occurances); i++)
      {
         prob = prob * pow(((double)markovbg[i]/(double)tot), occurances[i]);
      }
      markovCounts[i] = prob;
   }
}

template <>
void markov<>(String<double> & markovCounts, String<unsigned> const & kmerCounts,
            String<Dna5> const & sequence, unsigned const k, unsigned const markovOrder)
{
   // setup markovCounts
   Shape<Dna> myShape;
   resize(myShape, k);
   int kmerNumber = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(myShape));
   seqan::clear(markovCounts);
   seqan::resize(markovCounts, kmerNumber, 0);

   // Now create the background model
   String<unsigned> markovbg;
   countKmersNew(markovbg, sequence, markovOrder);
   unsigned tot = 0;

   // sum the occurances
   for(unsigned i = 0; i < length(markovbg); i++)
      tot = tot + markovbg[i];

   for(unsigned i = 0; i < length(markovCounts); i++)
   {
      String<Dna> inf;
      unhash(inf, i, k);
      String<unsigned> occurances;
      countKmersNew(occurances, inf, markovOrder);
      double prob = 1.0;
      for(unsigned i = 0; i < length(occurances); i++)
      {
         prob = prob * pow(((double)markovbg[i]/(double)tot), occurances[i]);
      }
      markovCounts[i] = prob;
   }
}
*/
#endif
