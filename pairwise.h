#include "distance.h"
#include "utils.h"
#include "print.h"
#include <seqan/alignment_free.h>

// for all others
template <typename TAlphabet>
void markov(String<double> & markovCounts, String<unsigned> const & kmerCounts,
            String<TAlphabet> const & sequence, unsigned const k, unsigned const markovOrder)
{
   // setup markovCounts
   Shape<TAlphabet> myShape;
   resize(myShape, k);
   int kmerNumber = _intPow((unsigned)ValueSize<TAlphabet>::VALUE, weight(myShape));

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

// for DNA sequences
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

template <typename TAlphabet>
int calcDistance(unsigned & rI, unsigned & cI,
              vector< vector<double> > & results,
              mutex & location,
              StringSet<String<unsigned> > & counts,
              StringSet<String<double> > & markovCounts,
              ModifyStringOptions options)
{
   unsigned row, column;

   while(1)
   {
      location.lock();
      row = rI;
      column = cI;

      if(cI < length(counts))
      {
         cI++;
      }
      else if(cI >= length(counts) && rI < length(counts))
      {
         cI = rI;
         rI++;
         location.unlock();
         continue;
      }
      else if(!(rI < length(counts)))
      {
         location.unlock();
         return 0;
      }
      location.unlock();

      if(row < length(counts) && column < length(counts))
      {
         double dist;

         if(options.type == "euclid")
            dist = euler(counts[row], counts[column]);
         else if(options.type == "d2")
            dist = d2(counts[row], counts[column]);
         else if(options.type == "manhattan")
            dist = manhattan(counts[row], counts[column]);
         else if(options.type == "chebyshev")
            dist = chebyshev(counts[row], counts[column]);
         else if(options.type == "bc")
            dist = bray_curtis_distance(counts[row], counts[column]);
         else if(options.type == "ngd")
            dist = normalised_google_distance(counts[row], counts[column]);
         else if(options.type == "d2s" || options.type == "D2S")
            dist = d2s(counts[row], counts[column], markovCounts[row], markovCounts[column]);
         else if(options.type == "hao")
            dist = hao(counts[row], counts[column], markovCounts[row], markovCounts[column]);
         else if(options.type == "d2star" || options.type == "D2Star")
            dist = d2star(counts[row], counts[column], markovCounts[row], markovCounts[column]);
         else if(options.type == "dai")
            dist = dai(counts[row], counts[column], markovCounts[row], markovCounts[column]);

         results[row][column] = dist;
         results[column][row] = dist;
      }
   }
};

template <typename TAlphabet>
int pairwise_matrix(ModifyStringOptions options, TAlphabet const & alphabetType)
{
   // Read in the records into memory
   SeqFileIn pairwiseFileIn;
   StringSet<CharString> pwids;
   StringSet<String<TAlphabet>> pwseqs;

   if(!open(pairwiseFileIn, (toCString(options.pairwiseFileName))))
   {
      cerr << "Error: could not open file ";
      cerr << toCString(options.pairwiseFileName) << endl;
      return 1;
   }

   try
   {
      readRecords(pwids, pwseqs, pairwiseFileIn);
   }
   catch(Exception const & e)
   {
      std::cout << "ERROR: The --sequence-type which was selected was \"";
      std::cout << options.sequenceType << "\" however when reading the \"";
      std::cout << options.pairwiseFileName << "\" file we get the following error;" << endl;
      std::cout << e.what() << std::endl;
      return 1;
   }

   // set up elements for thread
   vector<thread> vectorOfThreads;
   mutex location;
   unsigned rI = 0;
   unsigned cI = 0;
   //unsigned int n = std::thread::hardware_concurrency();
   unsigned int cores = options.num_threads;

   // store the distance calculation results
   vector< vector<double> > results;
   results.resize(length(pwseqs), vector<double>(length(pwseqs), 0.0));

   // store the kmer counts in RAM
   StringSet<String<unsigned> > counts;
   resize(counts, length(pwseqs));

   StringSet<String<double> > markovCounts;

   for(unsigned i = 0; i < length(pwseqs); i++)
   {
      String<TAlphabet> seq = pwseqs[i];

      // if it's DNA, then we usually do revC which means we need
      // need to specifically define String<Dna5> because reverseComplement
      // is only implemented for Dna/Dna5
      if(options.sequenceType == "dna")
      {
         String<Dna5> seqf = seq;

         if(options.noreverse == false)
         {
            String<Dna5> seqrc = seq;
            reverseComplement(seqrc);
            append(seqf, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"); // this should probably the same size as options.klen
            append(seqf, seqrc);
         }
         countKmers(counts[i], seqf, options.klen);

         if(options.type == "d2s" || options.type == "D2S" ||
            options.type == "d2star" || options.type == "D2Star" ||
            options.type == "hao" || options.type == "dai")
         {
            resize(markovCounts, length(pwseqs));
            markov(markovCounts[i], counts[i], seqf, options.klen, options.markovOrder);
         }
      }
      else // when doing aa/raa we don't (and can't) do reverse complement
      {
         countKmers(counts[i], seq, options.klen);
         if(options.type == "d2s" || options.type == "D2S" ||
            options.type == "d2star" || options.type == "D2Star" ||
            options.type == "hao" || options.type == "dai")
         {
            resize(markovCounts, length(pwseqs));
            markov(markovCounts[i], counts[i], seq, options.klen, options.markovOrder);
         }
      }
   }

   // do the distance calculations
   for(unsigned i = 0; i < cores; i++)
   {
      vectorOfThreads.push_back(thread(calcDistance<TAlphabet>, ref(rI), ref(cI), 
                                       ref(results), ref(location), ref(counts), 
                                       ref(markovCounts), options));
   }

   for(auto &thread : vectorOfThreads)
   {
      thread.join();
   }

   // Print results
   printPhylyp(options, pwids, results);

   return 0;
};

template <typename TAlphabet>
int pairwise_all_matrix(ModifyStringOptions options, TAlphabet const & alphabetType)
{
   // Read in the records into memory
   SeqFileIn pairwiseFileIn;
   StringSet<CharString> pwids;
   StringSet<String<TAlphabet>> pwseqs;

   if(!open(pairwiseFileIn, (toCString(options.pairwiseFileName))))
   {
      cerr << "Error: could not open file ";
      cerr << toCString(options.pairwiseFileName) << endl;
      return 1;
   }

   try
   {
      readRecords(pwids, pwseqs, pairwiseFileIn);
   }
   catch(Exception const & e)
   {
      std::cout << "ERROR: The --sequence-type which was selected was \"";
      std::cout << options.sequenceType << "\" however when reading the \"";
      std::cout << options.pairwiseFileName << "\" file we get the following error;" << endl;
      std::cout << e.what() << std::endl;
      return 1;
   }

   cout << "Q1\tQ2\tEuclid\td2\tManhattan\tBC\tNGD\tHao\tdai\tD2S\tD2Star\n";
   for(unsigned int i = 0; i < length(pwids); i++)
   {
      for(unsigned int j = 0; j < length(pwids); j++)
      {
         String<unsigned> counts_i, counts_j;
         countKmers(counts_i, pwseqs[i], options.klen);
         countKmers(counts_j, pwseqs[j], options.klen);
         String<double> markov_i, markov_j;
         markov(markov_i, counts_i, pwseqs[i], options.klen, options.markovOrder);
         markov(markov_j, counts_j, pwseqs[j], options.klen, options.markovOrder);
         
         cout << pwids[i] << "\t" << pwids[j] << "\t" << euler(counts_i, counts_j);
         cout << "\t" << d2(counts_i, counts_j) << "\t" << manhattan(counts_i, counts_j);
         cout << "\t" << bray_curtis_distance(counts_i, counts_j);
         cout << "\t" << normalised_google_distance(counts_i, counts_j);
         cout << "\t" << hao(counts_i, counts_j, markov_i, markov_j);
         cout << "\t" << dai(counts_i, counts_j, markov_i, markov_j);
         cout << "\t" << d2s(counts_i, counts_j, markov_i, markov_j);
         cout << "\t" << d2star(counts_i, counts_j, markov_i, markov_j);

         cout << endl;
      }
   }

   return 0;
};
