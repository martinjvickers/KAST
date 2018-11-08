#include "distance.h"
//#include "utils.h"
#include "print.h"
#include <seqan/alignment_free.h>
#include <sys/sysinfo.h> 


/*
I need to check that if we are using skip-mers, then we need to check that 
these are sensible.
  * skip-mer mask must all be 0/1's
  * skip-mer mask should be the same number of characters as the kmer size
  * all skipmers shoud have the same number of 1's across masks
*/
/*
int parseMask(ModifyStringOptions options, int &effectiveKlen)
{
   bool first = true;

   //vector<CharString> mask;
   for(auto m : options.mask)
   {
      // all should be the same number of characters as the klen
      if(length(m) != options.klen)
      {
         cerr << "ERROR: Mask sizes should be the same size ";
         cerr << "as the K-Mer length." << endl;
         return 1;
      }

      int counter = 0;

      // checks to see that the mask is only made of 0/1's
      for(int i = 0; i < length(m); i++)
      {
         if(m[i] != '0' && m[i] != '1')
         {
            cerr << "ERROR: Masks should only contain 0's or 1's." << endl;
            return 1;
         }

         if(m[i] == '1')
            counter++;
      }

      if(first == true)
      {
         effectiveKlen = counter;
         first = false;
      }
      else
      {
         if(counter != effectiveKlen)
         {
            cerr << "ERROR: The number of 0's and 1's in each mask ";
            cerr << "should be the same e.g. 10001, 11000, 00011" << endl;
            return 1;
         }
      }
   }
   return 0;
};
*/

/*markov is the alfsc implementation of the markov calculation;
  A much simplier version than the other, but I believe this
  is simply called a Background Model but i could be wrong.

*/

// for all others
/*
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
*/

// for DNA sequences
/*
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

template <typename TAlphabet>
void markov_again(String<double> & markovCounts, String<unsigned> const & kmerCounts,
                  String<TAlphabet> const & sequence, unsigned const k, unsigned const markovOrder)
{
}

/*
The alternative markov implementation seen in d2tools

I think there is an issue here where when the markovOrder is 0, then we 
don't reverse complement the count.

0 1 == true
1 2 == false
2 3 == false
3 4 == false

*/
template <>
void markov_again<>(String<double> & markovCounts, String<unsigned> const & kmerCounts,
                    String<Dna5> const & sequence, unsigned const k, unsigned const markovOrder)
{
   // setup markovCounts
   Shape<Dna> myShape;
   resize(myShape, k);
   int kmerNumber = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(myShape));
   seqan::clear(markovCounts);
   seqan::resize(markovCounts, kmerNumber, 0);

   if(markovOrder == 0)
   {
      String<unsigned> markovbg;
      // this should be run on the original sequence 
      countKmersNew(markovbg, sequence, markovOrder+1);
      unsigned tot = 0;

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
   else
   {
      String<unsigned> markovbg, markovbg_1;
      countKmersNew(markovbg, sequence, markovOrder);
      countKmersNew(markovbg_1, sequence, markovOrder+1);
      unsigned tot = 0;
      for(unsigned i = 0; i < length(markovbg); i++)
         tot = tot + markovbg[i];

      
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

/*
   Estimate how much RAM is needed to calculate this
*/
template <typename TAlphabet>
int mem_check(ModifyStringOptions options, int numRecord, TAlphabet const & alphabetType)
{
   struct sysinfo myinfo; 
   unsigned long total_bytes; 
   sysinfo(&myinfo); 

   total_bytes = myinfo.mem_unit * myinfo.totalram; 

   const unsigned long GIGABYTE = 1024 * 1024 * 1024;

   // calculate kmer count array needed
   unsigned alphSize = ValueSize<TAlphabet>::VALUE;
   if(alphSize == 5)
      alphSize = 4; // remember we initially work with Dna5 but only storage Dna
   
   // to store the counts of a single fasta entry
   unsigned long long int counts_mem = pow(alphSize, options.klen) * sizeof(unsigned);

   if((counts_mem * numRecord) > total_bytes)
   {
      cerr << "ERROR: Your machine has " << total_bytes/1024/1024 << " MB RAM but it will require approximately " << (counts_mem*numRecord)/1024/1024 << " MB to calculate." << endl;
      return 1;
   }
   return 0;
}

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

   // mem checker
   if(mem_check(options, length(pwseqs), alphabetType) == 1)
      return 1;

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

         if(options.mask.size() > 0)
            countKmersNew(counts[i], seqf, options.klen, options.effectiveLength, options.mask);
         else
            countKmersNew(counts[i], seqf, options.klen);

         if(options.type == "d2s" || options.type == "d2star" ||
            options.type == "hao" || options.type == "dai")
         {
            resize(markovCounts, length(pwseqs));
            markov(markovCounts[i], counts[i], seqf, options.klen, options.markovOrder);
         }
      }
      else // when doing aa/raa we don't (and can't) do reverse complement
      {
         countKmersNew(counts[i], seq, options.klen);

         if(options.type == "d2s" || options.type == "d2star" ||
            options.type == "hao")
         {
            resize(markovCounts, length(pwseqs));
            markov(markovCounts[i], counts[i], seq, options.klen, options.markovOrder);
         }
         else if(options.type == "D2S" || options.type == "D2Star" || options.type == "dai")
         {
            resize(markovCounts, length(pwseqs));
            markov_again(markovCounts[i], counts[i], seq, options.klen, options.markovOrder);
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
         countKmersNew(counts_i, pwseqs[i], options.klen);
         countKmersNew(counts_j, pwseqs[j], options.klen);
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
