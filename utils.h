#ifndef UTILS_H
#define UTILS_H

#include "common.h"

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <unordered_map>

// Parse our commandline options
ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options,
                                             int argc, char const ** argv)
{
   ArgumentParser parser("kast");
   addOption(parser, ArgParseOption("k", "klen", "Kmer Length.",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "klen", "3");
   addOption(parser, ArgParseOption("d", "debug", "Debug Messages."));
   addOption(parser, ArgParseOption("q", "query-file",
                                    "Path to the file containing your query \
                                    sequence data.\n",
                                    ArgParseArgument::INPUT_FILE, "IN"));
   addOption(parser, ArgParseOption("r", "reference-file",
                                    "Path to the file containing your reference\
                                     sequence data.",
                                     ArgParseArgument::INPUT_FILE, "IN"));
   addOption(parser, ArgParseOption("p", "pairwise-file",
                                    "Path to the file containing your sequence \
                                    data which you will perform pairwise \
                                    comparison on.",
                                    ArgParseArgument::INPUT_FILE, "IN"));
   addOption(parser, ArgParseOption("m", "markov-order", "Markov Order",
             ArgParseArgument::INTEGER, "INT"));
   addOption(parser, ArgParseOption("o", "output-file", "Output file.",
             ArgParseArgument::OUTPUT_FILE, "OUT"));

   setDefaultValue(parser, "markov-order", "0");
   addOption(parser, ArgParseOption("n", "num-hits",
             "Number of top hits to return when running a Ref/Query search.\
              If you want all the result, enter 0.", ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "num-hits", "10");
   addOption(parser, ArgParseOption("t", "distance-type",
                                    "The method of calculating the distance \
                                    between two sequences. For descriptions of \
                                    distance please refer to the wiki. ",
                                    ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "distance-type",
                  "d2 euclid d2s d2star manhattan chebyshev dai bc ngd all");
   //            "d2 euclid d2s D2S d2s-opt d2star D2Star manhattan chebyshev hao dai bc ngd all new");
   setDefaultValue(parser, "distance-type", "d2");
   addOption(parser, ArgParseOption("sc", "score-cutoff", "Score Cutoff for search mode.",
             ArgParseArgument::DOUBLE, "DOUBLE"));
   if(isSet(parser, "score-cutoff") == false)
   {
      options.score_cutoff = std::numeric_limits<double>::quiet_NaN();
   }
   addOption(parser, ArgParseOption("s", "sequence-type",
             "Define the type of sequence data to work with.",
             ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "sequence-type", "dna aa raa");
   setDefaultValue(parser, "sequence-type", "dna");
   addOption(parser, ArgParseOption("f", "output-format",
             "For Reference/query based usage you can select your output type.",
             ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "output-format", "default tabular blastlike");
   setDefaultValue(parser, "output-format", "default");
   addOption(parser, ArgParseOption("nr", "no-reverse",
             "Do not use reverse compliment."));
   addOption(parser, ArgParseOption("gc", "calc-gc",
             "Calculate GC content of query/ref in search mode."));
   addOption(parser, ArgParseOption("nh", "no-header",
             "Do not print header when performing search mode."));
   addOption(parser, ArgParseOption("mask", "skip-mer",
             "Specify binary masks where a zero indicates \
             skipping that base and one keeps it. e.g. 01110.",
             ArgParseArgument::STRING, "TEXT", true));
   addOption(parser, ArgParseOption("c", "num-cores", "Number of Cores.",
             ArgParseArgument::INTEGER, "INT"));
/*   addOption(parser, ArgParseOption("l", "low-ram", 
             "Does not store the reference in RAM. \
              As long as you're not using a very \
              large kmer size, this option will \
              allow you to run kast with a large \
              reference, however it will take much \
              longer."));*/
   setDefaultValue(parser, "num-cores", "1");
   setShortDescription(parser, "Kmer Alignment-free Search Tool.");
   setVersion(parser, "0.0.29");
   setDate(parser, "March 2019");
   addUsageLine(parser, "-q query.fasta -r reference.fasta -o results.txt [\\fIOPTIONS\\fP] ");
   addUsageLine(parser, "-p mydata.fasta -o results.txt [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Perform Alignment-free k-tuple frequency \
                           comparisons from sequences. This can be in \
                           the form of two input files (e.g. a \
                           reference and a query) or a single file for \
                           pairwise comparisons to be made.");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   // Only extract options if the program will continue after parseCommandLine
   if (res != seqan::ArgumentParser::PARSE_OK)
      return res;

   //begin extracting options
   getOptionValue(options.klen, parser, "klen");
   getOptionValue(options.score_cutoff, parser, "score-cutoff");
   getOptionValue(options.nohits, parser, "num-hits");
   getOptionValue(options.markovOrder, parser, "markov-order");
   getOptionValue(options.type, parser, "distance-type");
   getOptionValue(options.sequenceType, parser, "sequence-type");
   options.noreverse = isSet(parser, "no-reverse");
   options.calcgc = isSet(parser, "calc-gc");
   options.noheader = isSet(parser, "no-header");
   options.debug = isSet(parser, "debug");
   //options.lowram = isSet(parser, "low-ram");
   getOptionValue(options.queryFileName, parser, "query-file");
   getOptionValue(options.referenceFileName, parser, "reference-file");
   getOptionValue(options.pairwiseFileName, parser, "pairwise-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   getOptionValue(options.num_threads, parser, "num-cores");
   getOptionValue(options.output_format, parser, "output-format");

   for(int i = 0; i < getOptionValueCount(parser, "skip-mer"); i++)
   {
      CharString tmpVal;
      getOptionValue(tmpVal, parser, "skip-mer", i);
      options.mask.push_back(tmpVal);
   }

/*
   Check to see if markov order is correct. 
*/
   if(options.type == "d2s" || options.type == "d2star" || 
      options.type == "dai" || options.type == "hao" || 
      options.type == "all")
   {
      if(options.markovOrder < 0)
      {
         cerr << "Markov order must be >= 0 " << endl;
         return ArgumentParser::PARSE_ERROR;
      }

      if(options.mask.size() > 0)
      {
         if(options.markovOrder >= options.effectiveLength-1)
         {
            cerr << "Markov order must be < effectiveLength-1 " << endl;
            return ArgumentParser::PARSE_ERROR;
         }
      }
      else
      {
         if(options.markovOrder >= options.klen-1)
         {
            cerr << "Markov order must be < klen-1 " << endl;
            return ArgumentParser::PARSE_ERROR;
         }
      }
   }

   if(isSet(parser, "pairwise-file"))
   {
      if(isSet(parser, "reference-file") == true ||
         isSet(parser, "query-file") == true)
      {
         cerr << "If you are performing a pairwise comparison, ";
         cerr << "you do not need to specify a query (-q) and a ";
         cerr << "reference (-r) file. If you are performing a ";
         cerr << "reference/query based search you do not need to ";
         cerr << "specify a pairwise-file (-p). See kast -h for ";
         cerr << "details." << endl;
         return ArgumentParser::PARSE_ERROR;
      }
   }

   if(isSet(parser, "reference-file") == true &&
      isSet(parser, "query-file") == false)
   {
      cerr << "You have specified a reference (-r) file but ";
      cerr << "not a query (-q) file. See kast -h for details." << endl;
      printHelp(parser);
      return ArgumentParser::PARSE_ERROR;
   }

   if(isSet(parser, "reference-file") == false &&
      isSet(parser, "query-file") == true)
   {
      cerr << "You have specified a query (-q) file but ";
      cerr << "not a reference (-r) file. See kast -h for details." << endl;
      printHelp(parser);
      return ArgumentParser::PARSE_ERROR;
   }
   if(isSet(parser, "reference-file") == false &&
      isSet(parser, "query-file") == false &&
      isSet(parser, "pairwise-file") == false)
   {
      cerr << "You have not specifed any input file. ";
      cerr << "See kast -h for details." << endl;
      printHelp(parser);
      return ArgumentParser::PARSE_ERROR;
   }
   return ArgumentParser::PARSE_OK;
}

/*
I need to check that if we are using skip-mers, then we need to check that 
these are sensible.
  * skip-mer mask must all be 0/1's
  * skip-mer mask should be the same number of characters as the kmer size
  * all skipmers shoud have the same number of 1's across masks
*/
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

template <typename TAlphabet>
double gc_ratio(String<TAlphabet> sequence)
{
   int gc = 0;
   int agct = 0;

   for(int i = 0; i < length(sequence); i++)
   {
      if(sequence[i] == 'G' || sequence[i] == 'C')
         gc++;
      if(sequence[i] != 'N')
         agct++;
   }
   return (double)gc/(double)agct;
}

/*
Perform regular counting
*/
template <typename TAlphabet>
int countKmersNew(String<unsigned> & kmerCounts, String<TAlphabet> const & sequence, unsigned const k)
{
   Shape<TAlphabet> myShape;
   resize(myShape, k);
   unsigned long long int kmerNumber = _intPow((unsigned)ValueSize<TAlphabet>::VALUE, weight(myShape));

   seqan::clear(kmerCounts);
   seqan::resize(kmerCounts, kmerNumber, 0);

   auto itSequence = begin(sequence);

   for (; itSequence <= (end(sequence) - k); ++itSequence)
   {
      long long unsigned int hashValue = seqan::hash(myShape, itSequence);
      safe_increment(kmerCounts[hashValue]);
   }
   return 0;
};

template <>
int countKmersNew(String<unsigned> & kmerCounts, String<Dna5> const & sequence, unsigned const k)
{
   Shape<Dna> myShape;
   resize(myShape, k);
   int kmerNumber = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(myShape));
   seqan::clear(kmerCounts);
   seqan::resize(kmerCounts, kmerNumber, 0);

   auto itSequence = begin(sequence);
   int counterN = 0;

   // Check for any N that destroys the first kmers
   unsigned j = k - 1;
   for (auto i = position(itSequence); i <= j; ++i)
   {
       if (_repeatMaskValue(sequence[i]))
       {
           counterN = i + 1;
       }
   }

   for (; itSequence <= (end(sequence) - k); ++itSequence)
   {
       // Check if there is a "N" at the end of the new kmer
       if (_repeatMaskValue(value(itSequence + (k - 1))))
           counterN = k;  // Do not consider any kmer covering this "N"
       // If there is no "N" overlapping with the current kmer, count it
       if (counterN <= 0)
       {
           unsigned hashValue = seqan::hash(myShape, itSequence);
           safe_increment(kmerCounts[hashValue]);
       }
       counterN--;
   }

   return 0;
};


/*
Perform counting with a mask array
For AminoAcid and ReducedAminoAcid
*/
template <typename TAlphabet>
int countKmersNew(String<unsigned> & kmerCounts, String<TAlphabet> const & sequence, 
                  unsigned const k, unsigned const effectiveK,
                  vector<CharString> const & mask)
{
   /*Build the shape to traverse the sequence (this is the full kmer size)*/
   Shape<TAlphabet> myShape;
   resize(myShape, k);

   auto itSequence = begin(sequence);
   int counterN = 0;

   // Now create the effective size shape
   Shape<TAlphabet> maskShape;
   resize(maskShape, effectiveK);

   // resize the kmerCounts vector for this entry which 
   int kmerNumber = _intPow((unsigned)ValueSize<TAlphabet>::VALUE, weight(maskShape));
   seqan::clear(kmerCounts);
   seqan::resize(kmerCounts, kmerNumber, 0);

   // go through the sequence using k=X at a time
   for (; itSequence <= (end(sequence) - k); ++itSequence)
   {
      unsigned hashValue = seqan::hash(myShape, itSequence);
      String<TAlphabet> orig;
      unhash(orig, hashValue, k);

      //cout << orig << endl;

      // here, I need to loop
      for(int i = 0; i < mask.size(); i++)
      {
         String<TAlphabet> dnaSeq;
         for(int m = 0; m < length(mask[i]); m++)
         {
            if(mask[i][m] == '1')
            {
               dnaSeq += orig[m];
            }
         }
         //cout << dnaSeq << endl;
         auto it = begin(dnaSeq);
         unsigned hashMask = seqan::hash(maskShape, it);
         safe_increment(kmerCounts[hashMask]);
      }
   }
   return 0;
};

// For Dna
template <>
int countKmersNew(String<unsigned> & kmerCounts, String<Dna5> const & sequence,
                  unsigned const k, unsigned const effectiveK,
                  vector<CharString> const & mask)
{
   /*Build the shape to traverse the sequence (this is the full kmer size)*/
   Shape<Dna> myShape;
   resize(myShape, k);

   auto itSequence = begin(sequence);
   int counterN = 0;

   // Now create the effective size shape
   Shape<Dna> maskShape;
   resize(maskShape, effectiveK);

   int kmerNumber = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(maskShape));
   seqan::clear(kmerCounts);
   seqan::resize(kmerCounts, kmerNumber, 0);

   // Check for any N that destroys the first kmers
   unsigned j = k - 1;
   for (auto i = position(itSequence); i <= j; ++i)
   {
       if (_repeatMaskValue(sequence[i]))
       {
           counterN = i + 1;
       }
   }

   for (; itSequence <= (end(sequence) - k); ++itSequence)
   {
       // Check if there is a "N" at the end of the new kmer
       if (_repeatMaskValue(value(itSequence + (k - 1))))
           counterN = k;  // Do not consider any kmer covering this "N"

       // If there is no "N" overlapping with the current kmer, count it
       if (counterN <= 0)
       {
           unsigned hashValue = seqan::hash(myShape, itSequence);
           DnaString orig;
           unhash(orig, hashValue, k);

           // here, I need to loop
           for(int i = 0; i < mask.size(); i++)
           {
              String<Dna> dnaSeq;
              for(int m = 0; m < length(mask[i]); m++)
              {
                 if(mask[i][m] == '1')
                 {
                    dnaSeq += orig[m];
                 }
              }

              auto it = begin(dnaSeq);
              unsigned hashMask = seqan::hash(maskShape, it);
              safe_increment(kmerCounts[hashMask]);
           }
       }
       counterN--;
   }
   return 0;
};

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
   countKmersNew(markovbg, sequence, markovOrder+1);
   unsigned tot = 0;

   // sum the occurances
   for(unsigned i = 0; i < length(markovbg); i++)
      tot = tot + markovbg[i];

   for(unsigned i = 0; i < length(markovCounts); i++)
   {
      String<TAlphabet> inf;
      unhash(inf, i, k);
      String<unsigned> occurances;
      countKmersNew(occurances, inf, markovOrder+1);
      double prob = 1.0;
      for(unsigned i = 0; i < length(occurances); i++)
      {
         prob = prob * pow(((double)markovbg[i]/(double)tot), occurances[i]);
      }
      markovCounts[i] = prob;
   }
};

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
   countKmersNew(markovbg, sequence, markovOrder+1);
   unsigned long long int tot = 0;

   // sum the occurances
   for(unsigned i = 0; i < length(markovbg); i++)
      tot = tot + markovbg[i];

   //cout << "Markov " << endl;

   for(unsigned i = 0; i < length(markovCounts); i++)
   {
      String<Dna> inf;
      unhash(inf, i, k);
      String<Dna5> meh = inf; // this conversion is important, as countKmers requires Dna5
      String<unsigned> occurances;
      countKmersNew(occurances, meh, markovOrder+1);
      double prob = 1.0;
      for(unsigned i = 0; i < length(occurances); i++)
      {
         prob = prob * pow(((double)markovbg[i]/(double)tot), occurances[i]);
      }
      markovCounts[i] = prob;
      //cout << inf << "\t" << prob << "\t" << endl;
   }
};

#endif
