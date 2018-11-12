//#include "utils.h"

// Parse our commandline options
/*
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

   setDefaultValue(parser, "markov-order", "1");
   addOption(parser, ArgParseOption("n", "num-hits", 
             "Number of top hits to return", ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "num-hits", "10");
   addOption(parser, ArgParseOption("t", "distance-type", 
                                    "The method of calculating the distance \
                                    between two sequences. For descriptions of \
                                    distance please refer to the wiki. ",
                                    ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "distance-type", 
               "d2 euclid d2s D2S d2s-opt d2star D2Star manhattan chebyshev hao dai bc ngd all new");
   setDefaultValue(parser, "distance-type", "d2");
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
   addOption(parser, ArgParseOption("mask", "skip-mer", 
             "Specify binary masks where a zero indicates \
             skipping that base and one keeps it. e.g. 01110.", 
             ArgParseArgument::STRING, "TEXT", true));
   addOption(parser, ArgParseOption("c", "num-cores", "Number of Cores.", 
             ArgParseArgument::INTEGER, "INT"));
*/
/*   addOption(parser, ArgParseOption("l", "low-ram", 
             "Does not store the reference in RAM. \
              As long as you're not using a very \
              large kmer size, this option will \
              allow you to run kast with a large \
              reference, however it will take much \
              longer."));*/
/*
   setDefaultValue(parser, "num-cores", "1");
   setShortDescription(parser, "Kmer Alignment-free Search Tool.");
   setVersion(parser, "0.0.21");
   setDate(parser, "November 2018");
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
   getOptionValue(options.nohits, parser, "num-hits");
   getOptionValue(options.markovOrder, parser, "markov-order");
   getOptionValue(options.type, parser, "distance-type");
   getOptionValue(options.sequenceType, parser, "sequence-type");
   options.noreverse = isSet(parser, "no-reverse");
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

   if((options.markovOrder < 1) || (options.markovOrder > 3))
   {
      cerr << "Markov Order --markov-order should be ";
      cerr << "an integer 1, 2 or 3." << endl;
      return ArgumentParser::PARSE_ERROR;
   }

   // this is stopping me using k=1 for d2/euclid etc
   if(options.markovOrder > options.klen-1)
   {
      cerr << "Markov Order needs to be smaller than klen-1 and ";
      cerr << "an integer 1, 2 or 3. If you're using a small klen size, ";
      cerr << "decrease the size of --markov-order." << endl;
      return ArgumentParser::PARSE_ERROR;
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
*/

/*
int countKmersNew(String<unsigned> & kmerCounts, Dna5String const & sequence, unsigned const k)
{
   cout << "Meh " << endl;
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
           //++kmerCounts[hashValue];
           safe_increment(kmerCounts[hashValue]);

           DnaString orig;
           unhash(orig, hashValue, k);
           cout << "Orig: " << orig << endl;
       }
       counterN--;
   }

   return 0;
}

int countKmersNew(String<unsigned> & kmerCounts, Dna5String const & sequence, 
                  unsigned const k, unsigned const effectiveK, 
                  vector<CharString> const & mask)
{
   Shape<Dna> myShape;
   resize(myShape, k);
   int kmerNumber = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(myShape));
   seqan::clear(kmerCounts);
   seqan::resize(kmerCounts, kmerNumber, 0);

   auto itSequence = begin(sequence);
   int counterN = 0;

   Shape<Dna> maskShape;
   resize(maskShape, effectiveK);

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
                 //cout << mask[i] << "\t" << length(mask[i]) << endl;
                 if(mask[i][m] == '1')
                 {
                    dnaSeq += orig[m];
                    //cout << "Adding " << orig[m] << endl;
                 }
              }

              auto it = begin(dnaSeq);
              unsigned hashMask = seqan::hash(maskShape, it);
              //cout << "Orig: " << orig << "\tKmer: " << dnaSeq << "\tMask: " << mask[i] << endl;
              safe_increment(kmerCounts[hashMask]);
           }
       }
       counterN--;
   }
   return 0;
}
*/
