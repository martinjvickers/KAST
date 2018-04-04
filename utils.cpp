#include "utils.h"

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
   setValidValues(parser, "query-file", 
                  toCString(concat(getFileExtensions(SeqFileIn()), ' ')));
   addOption(parser, ArgParseOption("r", "reference-file", 
                                    "Path to the file containing your reference\
                                     sequence data.", 
                                     ArgParseArgument::INPUT_FILE, "IN"));
   setValidValues(parser, "reference-file", 
                  toCString(concat(getFileExtensions(SeqFileIn()), ' ')));
   addOption(parser, ArgParseOption("p", "pairwise-file", 
                                    "Path to the file containing your sequence \
                                    data which you will perform pairwise \
                                    comparison on.", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setValidValues(parser, "pairwise-file", 
                  toCString(concat(getFileExtensions(SeqFileIn()), ' ')));
   addOption(parser, ArgParseOption("m", "markov-order", "Markov Order", 
             ArgParseArgument::INTEGER, "INT"));
   addOption(parser, ArgParseOption("o", "output-file", "Output file.", 
             ArgParseArgument::OUTPUT_FILE, "OUT"));
   
   //setRequired(parser, "output-file");
   setDefaultValue(parser, "markov-order", "1");
   addOption(parser, ArgParseOption("n", "num-hits", 
             "Number of top hits to return", ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "num-hits", "10");
   addOption(parser, ArgParseOption("t", "distance-type", 
                                    "The method of calculating the distance \
                                    between two sequences. For descriptions of \
                                    distance please refer to the wiki ",
                                    ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "distance-type", 
                  "d2 kmer d2s d2s-opt d2star manhattan \
                   chebyshev hao dai bc ngd");
   setDefaultValue(parser, "distance-type", "d2");
   addOption(parser, ArgParseOption("s", "sequence-type", 
             "Define the type of sequence data to work with.", 
             ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "sequence-type", "nucl prot");
   setDefaultValue(parser, "sequence-type", "nucl");
   addOption(parser, ArgParseOption("f", "output-format",
             "For Reference/query based usage you can select your output type.",             ArgParseArgument::STRING, "STR"));
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
   addOption(parser, ArgParseOption("l", "low-ram", 
             "Does not store the reference in RAM. \
              As long as you're not using a very \
              large kmer size, this option will \
              allow you to run kast with a large \
              reference, however it will take much \
              longer."));
   setDefaultValue(parser, "num-cores", "1");
   setShortDescription(parser, "Kmer Alignment-free Search Tool.");
   setVersion(parser, "0.0.17");
   setDate(parser, "April 2018");
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
   options.lowram = isSet(parser, "low-ram");
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
}

AminoAcid getRevCompl(AminoAcid const & nucleotide)
{
   if(nucleotide == 'A')
      return 'T';
   if(nucleotide == 'T')
      return 'A';
   if(nucleotide == 'C')
      return 'G';
   if(nucleotide == 'G')
      return 'C';
   return 'N';
}

String<AminoAcid> doRevCompl(String<AminoAcid> seq)
{
   String<AminoAcid> allSeq;
   append(allSeq,seq);
   allSeq += "NNN";
   String<AminoAcid> revComplGenome;
   resize(revComplGenome, length(seq));
   
   for (unsigned i = 0; i < length(seq); ++i)
      revComplGenome[length(seq) - 1 - i] = getRevCompl(seq[i]);

   allSeq += revComplGenome;
   return allSeq;
}

/*
   Return a map of all the kmers in a sequence
*/
map<string, unsigned int> count(String<AminoAcid> sequence, int klen)
{
   int total = 0;
   map<string, unsigned int> map;
   for(int i = 0; i <= length(sequence)-klen; i++)
   {
      string kmer;
      assign(kmer,infix(sequence, i, i+klen));
      size_t found = kmer.find("N");

      if(found > kmer.size())
      {
         long long int count = map[kmer];
         map[kmer] = count + 1;
         total++;
      }
   }
   return map;
}

map<string, unsigned int> count(String<AminoAcid> sequence, int klen, 
                                vector<CharString> mask)
{
   int total = 0;
   map<string, unsigned int> map;
   for(int i = 0; i <= length(sequence)-klen; i++)
   {
      string kmer;
      assign(kmer,infix(sequence, i, i+klen));

      //so kmer has the wider bit we care about
      //now go through the mask and append kmers
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
         //total++;
      }
   }
   return map;
}

double gc_ratio(String<AminoAcid> sequence)
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

map<string, double> markov(int klen, String<AminoAcid> sequence, 
                           int markovOrder, map<string, bool> kmer_count_map)
{
   map<string, double> markovmap;
   double sum_prob = 0.0;

   map<string, unsigned int> markovcounts = count(sequence, markovOrder);
   double tot = 0;

   for(pair<string, unsigned int> p: markovcounts)
   {
      tot = tot + p.second;
   }

   for(pair<string, bool> p: kmer_count_map)
   {
      double prob = 1.0;
      Dna5String kmer = p.first;

      for(int l = 0; l <= (length(kmer))-markovOrder; l++)
      {
         int j = l + markovOrder;
         string inf;
         assign(inf,infix(kmer, l, j));
         prob = prob * ((double)markovcounts[inf] / (double)tot);
      }
      markovmap[p.first] = prob;
   }
   return markovmap;
}

map<string, bool> makequick(ModifyStringOptions options, 
                            StringSet<String<AminoAcid>> referenceseqs)
{
   map<string, bool> quickmers;

   for(int i = 0; i < length(referenceseqs); i++)
   {
      for(int j = 0; j <= length(referenceseqs[i])-options.effectiveLength; j++)
      {
         // end quickly if done	
         //int max = ipow(4, options.klen);
         int max = ipow(4, options.effectiveLength);
         if( quickmers.size() == max )
            return quickmers;

         string kmer;
         assign(kmer, infix(referenceseqs[i], j, j + options.effectiveLength));
         size_t found = kmer.find("N");

         if(found > kmer.size())
            quickmers[kmer] = true;
      }
   }
   return quickmers;
}

map<string, bool> makecomplete(ModifyStringOptions options)
{
   map<string, bool> finkmers;

   String<Dna5String> bases;
   appendValue(bases, "A");
   appendValue(bases, "G");
   appendValue(bases, "C");
   appendValue(bases, "T");

   String<Dna5String> kmers;
   kmers = bases;

   for(int j = 0; j < options.effectiveLength-1; j++)
   {
      String<Dna5String> temp;

      for(int k = 0; k < length(kmers); k++)
      {
         for(int m = 0; m < length(bases); m++)
         {
            String<Dna> kmer = bases[m];
            kmer += kmers[k];
            appendValue(temp,kmer);
         }
      }
      kmers = temp;
      clear(temp);
   }

   for(int i = 0; i < length(kmers); i++)
   {
      string kmer2;
      assign(kmer2,infix(kmers[i], 0, length(kmers[i])));
      finkmers[kmer2] = true;
   }
   return finkmers;
}

int ipow(int base, int exp)
{
   int result = 1;
   while(exp)
   {
      if(exp & 1)
      {
         result *= base;
      }
      exp >>= 1;
      base *= base;
   }
   return result;
}

CharString namecut(CharString seq, int val)
{
   if(length(seq) == val)
      return seq;

   if(length(seq) < val)
   {
      CharString charSeq;
      resize(charSeq, val, Exact());
      assign(charSeq, seq);

      for(int i = 0; i < val-length(seq); i++)
         insert(charSeq,length(seq)+i," ");

      return charSeq;
   }

   if(length(seq) > val)
   {
      CharString charSeq;
      resize(charSeq, val, Exact());
      assign(charSeq, seq, Limit());

      return charSeq;
   }
}

int printPhylyp(ModifyStringOptions options, 
                vector<pair<CharString, map<string, unsigned int>>> pw_counts,
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
}
