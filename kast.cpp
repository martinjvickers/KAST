/*
KAST - Kmer Alignment-free Search Tool
Version 0.0.17
Written by Dr. Martin Vickers (martin.vickers@jic.ac.uk)

MIT License

Copyright (c) 2018 Martin James Vickers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>
#include <seqan/store.h>
#include <string>
#include <thread>
#include <map>
#include <vector>
#include "common.h"
#include "distance.h"
#include "utils.h"

using namespace seqan;
using namespace std;

map<string, bool> kmer_count_map;
SeqFileIn queryFileIn;
mutex m, n, r;
vector<map<string, unsigned int>> ref_counts_vec;
vector<map<string, double>> ref_markov_vec;
StringSet<CharString> referenceids;
StringSet<String<AminoAcid>> referenceseqs;
ofstream outfile; //output file
int current_row = 0;
vector< vector<double> > array_threaded;

/*
   TODO: mjv08

   The way I think I'm going to tackle this, for now (brain dump)...

   I don't want to have to care at the readRecord level what type of 
   alphabet we have. But, we know the type of distance measure so in 
   theory we should be able to say, "you're using the wrong measure
   on your dataset", e.g ngd on a protein input. However, this may be
   something someone wishes to do, so it might not be a good idea.

   Also, defaults are supposed to be different, e.g. we don't do
   reverse compl on AminoAcids.

   The dream would be the following;

   Read the sequence into a TPattern or THost or TSequence (is this correct?!, 
   can I make anything up here?) I define TSeqquence like this;

   typedef Dna5 TAlphabet;
   typedef String<TAlphabet> TSequence;
   TSequence seq;

   // this allows me to check what the sequence alphabet is
   if(IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)

   typedef typename Size<TAlphabet>::Type TSize;
   TSize alphSize = ValueSize<TAlphabet>::VALUE;

   typedef typename Size<Dna5>::Type TSizeDna5;
   TSizeDna5 dnaalphSize = ValueSize<Dna5>::VALUE;

   cout << "HiYa"<< endl;
   cout << "Current Alphabet\t" << alphSize << "\tDna5:\t" << dnaalphSize << "\tAA:\t" << ValueSize<AminoAcid>::VALUE << endl;


*/
void distance_thread(map<string, unsigned int> q1, map<string, unsigned int> q2, int i1, int i2, vector< vector<double> > & array_threaded_internal, ModifyStringOptions options)
{
   double dist;

   if(options.type == "kmer")
      dist = euler(options, q1, q2);
   else if(options.type == "d2")
      dist = d2(options, q1, q2);
   else if(options.type == "manhattan")
      dist = manhattan(options, q1, q2);
   else if(options.type == "chebyshev")
      dist = chebyshev(options, q1, q2);
   else if(options.type == "bc")
      dist = bray_curtis_distance(options, q1, q2);
   else if(options.type == "ngd")
      dist = normalised_google_distance(options, q1, q2);

   array_threaded_internal[i1][i2] = dist;
   array_threaded_internal[i2][i1] = dist;
}

int pairwise_matrix_test(ModifyStringOptions options)
{
   /* Read in the input file */
   SeqFileIn pairwiseFileIn;
   CharString pwid;
   String<AminoAcid> pwseq;
   vector< vector<double> > array_threaded_internal;

   if(!open(pairwiseFileIn, (toCString(options.pairwiseFileName))))
   {
      cerr << "Error: could not open file ";
      cerr << toCString(options.pairwiseFileName) << endl;
      return 1;
   }

   /* Calculate counts  */
   vector<pair<CharString, map<string, unsigned int>>> pw_counts;

   while(!atEnd(pairwiseFileIn))
   {
      readRecord(pwid, pwseq, pairwiseFileIn);
      String<AminoAcid> refseq;
      if(options.noreverse == false)
         refseq = doRevCompl(pwseq);
      pw_counts.push_back(make_pair(pwid, count(refseq, options.klen)));
   }

   close(pairwiseFileIn);

   /* Store the results */
   array_threaded_internal.resize(pw_counts.size(), 
                                  vector<double>(pw_counts.size(), 0.0));

   // create kmer_count_map if doing a markov model
   if(options.type == "d2s" || options.type == "hao" ||
      options.type == "d2star" || options.type == "dai")
   {
      kmer_count_map = makecomplete(options);
   }

   vector<thread> vectorOfThreads;

   /* Calculate the distances */
   for(unsigned rI = 0; rI < pw_counts.size(); ++rI)
      for(unsigned cI = rI; cI < pw_counts.size(); ++cI)
         vectorOfThreads.push_back(thread(distance_thread, pw_counts[rI].second, pw_counts[cI].second, rI, cI, std::ref(array_threaded_internal), options));

   for(auto &thread : vectorOfThreads)
      thread.join();

   /* Print the results */
   printPhylyp(options, pw_counts, array_threaded_internal);

   return 0;

}

//main threaded loop
int mainloop(ModifyStringOptions options)
{
   while(1)
   {
      String<AminoAcid> queryseq;
      CharString queryid;

      m.lock();
      if(!atEnd(queryFileIn))
      {
         readRecord(queryid, queryseq, queryFileIn);
      }
      else
      {
         m.unlock();
         return 0;
      }
      m.unlock();

      map<double, int> results;

      String<AminoAcid> seq;
      if(options.noreverse == false)
         seq = doRevCompl(queryseq);
      else
         seq = queryseq;

      map<string, unsigned int> querycounts;
      if(options.mask.size() > 0)
         querycounts = count(seq, options.klen, options.mask);
      else
         querycounts = count(seq, options.klen);

      map<string, double> querymarkov;
      if(options.type == "d2s" || options.type == "hao" || 
         options.type == "d2star" || options.type == "dai" ||
         options.type == "d2s-opt")
      {
         querymarkov = markov(options.effectiveLength, seq, 
                              options.markovOrder, kmer_count_map);
      }

      if(options.lowram != true)
      {
         for(int j = 0; j < length(referenceids); j++)
         {
            double dist;

            if(options.type == "kmer")
               dist = euler(options, ref_counts_vec[j], querycounts);
            else if(options.type == "d2")
               dist = d2(options, ref_counts_vec[j], querycounts);
            else if(options.type == "d2s" || options.type == "d2s-opt")
               dist = d2s(options, kmer_count_map, ref_counts_vec[j], 
                          ref_markov_vec[j], querycounts, querymarkov);
            else if(options.type == "d2star")
               dist = d2star(options, kmer_count_map, ref_counts_vec[j], 
                             ref_markov_vec[j], querycounts, querymarkov);
            else if(options.type == "dai")
               dist = dai(options, kmer_count_map, ref_counts_vec[j], 
                          ref_markov_vec[j], querycounts, querymarkov);
            else if(options.type == "hao")
               dist = hao(options, kmer_count_map, ref_counts_vec[j], 
                          ref_markov_vec[j], querycounts, querymarkov);
            else if(options.type == "manhattan")
               dist = manhattan(options, ref_counts_vec[j], querycounts);
            else if(options.type == "chebyshev")
               dist = chebyshev(options, ref_counts_vec[j], querycounts);
            else if(options.type == "bc")
               dist = bray_curtis_distance(options, ref_counts_vec[j], 
                                           querycounts);
            else if(options.type == "ngd")
               dist = normalised_google_distance(options, ref_counts_vec[j], 
                                                 querycounts);

            results.insert(pair<double, int> (dist, j));
            if(results.size() > options.nohits)
            {
               map<double, int>::iterator it = results.end();
               results.erase(--it);
            }
         }
      }
      else
      {
         CharString refid;
         String<AminoAcid> refseq;
         SeqFileIn refFileIn;
         int counter = 0;
         if(!open(refFileIn, (toCString(options.referenceFileName))))
         {
            cerr << "Error: could not open file ";
            cerr << toCString(options.referenceFileName) << endl;
            return 1;
         }

         while(!atEnd(refFileIn))
         {
            readRecord(refid, refseq, refFileIn);

            String<AminoAcid> rseq;
            if(options.noreverse == false)
               rseq = doRevCompl(refseq);
            else
               rseq = refseq;

            map<string, unsigned int> refcounts;

            if(options.mask.size() > 0)
               count(rseq, options.klen, options.mask);
            else
               count(rseq, options.klen);

            map<string, double> refmarkov;
            if(options.type == "d2s" || options.type == "hao" || 
               options.type == "d2star" || options.type == "dai" ||
               options.type == "d2s-opt")
            {
               refmarkov = markov(options.effectiveLength, rseq, 
                                  options.markovOrder, kmer_count_map);
            }

            double dist;
            if(options.type == "kmer")
               dist = euler(options, refcounts, querycounts);
            else if(options.type == "d2")
               dist = d2(options, refcounts, querycounts);
            else if(options.type == "d2s" || options.type == "d2s-opt")
               dist = d2s(options, kmer_count_map, refcounts, refmarkov, 
                          querycounts, querymarkov);
            else if(options.type == "d2star")
               dist = d2star(options, kmer_count_map, refcounts, refmarkov, 
                             querycounts, querymarkov);
            else if(options.type == "dai")
               dist = dai(options, kmer_count_map, refcounts, refmarkov, 
                          querycounts, querymarkov);
            else if(options.type == "hao")
               dist = hao(options, kmer_count_map, refcounts, refmarkov, 
                          querycounts, querymarkov);
            else if(options.type == "manhattan")
               dist = manhattan(options, refcounts, querycounts);
            else if(options.type == "chebyshev")
               dist = chebyshev(options, refcounts, querycounts);
            else if(options.type == "bc")
               dist = bray_curtis_distance(options, refcounts, querycounts);
            else if(options.type == "ngd")
               dist = normalised_google_distance(options, refcounts, 
                                                 querycounts);

            results.insert(pair<double, int> (dist, counter));
            if(results.size() > options.nohits)
            {
               map<double, int>::iterator it = results.end();
               results.erase(--it);
            }
            counter++;
         }
      }

      n.lock();
      printResult(options, queryid, outfile, queryseq, results,
                  referenceids, referenceseqs);
      n.unlock();
   }
   return 0;
}

/*
int pwthread(ModifyStringOptions options, StringSet<CharString> pairwiseid,
             StringSet<String<AminoAcid>> pairwiseseq)
{
   //not sure about this loop condition
   while(current_row < length(pairwiseid))
   {
      //get row locally and increment it
      int i;
      r.lock();
      i = current_row;
      if(i >= length(pairwiseid))
         return 0;
      current_row++;
      r.unlock();

      String<AminoAcid> seq;
      if(options.noreverse == false)
         seq = doRevCompl(pairwiseseq[i]);

      map<string, unsigned int> querycounts;
      if(options.mask.size() > 0)
         querycounts = count(seq, options.klen, options.mask);
      else
         querycounts = count(seq, options.klen);

      map<string, double> querymarkov;
      if(options.type == "d2s" || options.type == "hao" || 
         options.type == "d2star" || options.type == "dai" ||
         options.type == "d2s-opt")
      {
         querymarkov = markov(options.effectiveLength, seq, 
                              options.markovOrder, kmer_count_map);
      }

      for(int j = 0; j < i; j++)
      {
         double dist;
         String<AminoAcid> refseq;
         if(options.noreverse == false)
            refseq = doRevCompl(pairwiseseq[j]);

         map<string, unsigned int> refcounts;
         if(options.mask.size() > 0)
            refcounts = count(refseq, options.klen, options.mask);
         else
            refcounts = count(refseq, options.klen);

         map<string, double> refmarkov;
         if(options.type == "d2s" || options.type == "hao" || 
            options.type == "d2star" || options.type == "dai" ||
            options.type == "d2s-opt")
         {
            refmarkov = markov(options.effectiveLength, refseq, 
                               options.markovOrder, kmer_count_map);
         }

         if(options.type == "kmer")
            dist = euler(options, refcounts, querycounts);
         else if(options.type == "d2")
            dist = d2(options, refcounts, querycounts);
         else if(options.type == "d2s" || options.type == "d2s-opt")
            dist = d2s(options, kmer_count_map, refcounts, refmarkov, 
                       querycounts, querymarkov);
         else if(options.type == "d2star")
            dist = d2star(options, kmer_count_map, refcounts, refmarkov, 
                          querycounts, querymarkov);
         else if(options.type == "dai")
            dist = dai(options, kmer_count_map, refcounts, refmarkov, 
                       querycounts, querymarkov);
         else if(options.type == "hao")
            dist = hao(options, kmer_count_map, refcounts, refmarkov, 
                       querycounts, querymarkov);
         else if(options.type == "manhattan")
            dist = manhattan(options, refcounts, querycounts);
         else if(options.type == "chebyshev")
            dist = chebyshev(options, refcounts, querycounts);
         else if(options.type == "bc")
            dist = bray_curtis_distance(options, refcounts, querycounts);
         else if(options.type == "ngd")
            dist = normalised_google_distance(options, refcounts, querycounts);
		
         array_threaded[i][j] = dist;
         array_threaded[j][i] = dist;
      }
   }
   return 0;
}

int threaded_pw(ModifyStringOptions options)
{
   //so I'd need a pairwise set of records
   SeqFileIn pairwiseFileIn;
   StringSet<String<AminoAcid>> pairwiseseq;
   StringSet<CharString> pairwiseid;

   if(!open(pairwiseFileIn, (toCString(options.pairwiseFileName))))
   {
      cerr << "Error: could not open file ";
      cerr << toCString(options.pairwiseFileName) << endl;
      return 1;
   }

   readRecords(pairwiseid, pairwiseseq, pairwiseFileIn);

   //sort out our global matrix
   const int size = length(pairwiseid);
   int last = 1;
   array_threaded.resize(size);
   for(int i = 0; i < size; i++)
      array_threaded[i].resize(size);

   // create kmer_count_map
   if(options.type == "d2s" || options.type == "hao" ||
      options.type == "d2star" || options.type == "dai")
   {
      kmer_count_map = makecomplete(options);
   }
   else if(options.type == "d2s-opt")
   {
      kmer_count_map = makequick(options, pairwiseseq);
   }

   //run our threads, this is where we do the work
   int arraySize = options.num_threads;
   thread * workers = new thread[arraySize];
   for(int w = 0; w < options.num_threads; w++)
   {
      workers[w] = thread(pwthread, options, pairwiseid, pairwiseseq);
   }

   //do not exit until all the threads have finished
   for(int w = 0; w < options.num_threads; w++)
   {
      workers[w].join();
   }

   //print out in phylyp mode
   if(options.phylyp == true)
   {
      if(options.outputFileName == NULL)
      {
         cout << length(pairwiseid) << endl;
      }
      else
      {
         outfile << length(pairwiseid) << endl;
      }

      for(int i = 0; i < length(pairwiseid); i++)
      {
         StringSet<CharString> split;
         strSplit(split, pairwiseid[i]);
         int cutsize = 10;
         CharString qName = split[0];
         qName = namecut(qName, cutsize);

         if(options.outputFileName == NULL)
         {
            cout << qName << "\t";
            for(int j = 0; j < length(pairwiseid); j++)
               cout << array_threaded[i][j] << "\t";
            cout << endl;
         }
         else
         {
               outfile << qName << "\t";
               for(int j = 0; j < length(pairwiseid); j++)
                  outfile << array_threaded[i][j] << "\t";
               outfile << endl;
         }
      }   
   }
   return 0;
}
*/

int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

   try
   {
      outfile.open(toCString(options.outputFileName), std::ios_base::out);
   }
   catch (const ifstream::failure& e)
   {
      cout << "Error: could not open output file ";
      cout << toCString(options.outputFileName) << endl;
      return 1; 
   }

   options.effectiveLength = options.klen;
   if(parseMask(options, options.effectiveLength) == 1)
      return 1;

   if(options.pairwiseFileName != NULL)
   {
      // maybe I should make secondry version
      clock_t start;

/*
      start = clock();
      threaded_pw(options);
      cout << "Time: ";
      cout << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
      cout << " ms" << endl;
*/

      // new one
      start = clock();
      pairwise_matrix_test(options);
      cout << "Time: ";
      cout << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
      cout << " ms" << endl;
   }
   else if (options.referenceFileName != NULL && options.queryFileName != NULL)
   {
      //read in reference
      try 
      {
         SeqFileIn seqRefFileIn(toCString(options.referenceFileName));
         readRecords(referenceids, referenceseqs, seqRefFileIn);
      } 
      catch (IOError const & e)
      {
         cout << "Could not read in records " << endl;
         cout << e.what() << endl;
         return 1;
      }
      catch (ParseError const & e)
      {
         cout << "There is a formating error in ";
         cout << options.referenceFileName << endl << e.what() << endl;
         return 1;
      }
      catch (...)
      {
         return 1;
      }

      if(options.type == "d2s" || options.type == "hao" || 
         options.type == "d2star" || options.type == "dai")
      {
         kmer_count_map = makecomplete(options);
      }
      else if(options.type == "d2s-opt")
      {
         kmer_count_map = makequick(options, referenceseqs);
      }

      if(options.debug == true)
         cout << " Getting ready to read reference" << endl;

      //read reference and counts into memory
      for(int i = 0; i < length(referenceids); i++)
      {
         String<AminoAcid> seq;
         if(options.noreverse == false)
            seq = doRevCompl(referenceseqs[i]);
         else
            seq = referenceseqs[i];

         // ensure counts are done if mask is done
         if(options.mask.size() > 0)
            ref_counts_vec.push_back(count(seq, options.klen, options.mask));
         else
            ref_counts_vec.push_back(count(seq, options.klen));

         if(options.type == "d2s" || options.type == "hao" || 
            options.type == "d2star" || options.type == "d2s-opt" || 
            options.type == "dai")
         {
            ref_markov_vec.push_back(markov(options.effectiveLength, seq,
                                     options.markovOrder, kmer_count_map));
         }
      }

      if(!open(queryFileIn, (toCString(options.queryFileName))))
      {
         cerr << "Error: could not open file ";
         cerr << toCString(options.queryFileName) << endl;
         return 1;
      }

      int arraySize = options.num_threads;
      thread * workers = new thread[arraySize];

      for(int w = 0; w < options.num_threads; w++)
      {
         workers[w] = thread(mainloop, options);
      }

      for(int w = 0; w < options.num_threads; w++)
      {
         workers[w].join();
      }
   }
   return 0;
}
