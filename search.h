template <typename TAlphabet>
int printResult(ModifyStringOptions options, CharString &queryid,
                ofstream &outfile, String<TAlphabet> &queryseq,
                multimap<double, int> &results, StringSet<CharString> const & referenceids,
                StringSet<String<TAlphabet>> const & referenceseqs)
{
      if(options.output_format == "tabular")
      {
         StringSet<CharString> split;
         strSplit(split, queryid);
         CharString qName = split[0];
         if(options.outputFileName == NULL)
         {
            cout << "############################ ";
            //cout << length(queryseq) << "\t" << gc_ratio(queryseq);
            cout << "\t" << qName << endl;
         }
         else
         {
            outfile << "############################ ";
            //outfile << length(queryseq) << "\t" << gc_ratio(queryseq);
            outfile << "\t" << qName << endl;
         }

         for(pair<double, int> p: results)
         {
            StringSet<CharString> split2;
            strSplit(split2, referenceids[p.second]);
            if(options.outputFileName == NULL)
            {
               cout << p.first << "\t" << length(referenceseqs[p.second]);
               //cout << "\t" << gc_ratio(referenceseqs[p.second]);
               cout << "\t" << split2[0] << endl;
            }
            else
            {
               outfile << p.first << "\t" << length(referenceseqs[p.second]);
               //outfile << "\t" << gc_ratio(referenceseqs[p.second]);
               outfile << "\t" << split2[0] << endl;
            }
         }
      }
      else if(options.output_format == "blastlike")
      {
         if(options.outputFileName == NULL)
         {
            cout << "RefID\tQryID\tRefLen\tQryLen\t";
            cout << "RefGC\tQryGC\tHitRank\tScore" << endl;
         }
         else
         {
            outfile << "RefID\tQryID\tRefLen\tQryLen\t";
            outfile << "RefGC\tQryGC\tHitRank\tScore" << endl;
         }

         int count = 1;
         for(pair<double, int> p: results)
         {
            if(options.outputFileName == NULL)
            {
               cout << referenceids[p.second] << "\t" << queryid << "\t";
               cout << length(referenceseqs[p.second]) << "\t";
               cout << length(queryseq) << "\t";
               //cout << gc_ratio(referenceseqs[p.second]) << "\t";
               //cout << gc_ratio(queryseq) << "\t" << count << "\t";
               cout << p.first << endl;
            }
            else
            {
               outfile << referenceids[p.second] << "\t" << queryid << "\t";
               outfile << length(referenceseqs[p.second]) << "\t";
               outfile << length(queryseq) << "\t";
               //outfile << gc_ratio(referenceseqs[p.second]) << "\t";
               //outfile << gc_ratio(queryseq) << "\t" << count << "\t";
               outfile << p.first << endl;
            }
            count++;
         }

      }
      else
      {
         if(options.outputFileName == NULL)
         {
            cout << "############################ " << queryid << endl;
         }
         else
         {
            outfile << "############################ " << queryid << endl;
         }

         for(pair<double, int> p: results)
         {
            if(options.outputFileName == NULL)
            {
               cout << referenceids[p.second] << " " << p.first << endl;
            }
            else
            {
               outfile << referenceids[p.second] << " " << p.first << endl;
            }
         }
      }
}



template <typename TAlphabet>
int search_thread(ModifyStringOptions options, SeqFileIn & qrySeqFileIn,
                  StringSet<CharString> const & refids, 
                  StringSet<String<TAlphabet>> const & refseqs,
                  StringSet<String<unsigned>> const & refcounts,
                  StringSet<String<double>> const & markovcounts,
                  mutex & read, mutex & write, ofstream & outfile)
{
   while(1)
   {
      CharString queryid;
      String<TAlphabet> queryseq;
      CharString tempseq;
      read.lock();
      if(!atEnd(qrySeqFileIn))
      {
         readRecord(queryid, tempseq, qrySeqFileIn);
      }
      else
      {
         read.unlock();
         return 0;
      }
      read.unlock();

      // this gets past readRecord not being able to 
      // work with reducedAminoAcid
      queryseq = tempseq;

      String<unsigned> querycounts;
      String<double> querymarkov;

      if(options.sequenceType == "dna")
      {
         String<Dna5> qseq = queryseq;
         if(options.noreverse == false)
         {
            String<Dna5> qseqrc = queryseq;
            reverseComplement(qseqrc);
            append(qseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"); // this should probably the same size as options.klen
            append(qseq, qseqrc);
         }
         countKmers(querycounts, qseq, options.klen);

         if(options.type == "d2s" || options.type == "D2S" ||
            options.type == "d2star" || options.type == "D2Star" ||
            options.type == "hao" || options.type == "dai")
         {
            markov(querymarkov, querycounts, qseq, options.klen, options.markovOrder);
         }
      }
      else
      {
         countKmers(querycounts, queryseq, options.klen);
         if(options.type == "d2s" || options.type == "D2S" ||
            options.type == "d2star" || options.type == "D2Star" ||
            options.type == "hao" || options.type == "dai")
         {
            markov(querymarkov, querycounts, queryseq, options.klen, options.markovOrder);
         }
      }

      // store the results
      multimap<double, int> results;

      for(int i = 0; i < length(refids); i++)
      {
         double dist = 0;
         if(options.type == "euclid")
            dist = euler(querycounts, refcounts[i]);
         else if(options.type == "d2")
            dist = d2(querycounts, refcounts[i]);
         else if(options.type == "manhattan")
            dist = manhattan(querycounts, refcounts[i]);
         else if(options.type == "chebyshev")
            dist = chebyshev(querycounts, refcounts[i]);
         else if(options.type == "bc")
            dist = bray_curtis_distance(querycounts, refcounts[i]);
         else if(options.type == "ngd")
            dist = normalised_google_distance(querycounts, refcounts[i]);
         else if(options.type == "d2s" || options.type == "D2S")
            dist = d2s(querycounts, refcounts[i], querymarkov, markovcounts[i]);
         else if(options.type == "hao")
            dist = hao(querycounts, refcounts[i], querymarkov, markovcounts[i]);
         else if(options.type == "d2star" || options.type == "D2Star")
            dist = d2star(querycounts, refcounts[i], querymarkov, markovcounts[i]);
         else if(options.type == "dai")
            dist = dai(querycounts, refcounts[i], querymarkov, markovcounts[i]);

         // stores the smallest distance results and corresponding location in refSeq
         results.insert(pair<double, int> (dist, i));
         if(results.size() > options.nohits)
         {
            map<double, int>::iterator it = results.end();
            results.erase(--it);
         }
      }

      write.lock();
      printResult<TAlphabet>(options, queryid, outfile, queryseq, results, refids, refseqs);
      write.unlock();
   }

   return 0;
}

template <typename TAlphabet>
int query_ref_search(ModifyStringOptions options, TAlphabet const & alphabetType)
{
   // read the reference into memory
   SeqFileIn refSeqFileIn;
   if(!open(refSeqFileIn, (toCString(options.referenceFileName))))
   {
      cerr << "Error: could not open reference file ";
      cerr << toCString(options.referenceFileName) << endl;
      return 1;
   }

   StringSet<CharString> refids;
   StringSet<String<TAlphabet>> refseqs;
   
   while(!atEnd(refSeqFileIn))
   {
      CharString id;
      CharString seq;
      readRecord(id, seq, refSeqFileIn);
      appendValue(refids, id);
      String<TAlphabet> convseq = seq;
      appendValue(refseqs, convseq);
   }

   // check how much RAM is required to store the reference
   if(mem_check(options, length(refseqs), alphabetType) == 1)
      return 1;

   StringSet<String<unsigned> > counts;
   resize(counts, length(refseqs));

   //if markov we don't need to do this
   StringSet<String<double> > markovCounts;
   if(options.type == "d2s" || options.type == "D2S" ||
      options.type == "d2star" || options.type == "D2Star" ||
      options.type == "hao" || options.type == "dai")
      {
         resize(markovCounts, length(refseqs));
      }

   // populate the counts
   for(int i = 0; i < length(refseqs); i++)
   {

      if(options.sequenceType == "dna")
      {
         String<Dna5> seq = refseqs[i];
         if(options.noreverse == false)
         {
            String<Dna5> seqrc = refseqs[i];
            reverseComplement(seqrc);
            append(seq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"); // this should probably the same size as options.klen
            append(seq, seqrc);
         }
 
         countKmers(counts[i], seq, options.klen);

         if(options.type == "d2s" || options.type == "D2S" ||
            options.type == "d2star" || options.type == "D2Star" ||
            options.type == "hao" || options.type == "dai")
         {
            markov(markovCounts[i], counts[i], seq, options.klen, options.markovOrder);
         }
      }
      else
      {
         String<TAlphabet> seq = refseqs[i];
         countKmers(counts[i], seq, options.klen);

         if(options.type == "d2s" || options.type == "D2S" ||
            options.type == "d2star" || options.type == "D2Star" ||
            options.type == "hao" || options.type == "dai")
         {
            markov(markovCounts[i], counts[i], seq, options.klen, options.markovOrder);
         }
      }
   }

   // search thread
   // open up query file
   SeqFileIn qrySeqFileIn;
   if(!open(qrySeqFileIn, (toCString(options.queryFileName))))
   {
      cerr << "Error: could not open query file ";
      cerr << toCString(options.queryFileName) << endl;
      return 1;
   }

   // open up output file
   ofstream outfile;

   try
   {
      outfile.open(toCString(options.outputFileName), std::ios_base::out);
   }
   catch (const ifstream::failure& e)
   {
      cerr << "Error: could not open output file ";
      cerr << toCString(options.outputFileName) << endl;
      return 1;
   }


   mutex read, write;
   vector<thread> vectorOfSeqs;

   for(unsigned i = 0; i < options.num_threads; i++)
      vectorOfSeqs.push_back(thread(search_thread<TAlphabet>, options, 
                                    ref(qrySeqFileIn), ref(refids), 
                                    ref(refseqs), ref(counts), ref(markovCounts),
                                    ref(read), ref(write), ref(outfile)
                                   )
                             );

   for(auto &thread : vectorOfSeqs)
      thread.join();

   return 0;
};
