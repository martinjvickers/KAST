
template <typename TAlphabet>
int printResult(ModifyStringOptions options, CharString &queryid,
                ofstream &outfile, String<TAlphabet> &queryseq,
                map<double, int> &results, StringSet<CharString> &referenceids,
                StringSet<String<TAlphabet>> &referenceseqs)
{
      if(options.output_format == "tabular")
      {
         StringSet<CharString> split;
         strSplit(split, queryid);
         CharString qName = split[0];
         if(options.outputFileName == NULL)
         {
            cout << "############################ ";
            cout << length(queryseq) << "\t" << gc_ratio(queryseq);
            cout << "\t" << qName << endl;
         }
         else
         {
            outfile << "############################ ";
            outfile << length(queryseq) << "\t" << gc_ratio(queryseq);
            outfile << "\t" << qName << endl;
         }

         for(pair<double, int> p: results)
         {
            StringSet<CharString> split2;
            strSplit(split2, referenceids[p.second]);
            if(options.outputFileName == NULL)
            {
               cout << p.first << "\t" << length(referenceseqs[p.second]);
               cout << "\t" << gc_ratio(referenceseqs[p.second]);
               cout << "\t" << split2[0] << endl;
            }
            else
            {
               outfile << p.first << "\t" << length(referenceseqs[p.second]);
               outfile << "\t" << gc_ratio(referenceseqs[p.second]);
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
               cout << gc_ratio(referenceseqs[p.second]) << "\t";
               cout << gc_ratio(queryseq) << "\t" << count << "\t";
               cout << p.first << endl;
            }
            else
            {
               outfile << referenceids[p.second] << "\t" << queryid << "\t";
               outfile << length(referenceseqs[p.second]) << "\t";
               outfile << length(queryseq) << "\t";
               outfile << gc_ratio(referenceseqs[p.second]) << "\t";
               outfile << gc_ratio(queryseq) << "\t" << count << "\t";
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
int create_ref_db(StringSet<String<TAlphabet>> &refseqs, 
                  vector<map<String<TAlphabet>, unsigned int>> &refcounts, 
                  vector<map<String<TAlphabet>, double>> &refmarkov,
                  unsigned int klen, bool noreverse, CharString distType,
                  unsigned int markovOrder, vector<String<TAlphabet>> &allKmers)
{

   // let's just do this sequentially to begin with and we 
   // can sort out threading afterwards.

   // go through reference and count
   for(unsigned int i = 0; i < length(refseqs); ++i)
   {
      refcounts.push_back(count_test(refseqs[i], klen, noreverse));

      if(distType == "d2s" || distType == "hao" ||
         distType == "d2star" || distType == "dai")
      {
         refmarkov.push_back(markov_test(klen, refseqs[i], markovOrder, allKmers, noreverse));
      }
      else if(distType == "D2Star" || distType == "D2S")
      {
         refmarkov.push_back(markov_old(klen, refseqs[i], markovOrder, allKmers, noreverse));
      }
   }

   return 0; 
}

template <typename TAlphabet>
int search_thread(ModifyStringOptions options,
                  vector<map<String<TAlphabet>, unsigned int>> &refcounts,
                  vector<map<String<TAlphabet>, double>> &refmarkov,
                  mutex &read, mutex &write, SeqFileIn &qrySeqFileIn,
                  vector<String<TAlphabet>> &allKmers,
                  StringSet<CharString> &refids,
                  StringSet<String<TAlphabet>> &refseqs,
                  ofstream &outfile)
{

   while(1)
   {
      CharString tempseq;
      CharString queryid;

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

      // Read in the sequence as a CharString and then convert to our
      // TAlphabet. This is because readRecord doesn't seem to like
      // reducedAminoAcid
      String<TAlphabet> queryseq = tempseq;

      map<String<TAlphabet>, unsigned int> querycounts;
      querycounts = count_test(queryseq, options.klen, options.noreverse);

      map<String<TAlphabet>, double> qrymarkov;

      // I now need to search the db to find the best match.
      if(options.type == "d2s" || options.type == "hao" ||
         options.type == "d2star" || options.type == "dai")
      {
         qrymarkov = markov_test(options.klen, queryseq, options.markovOrder, allKmers, options.noreverse);

      }
      else if(options.type == "D2Star" || options.type == "D2S")
      {
         qrymarkov = markov_old(options.klen, queryseq, options.markovOrder, allKmers, options.noreverse);
      }

      // to store the results 
      map<double, int> results;

      // here we can go about search for the nearest
      if((refcounts.size() != refmarkov.size()) && (options.type == "d2s" || options.type == "hao" ||
         options.type == "d2star" || options.type == "dai" || options.type == "D2Star" || options.type == "D2S"))
      {
         // should be the same size
         cerr << "ERROR: sizes are not the same " << endl;
         return 1;
      }

      for(unsigned int i = 0; i < refcounts.size(); i++)
      {
         double dist;

         if(options.type == "kmer")
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
            dist = d2s(allKmers, querycounts, qrymarkov, refcounts[i], refmarkov[i]);
         else if(options.type == "hao")
            dist = hao(allKmers, querycounts, qrymarkov, refcounts[i], refmarkov[i]);
         else if(options.type == "d2star" || options.type == "D2Star")
            dist = d2star(allKmers, querycounts, qrymarkov, refcounts[i], refmarkov[i]);
         else if(options.type == "dai")
            dist = dai(allKmers, querycounts, qrymarkov, refcounts[i], refmarkov[i]);

         // stores the smallest distance results and corresponding location in refSeq
         results.insert(pair<double, int> (dist, i));
         if(results.size() > options.nohits)
         {
            map<double, int>::iterator it = results.end();
            results.erase(--it);
         }
      }

      write.lock();
      printResult(options, queryid, outfile, queryseq, results, refids, refseqs);
      write.unlock();
   }

   return 0;
}

template <typename TAlphabet>
int search_db(ModifyStringOptions options,
              vector<map<String<TAlphabet>, unsigned int>> &refcounts, 
              vector<map<String<TAlphabet>, double>> &refmarkov,
              vector<String<TAlphabet>> &allKmers,
              StringSet<CharString> &refids,
              StringSet<String<TAlphabet>> &refseqs)
{
   // open up file
   SeqFileIn qrySeqFileIn;
   if(!open(qrySeqFileIn, (toCString(options.queryFileName))))
   {
      cerr << "Error: could not open query file ";
      cerr << toCString(options.queryFileName) << endl;
      return 1;
   }

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
   {
      vectorOfSeqs.push_back(thread(search_thread<TAlphabet>, options, ref(refcounts), ref(refmarkov), ref(read), ref(write), ref(qrySeqFileIn), ref(allKmers), ref(refids), ref(refseqs), ref(outfile) ));
   }

   for(auto &thread : vectorOfSeqs)
   {
      thread.join();
   }

   return 0;
}

// This is the main entry point
template <typename TAlphabet>
int query_ref_search(ModifyStringOptions options, TAlphabet const & alphabetType)
{
   cout << "Performing search based mode" << endl;

   // open up file
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
      CharString tmprefseq, refid;
      readRecord(refid, tmprefseq, refSeqFileIn);
      appendValue(refids, refid);
      appendValue(refseqs, tmprefseq);
   }

   cout << "Reference read in " << length(refseqs) << endl;

   // will store the counts/bg_model
   vector<map<String<TAlphabet>, unsigned int>> refcounts;
   vector<map<String<TAlphabet>, double>> refmarkov;

   vector<String<TAlphabet>> allKmers;
   if(options.type == "d2s" || options.type == "hao" ||
      options.type == "d2star" || options.type == "dai" ||
      options.type == "D2Star" || options.type == "D2S")
   {
      allKmers = makecomplete(options.klen, alphabetType);
   }

   // create reference database of counts
   create_ref_db(refseqs, refcounts, refmarkov, options.klen, options.noreverse, options.type, options.markovOrder, allKmers);

   // search
   search_db(options, refcounts, refmarkov, allKmers, refids, refseqs);

   return 0;
}
