
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
   }

   return 0; 
}

template <typename TAlphabet>
int search_thread(ModifyStringOptions options,
                  vector<map<String<TAlphabet>, unsigned int>> &refcounts,
                  vector<map<String<TAlphabet>, double>> &refmarkov,
                  mutex &read, mutex &write, SeqFileIn &qrySeqFileIn,
                  vector<String<TAlphabet>> &allKmers)
{

   while(1)
   {
      String<TAlphabet> queryseq;
      CharString queryid;

      read.lock();
      if(!atEnd(qrySeqFileIn))
      {
         readRecord(queryid, queryseq, qrySeqFileIn);
      }
      else
      {
         read.unlock();
         return 0;
      }
      read.unlock();

      map<String<TAlphabet>, unsigned int> querycounts;
      querycounts = count_test(queryseq, options.klen, options.noreverse);

      map<String<TAlphabet>, double> qrymarkov;

      // I now need to search the db to find the best match.
      if(options.type == "d2s" || options.type == "hao" ||
         options.type == "d2star" || options.type == "dai")
      {
         qrymarkov = markov_test(options.klen, queryseq, options.markovOrder, allKmers, options.noreverse);

      }

      // here we can go about search for the nearest
      if(refcounts.size() != refmarkov.size())
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
         else if(options.type == "d2s")
            dist = d2s(allKmers, querycounts, qrymarkov, refcounts[i], refmarkov[i]);
         else if(options.type == "hao")
            dist = d2s(allKmers, querycounts, qrymarkov, refcounts[i], refmarkov[i]);
         else if(options.type == "d2star")
            dist = d2s(allKmers, querycounts, qrymarkov, refcounts[i], refmarkov[i]);
         else if(options.type == "dai")
            dist = d2s(allKmers, querycounts, qrymarkov, refcounts[i], refmarkov[i]);

      }
   }

   return 0;
}

template <typename TAlphabet>
int search_db(ModifyStringOptions options,
              vector<map<String<TAlphabet>, unsigned int>> &refcounts, 
              vector<map<String<TAlphabet>, double>> &refmarkov,
              vector<String<TAlphabet>> &allKmers)
{
   // open up file
   SeqFileIn qrySeqFileIn;
   if(!open(qrySeqFileIn, (toCString(options.queryFileName))))
   {
      cerr << "Error: could not open query file ";
      cerr << toCString(options.queryFileName) << endl;
      return 1;
   }

   mutex read, write;
   vector<thread> vectorOfSeqs;

   for(unsigned i = 0; i < options.num_threads; i++)
   {
      vectorOfSeqs.push_back(thread(search_thread<TAlphabet>, options, ref(refcounts), ref(refmarkov), ref(read), ref(write), ref(qrySeqFileIn), ref(allKmers) ));
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

   // if it's a markov method, just make the complete map, it'll be easier

   // read sequences into RAM
   StringSet<CharString> refids;
   StringSet<String<TAlphabet>> refseqs;
   try
   {
      SeqFileIn refSeqFileIn(toCString(options.referenceFileName));
      readRecords(refids, refseqs, refSeqFileIn);
   }
   catch (...)
   {
      return 1;
   }

   // will store the counts/bg_model
   vector<map<String<TAlphabet>, unsigned int>> refcounts;
   vector<map<String<TAlphabet>, double>> refmarkov;

   vector<String<TAlphabet>> allKmers;
   if(options.type == "d2s" || options.type == "hao" ||
      options.type == "d2star" || options.type == "dai")
   {
      allKmers = makecomplete(options.klen, alphabetType);
   }

   // create reference database of counts
   create_ref_db(refseqs, refcounts, refmarkov, options.klen, options.noreverse, options.type, options.markovOrder, allKmers);

   // search
   search_db(options, refcounts, refmarkov, allKmers);

   // print result

   return 0;
}
