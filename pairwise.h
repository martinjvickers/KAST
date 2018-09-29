
template <typename TAlphabet>
int count_threads(ModifyStringOptions options, SeqFileIn &pairwiseFileIn,
                  vector<pair<CharString, map<String<TAlphabet>, unsigned int>>> &pw_counts,
                  mutex &read, mutex &write,
                  vector<pair<CharString, map<String<TAlphabet>, double>>> &mv_counts,
                  vector<String<TAlphabet>> &kmer_count_map)
{
   CharString pwid;
   CharString tmpseq;

   while(1)
   {
      read.lock();
      if(!atEnd(pairwiseFileIn))
      {
         readRecord(pwid, tmpseq, pairwiseFileIn);
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
      String<TAlphabet> refseq = tmpseq;

      pair< CharString, map<String<TAlphabet>, unsigned int> > meh = make_pair(pwid, count_test(refseq, options.klen, options.noreverse));
      pair< CharString, map<String<TAlphabet>, double>> meh_markov;

      // If markov type, calculate markov counts
      if(options.type == "d2s" || options.type == "hao" ||
         options.type == "d2star" || options.type == "dai")
      {
         meh_markov = make_pair(pwid, markov_test(options.klen, refseq, options.markovOrder, kmer_count_map, options.noreverse));
      }
      else if(options.type == "D2S" || options.type == "D2Star")
      {
         meh_markov = make_pair(pwid, markov_old(options.klen, refseq, options.markovOrder, kmer_count_map, options.noreverse));
      }

      // lock both because we want both the markov and counts to be in the same position
      // in their respective vectors
      write.lock();
      pw_counts.push_back(meh);
      if(options.type == "d2s" || options.type == "hao" ||
         options.type == "d2star" || options.type == "dai" ||
         options.type == "D2S" || options.type == "D2Star")
      {
         mv_counts.push_back(meh_markov);
      }
      write.unlock();

   }

}


template <typename TAlphabet>
int distance_thread(vector<pair<CharString, map<String<TAlphabet>, unsigned int>>> &pw_counts, 
                    unsigned &rI, unsigned &cI, vector< vector<double> > & array_threaded_internal, 
                    ModifyStringOptions options, mutex &location, 
                    vector<pair<CharString, map<String<TAlphabet>, double>>> &mv_counts, 
                    vector<String<TAlphabet>> &kmer_count_map)
{
   unsigned row, column;

   while(1)
   {
      location.lock();
      row = rI;
      column = cI;

      if(cI < pw_counts.size())
      {
         cI++;
      }
      else if(cI >= pw_counts.size() && rI < pw_counts.size())
      {
         cI = rI;
         rI++;
         location.unlock();
         continue;
      }
      else if(!(rI < pw_counts.size()))
      {
         location.unlock();
         return 0;
      }
      location.unlock();

      if(row < pw_counts.size() && column < pw_counts.size())
      {
         map<String<TAlphabet>, unsigned int> q1 = pw_counts[column].second;
         map<String<TAlphabet>, unsigned int> q2 = pw_counts[row].second;

         map<String<TAlphabet>, double> m1;
         map<String<TAlphabet>, double> m2;

         if(options.type == "d2s" || options.type == "hao" ||
            options.type == "d2star" || options.type == "dai" ||
            options.type == "D2S" || options.type == "D2Star")
         {
            m1 = mv_counts[column].second;
            m2 = mv_counts[row].second;
         }

         double dist;

         if(options.type == "euclid")
            dist = euler(q1, q2);
         else if(options.type == "d2")
            dist = d2(q1, q2);
         else if(options.type == "manhattan")
            dist = manhattan(q1, q2);
         else if(options.type == "chebyshev")
            dist = chebyshev(q1, q2);
         else if(options.type == "bc")
            dist = bray_curtis_distance(q1, q2);
         else if(options.type == "ngd")
            dist = normalised_google_distance(q1, q2);
         else if(options.type == "d2s" || options.type == "D2S")
            dist = d2s(kmer_count_map, q1, m1, q2, m2);
         else if(options.type == "hao")
            dist = hao(kmer_count_map, q1, m1, q2, m2);
         else if(options.type == "d2star" || options.type == "D2Star")
            dist = d2star(kmer_count_map, q1, m1, q2, m2);
         else if(options.type == "dai")
            dist = dai(kmer_count_map, q1, m1, q2, m2);

         array_threaded_internal[row][column] = dist;
         array_threaded_internal[column][row] = dist;

      }
   }
}

template <typename TAlphabet>
int pairwise_matrix(ModifyStringOptions options, TAlphabet const & alphabetType)
{
   // Read in the input file
   SeqFileIn pairwiseFileIn;
   CharString pwid;
   String<TAlphabet> pwseq;
   vector< vector<double> > array_threaded_internal; // stores the distance matrix
   vector<pair<CharString, map<String<TAlphabet>, unsigned int>>> pw_counts;

   // things for markov distances
   vector<String<TAlphabet>> kmer_count_map;
   vector<pair<CharString, map<String<TAlphabet>, double>>> mv_counts;

   // Open fasta/fastq file
   if(!open(pairwiseFileIn, (toCString(options.pairwiseFileName))))
   {
      cerr << "Error: could not open file ";
      cerr << toCString(options.pairwiseFileName) << endl;
      return 1;
   }

   // Create thread vector and mutex's
   mutex read, write;
   vector<thread> vectorOfSeqs;
   unsigned int n = std::thread::hardware_concurrency();
   unsigned int cores = options.num_threads;

   // create kmer_count_map if doing a markov model
   if(options.type == "d2s" || options.type == "hao" ||
      options.type == "d2star" || options.type == "dai" ||
      options.type == "D2Star" || options.type == "D2S")
   {
      kmer_count_map = makecomplete(options.klen, alphabetType);
   }

   cout << "About to count" << endl;

   // thread to count // here, I think, if it's markov based, the thread also does the markov count
   for(unsigned i = 0; i < cores; i++)
   {
      vectorOfSeqs.push_back(thread(count_threads<TAlphabet>, options, ref(pairwiseFileIn), 
                                    ref(pw_counts), ref(read), ref(write), ref(mv_counts), 
                                    ref(kmer_count_map) ));
   }

   for(auto &thread : vectorOfSeqs)
   {
      thread.join();
   }

   close(pairwiseFileIn);

   // Store the results 
   array_threaded_internal.resize(pw_counts.size(),
                                  vector<double>(pw_counts.size(), 0.0));

   // Calculate the distances
   vector<thread> vectorOfThreads;
   mutex location;
   unsigned rI = 0;
   unsigned cI = 0;

   for(unsigned i = 0; i < cores; i++)
   {
      vectorOfThreads.push_back(thread(distance_thread<TAlphabet>, ref(pw_counts), ref(rI), 
                                       ref(cI), std::ref(array_threaded_internal), options, 
                                       ref(location), ref(mv_counts), ref(kmer_count_map) ));
   }
   for(auto &thread : vectorOfThreads)
   {
      thread.join();
   }

   // Print results
   printPhylyp(options, pw_counts, array_threaded_internal);

   return 0;

}

template <typename TAlphabet>
int pairwise_all_matrix(ModifyStringOptions options, TAlphabet const & alphabetType)
{
   SeqFileIn pairwiseFileIn;
   StringSet<CharString> pwids;
   StringSet<String<TAlphabet>> pwseqs;
   
   // Open fasta/fastq file
   if(!open(pairwiseFileIn, (toCString(options.pairwiseFileName))))
   {
      cerr << "Error: could not open file ";
      cerr << toCString(options.pairwiseFileName) << endl;
      return 1;
   }
   
   readRecords(pwids, pwseqs, pairwiseFileIn);

   vector<String<TAlphabet>> kmer_count_map = makecomplete(options.klen, alphabetType);

   cout << "Q1\tQ2\tEuler\td2\tManhattan\tChebyshev\tBC\tNGD\td2s\tHao\td2Star\tdai\n";

   // let's go through every comparison
   for(unsigned int i = 0; i < length(pwids); i++)
   {
      for(unsigned int j = 0; j < length(pwids); j++)
      {
         map<String<TAlphabet>, unsigned int> count_p1 = count_test(pwseqs[i], options.klen, options.noreverse);
         map<String<TAlphabet>, unsigned int> count_p2 = count_test(pwseqs[j], options.klen, options.noreverse);

         map<String<TAlphabet>, double> markov_p1 = markov_test(options.klen, pwseqs[i], options.markovOrder, 
                                                                 kmer_count_map, options.noreverse);
         map<String<TAlphabet>, double> markov_p2 = markov_test(options.klen, pwseqs[j], options.markovOrder, 
                                                                 kmer_count_map, options.noreverse);

         cout << pwids[i] << "\t" << pwids[j] << "\t" << euler(count_p1, count_p2) << "\t";
         cout << d2(count_p1, count_p2) << "\t" << manhattan(count_p1, count_p2) << "\t";
         cout << bray_curtis_distance(count_p1, count_p2) << "\t";
         cout << normalised_google_distance(count_p1, count_p2) << "\t";
         cout << d2s(kmer_count_map, count_p1, markov_p1, count_p2, markov_p2) << "\t";
         cout << hao(kmer_count_map, count_p1, markov_p1, count_p2, markov_p2) << "\t";
         cout << d2star(kmer_count_map, count_p1, markov_p1, count_p2, markov_p2) << "\t";
         cout << dai(kmer_count_map, count_p1, markov_p1, count_p2, markov_p2) << "\t";

         cout << endl;
         
      }
   }
}

