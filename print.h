#ifndef PRINT_H
#define PRINT_H

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
};

template <typename TAlphabet>
int printPhylyp(ModifyStringOptions options,
                StringSet<String<TAlphabet>> const & pwseqs,
                vector<vector<double>> const & results)
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
      cout << length(pwseqs) << endl;
   }
   else
   {
      outfile_new << length(pwseqs) << endl;
   }

   for(int i = 0; i < length(pwseqs); i++)
   {
      StringSet<CharString> split;
      strSplit(split, pwseqs[i]);
      int cutsize = 10;
      CharString qName = split[0];
      qName = namecut(qName, cutsize);

      if(options.outputFileName == NULL)
      {
         cout << qName << "\t";
         for(int j = 0; j < length(pwseqs); j++)
            cout << results[i][j] << "\t";
         cout << endl;
      }
      else
      {
         outfile_new << qName << "\t";
         for(int j = 0; j < length(pwseqs); j++)
            outfile_new << results[i][j] << "\t";
         outfile_new << endl;
      }
   }

   close(outfile_new);
   return 0;
};

#endif
