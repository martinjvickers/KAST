#include <unordered_map>

/**/
Iupac getRevCompl(Iupac const & nucleotide)
{
        if (nucleotide == 'A')
                return 'T';
        if (nucleotide == 'T')
                return 'A';
        if (nucleotide == 'C')
                return 'G';
        if (nucleotide == 'G')
                return 'C';
        return 'N';
}

/**/
Dna5String doRevCompl(Dna5String seq)
{
        Dna5String allSeq;
        append(allSeq,seq);
        allSeq += "NNN";
        Dna5String revComplGenome;
        resize(revComplGenome, length(seq));
        for (unsigned i = 0; i < length(seq); ++i)
        {
                revComplGenome[length(seq) - 1 - i] = getRevCompl(seq[i]);
        }
        allSeq += revComplGenome;
        return allSeq;
}

//calculate the counts
//BUG!!! I need a check that if the klen is longer than the sequence, it should fall out gracefully
count_obj count(int klen, IupacString sequence)
{
        int total = 0;
        unordered_map<string, long long int> map;

        //iterate over the sequence
        for(int i = 0; i <= length(sequence)-klen; i++)
        {
                //get our kmer
                string kmer;
                assign(kmer,infix(sequence, i, i+klen));

                //need to drop if there is an N in it
                size_t found = kmer.find("N");

                if(found > kmer.size()){
                        long long int count = map[kmer];
                        map[kmer] = count + 1;
                        total++;
                }
        }
        count_obj obj;
        obj.kmer_counts = map;
        obj.total = total;
        return obj;
}

/*Will need to look at good ways of speeding this up*/
map<string,bool> makeall(ModifyStringOptions options)
{
        SeqFileIn refFileIn;
        CharString refid;
        IupacString refseq;
        open(refFileIn, (toCString(options.referenceFileName)));

        map<string,bool> kmermap;

        while(!atEnd(refFileIn))
        {
                if(kmermap.size() != pow(4,options.klen))
                {
                        readRecord(refid, refseq, refFileIn);

			Dna5String sequence;
			if(options.noreverse==1)
			{
				sequence = refseq;
			} else
			{
				sequence = doRevCompl(refseq);
			}

                        //iterate over the sequence
                        for(int i = 0; i <= length(sequence)-options.klen; i++)
                        {
                                string kmer;
                                assign(kmer,infix(sequence, i, i+options.klen));

                                //need to drop if there is an N in it
                                size_t found = kmer.find("N");

                                if(found > kmer.size())
                                {
                                        kmermap[kmer] = true;
                                }
                        }
                } else {
			cout << "There should be no more kmers to find." << endl;
			break;
		}
        }

        SeqFileIn qryFileIn;
        CharString qryid;
        IupacString qryseq;
        open(qryFileIn, (toCString(options.queryFileName)));

        while(!atEnd(qryFileIn))
        {
                if(kmermap.size() != pow(4,options.klen))
                {
                        readRecord(qryid, qryseq, qryFileIn);

                        Dna5String sequence;
                        if(options.noreverse==1)
                        {
                                sequence = qryseq;
                        } else
                        {
                                sequence = doRevCompl(qryseq);
                        }


                        //iterate over the sequence
                        for(int i = 0; i <= length(sequence)-options.klen; i++)
                        {
                                string kmer;
                                assign(kmer,infix(sequence, i, i+options.klen));

                                //need to drop if there is an N in it
                                size_t found = kmer.find("N");

                                if(found > kmer.size())
                                {
                                        kmermap[kmer] = true;
                                }
                        }
                } else {
			cout << "There should be no more kmers to find " << endl;
			break;
		}
        }
        return kmermap;
}

map<string,bool> makecomplete(ModifyStringOptions options)
{
	map<string,bool> finkmers;

        String<Dna5String> bases;
        appendValue(bases, "A");
        appendValue(bases, "G");
        appendValue(bases, "C");
        appendValue(bases, "T");

        String<Dna5String> kmers;
        kmers = bases;

        for(int j = 0; j < options.klen-1; j++)
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
//	Dna5String meh = "AGCT";
//	string kmer2;
  //      assign(kmer2,infix(meh, 0, length(meh)));
//	cout << kmer2 << endl;
//	finkmers[kmer2] = true;
	return finkmers;
}

/*Record*/
void record(ModifyStringOptions options)
{

}
