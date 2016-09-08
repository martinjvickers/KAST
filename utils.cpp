/*
MIT License

Copyright (c) 2016 Martin James Vickers

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

#include "common.h"
#include "utils.h"

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

/*

*/
void count(Dna5String seq, int klen, unordered_map<string, long long int> & map)
{

	//iterate over the sequence
        for(int i = 0; i <= length(seq)-klen; i++)
        {
		//get our kmer
                string kmer;
                assign(kmer,infix(seq, i, i+klen));

		//need to drop if there is an N in it
		size_t found = kmer.find("N");
		if(found>kmer.size()){
			long long int count = map[kmer];
			map[kmer] = count + 1;
		}
        }
}

/*

*/
void markov(Dna5String seq, int klen, int markovOrder, unordered_map<string,markov_dat> & markovmap)
{

	//get out markov counts
	int newmarkov = markovOrder + 1;
	unordered_map<string, long long int> markovcounts;
	count(seq, newmarkov, markovcounts);

	//get out regular kmer counts
	unordered_map<string, long long int> kmercounts;
	count(seq, klen, kmercounts);

	//calculate total number of markov kmers
	int sumCounts = 0;
	for(pair<string, long long int> p: markovcounts)
	{
		sumCounts = sumCounts + p.second;
	}

	//calculate frequency of kmers
	unordered_map<string, double> markovFreq;
	for(pair<string, long long int> p: markovcounts)
        {
		markovFreq[p.first] = (double)p.second / (double)sumCounts;
        }

	for(pair<string, long long int> p: kmercounts)
	{
		double prob = 1.0;
		Dna5String kmer = p.first;
		
		for(int l = 0; l < length(kmer); l++)
                {
			int j = l + newmarkov;
			string inf;
			assign(inf,infix(kmer, l, j));

			//this is where i need to decide the behaviour when nothing is returned
			double freq = markovFreq[inf];
                        if(freq != 0)
                        {
				prob = prob * freq;
                        } else {
                        }
	
			if(j == length(kmer))
			{
				break;
			}

		}
		markov_dat dat;
		dat.count = p.second;
		dat.prob = prob;
		markovmap[p.first] = dat;
	}
}

//I want it to return an array with element 0 is the smallest
void recordall(int nohits, double hits[], double value, int seqcurrpos, int hitpos[])
{
        //if value is smaller than the current largest, lets knock it off
        if(value < hits[nohits-1]){

                hits[nohits-1] = value;
                hitpos[nohits-1] = seqcurrpos; //copy position move
                int j;
                double temp;
                int temppos;

                //now, lets see if our new value is smaller than any of the others
                //iterate through the rest of the elements
                int i = nohits-1;
                while(i >= 0 && hits[i] < hits[i-1])
                {
                        temp = hits[i-1];
                        temppos = hitpos[i-1];

                        hits[i-1] = hits[i];
                        hitpos[i-1] = hitpos[i];

                        hits[i] = temp;
                        hitpos[i] = temppos;

                        i--;
                }
        }
}
