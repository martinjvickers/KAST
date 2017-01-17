#include <string>
#include <unordered_map>
//#include "utils.cpp"
//#include "common.h"
class Seq {

	public:
		Seq(IupacString seq, CharString seqid, bool noreverse, int klen, int markovOrder); // constructor
		~Seq(); // destructor
		IupacString getSeq(void); // return sequence
		CharString getID(void);
		int getTotalCounts(int klen); //return total counts
		int getTotalMarkov(int klen, int markov);
		double getSumProb(int klen, int markov);
		double missing(Dna5String kmer, int markovOrder);
		unordered_map<string, long long int> getCounts(int klen);
		unordered_map<string,markov_dat> getMarkov(int klen, int markov);
		//printers
		void printCounts(int klen);
		void printMarkov(int klen, int markov);
	private:
		//functions
		count_obj count(int klen);
		markov_obj markov(int klen, int markovOrder);
		Iupac getRevCompl(Iupac const & nucleotide);
                Dna5String doRevCompl(Dna5String seq);
		//values
		unordered_map<string, long long int> kmer_count_map;
		unordered_map<string, markov_dat> markovmap;
		unordered_map<string, long long int> r_count;
		int r_total = -1;
		unordered_map<string, long long int> r1_count;
		int r1_total = -1;
		int markov_total = -1;
		double markov_sum = -1;
		IupacString sequence;
		CharString id;
		IupacString orig_sequence;
		int tot = -1; //if -1 it means it's never been initialised
	protected:
};

//PUBLIC

//constructor, requires the sequence
Seq::Seq(IupacString seq, CharString seqid, bool noreverse, int klen, int markovOrder)
{
	id = seqid;
	orig_sequence = seq;
	if(noreverse==1)
	{
		sequence = seq;
	} else 
	{
		sequence = Seq::doRevCompl(seq);
	}
	id = seqid;
	count_obj obj = Seq::count(klen);
	tot = obj.total;
	kmer_count_map = obj.kmer_counts;
	/*markov_obj obj = Seq::markov(klen, markovOrder);
	markovmap = obj.markov_counts;
	markov_total = obj.total_count;
	markov_sum = obj.sum_prob;*/
}

//destructor
Seq::~Seq(void)
{
}

//returns the sequence
IupacString Seq::getSeq(void)
{
	return sequence;
}

//returns the sequence
CharString Seq::getID(void)
{
        return id;
}

//returns the total count
int Seq::getTotalCounts(int klen)
{
	//you can't return total counts if you've not done the counting
	if(tot < 0)
	{
	//	cout << "I'm doing this the first time for " << id << endl;
		count_obj obj = Seq::count(klen);
		tot = obj.total;
		kmer_count_map = obj.kmer_counts;
	} else {
	//	cout << "I've already done this for " << id << endl;
	}

	return tot;
}

int Seq::getTotalMarkov(int klen, int markov)
{
	//check to see if we need to calculate counts
        if(kmer_count_map.empty())
        {
                count_obj obj = Seq::count(klen);
                tot = obj.total;
                kmer_count_map = obj.kmer_counts;
        }

        if(markovmap.empty())
        {
                markov_obj obj = Seq::markov(klen, markov);
                markovmap = obj.markov_counts;
                markov_total = obj.total_count;
                markov_sum = obj.sum_prob;
        }
        return markov_total;
}

double Seq::getSumProb(int klen, int markov)
{
	//check to see if we need to calculate counts
        if(kmer_count_map.empty())
        {
                count_obj obj = Seq::count(klen);
                tot = obj.total;
                kmer_count_map = obj.kmer_counts;
        }

        if(markovmap.empty())
        {
                count_obj r_obj = Seq::count(markov);
                r_count = r_obj.kmer_counts;
		r_total  = r_obj.total;
                count_obj r1_obj = Seq::count(markov+1);
                r1_count = r1_obj.kmer_counts;
		r1_total  = r1_obj.total;
                markov_obj obj = Seq::markov(klen, markov);
                markovmap = obj.markov_counts;
                markov_total = obj.total_count;
                markov_sum = obj.sum_prob;
        }
        return markov_sum;
}

unordered_map<string, long long int> Seq::getCounts(int klen)
{
	//check to see if we need to calculate this
	if(kmer_count_map.empty())
	{
		count_obj obj = Seq::count(klen);
                tot = obj.total;
                kmer_count_map = obj.kmer_counts;
	}
	return kmer_count_map;
}

unordered_map<string, markov_dat> Seq::getMarkov(int klen, int markov)
{
        //check to see if we need to calculate counts
        if(kmer_count_map.empty())
        {
                count_obj obj = Seq::count(klen);
                tot = obj.total;
                kmer_count_map = obj.kmer_counts;
        }

	if(markovmap.empty())
	{
		count_obj r_obj = Seq::count(markov);
		//count_obj r_obj = count(markov, orig_sequence);
		r_count = r_obj.kmer_counts;
		r_total  = r_obj.total;
		count_obj r1_obj = Seq::count(markov+1);
		r1_count = r1_obj.kmer_counts;
		r1_total  = r1_obj.total;
		markov_obj obj = Seq::markov(klen, markov);
                markovmap = obj.markov_counts;
		markov_total = obj.total_count;
                markov_sum = obj.sum_prob;
	}
        return markovmap;
}

void Seq::printCounts(int klen)
{
	//if not already calculated
	if(kmer_count_map.empty())
	{
		count_obj obj = Seq::count(klen);
                tot = obj.total;
                kmer_count_map = obj.kmer_counts;
        }

	//print         
	for(auto const& p: kmer_count_map)
        {
      	        cout << p.first << " " << p.second << endl;
        }
}

void Seq::printMarkov(int klen, int markov)
{
        //if not already calculated
        if(kmer_count_map.empty())
        {
	//	cout << "kmer_count_map.empty()" << endl;
                count_obj obj = Seq::count(klen);
                tot = obj.total;
                kmer_count_map = obj.kmer_counts;
        }

	if(markovmap.empty())
        {
	//	cout << "markovmap.empty()" << endl;
                count_obj r_obj = Seq::count(markov);
                r_count = r_obj.kmer_counts;
		r_total  = r_obj.total;
                count_obj r1_obj = Seq::count(markov+1);
                r1_count = r1_obj.kmer_counts;
		r1_total = r1_obj.total;
		markov_obj obj = Seq::markov(klen, markov);
		markovmap = obj.markov_counts;
		markov_total = obj.total_count;
                markov_sum = obj.sum_prob;
        }

        //print         
        for(auto const& p: markovmap)
        {
                cout << p.first << " " << p.second.count << " " << p.second.prob <<  endl;
        }
	for(auto const& p: r_count)
	{
		cout << p.first << " " << p.second << endl;
	}
}

//PRIVATE

//calculate the counts
//BUG!!! I need a check that if the klen is longer than the sequence, it should fall out gracefully
count_obj Seq::count(int klen)
{
	int total = 0;
        unordered_map<string, long long int> map;

	/*
	if(klen==1)
	{
		cout << "testing count k=1 " << endl;
	}
	*/

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
		/*	if(klen==1)
			{
				cout << kmer << " " << kmer.size() << " " << count << endl;
			}
		*/
                        total++;
                }
        }

	count_obj obj;
	obj.kmer_counts = map;
	obj.total = total;
	return obj;
}

double Seq::missing(Dna5String kmer, int markovOrder)
{
	count_obj markovcounts = Seq::count(markovOrder);
	double prob = 1.0;
	double tot = markovcounts.total;
	for(int l = 0; l < (length(kmer)); l++)
	{
		int j = l + markovOrder;
                string inf;
                assign(inf,infix(kmer, l, j));
                prob = prob * (markovcounts.kmer_counts[inf] / tot);
	}
	return prob;
}

markov_obj Seq::markov(int klen, int markovOrder)
{
	unordered_map<string,markov_dat> markovmap;
	double sum_prob = 0.0;
	int total_count = 0;

	count_obj markovcounts = Seq::count(markovOrder);
	double tot = markovcounts.total;

	for(pair<string, long long int> p: kmer_count_map)
        {
		double prob = 1.0;
		Dna5String kmer = p.first;

		for(int l = 0; l < (length(kmer)); l++)
                {
			int j = l + markovOrder;
        	        string inf;
			assign(inf,infix(kmer, l, j));
			prob = prob * (markovcounts.kmer_counts[inf] / tot);
		}
		markov_dat dat;
		dat.count = p.second;
		dat.prob = prob;
		sum_prob = sum_prob + prob;
		markovmap[p.first] = dat;
	}

	markov_obj obj;
	obj.markov_counts = markovmap;
	obj.total_count = total_count;
	obj.sum_prob = sum_prob;
	return obj;

}

/**/
Iupac Seq::getRevCompl(Iupac const & nucleotide)
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
Dna5String Seq::doRevCompl(Dna5String seq)
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
