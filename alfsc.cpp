/*
 *	ALFSC - Alignment-free Sequence Comparison
 *	Version 0.1
 *	Written by Dr. Martin Vickers (mjv08@aber.ac.uk)
 *
 *	Todo:
 *	1) Ensure d2s and d2star are working as desired
 *	4) Improve performance of kmer/markovKmer lookup
 *		e.g. use hash tables
 *	5) Multithread the program, probably best bet is
 *	   to simply run each query input as a thread.
 *	6) Create nice output files for the results
 *	   	Would be nice if they were compatible 
 *	   	with other software. e.g. blastoutputs
 *	*) Memory management, need to ensure decent use of memory
 *	*) sort out input to allow a single input fasta to work as a pairwise comparison.
 */

#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include "boost/multi_array.hpp"
#include <cassert>
#include <boost/unordered_map.hpp>
#include <string>
#include <thread>
#include <mutex>
#include <chrono>

using namespace seqan;
using namespace std;

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array<double, 2> array_type2;

int num_threads;
mutex m;
queue<thread> thread_queue;
SeqFileIn queryFileIn;


//overload the SeqFileBuffer_ so that it uses Iupac String. In this way 
//the input file is checked against Iupac and any non-A/C/G/T is silently 
//converted into a N.
namespace seqan {
	template <typename TString, typename TSSetSpec, typename TSpec>
	struct SeqFileBuffer_<StringSet<TString, TSSetSpec>, TSpec>
	{
		typedef String<Iupac> Type;
	};
}

struct ModifyStringOptions
{
        unsigned klen;
        int nohits;
        int markovOrder;
        CharString type;
        bool noreverse;
        CharString queryFileName;
	CharString referenceFileName;
};

String<Dna5String> defineKmers(int kmerlength)
{
	String<Dna5String> bases;
	appendValue(bases, "A");
	appendValue(bases, "G");
	appendValue(bases, "C");
	appendValue(bases, "T");

        String<Dna5String> kmers;

	kmers = bases;

	for(int j = 0; j < kmerlength-1; j++)
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

	return kmers;
}

//this isn't working correction. What's happening?
void count(String<Dna5String> kmers, Dna5String seq, int klen, int counts[])
{
        typedef boost::unordered_map<string, int> unordered_map;
        unordered_map map;
	//puts all the kmers into an unordered map
        for(int i = 0; i < length(kmers); i++)
        {
		string meh;
		assign(meh,kmers[i]);
		map.emplace(meh, 0);
        }

        for(int i = 0; i <= length(seq)-klen; i++)
        {
		string meh;
		assign(meh,infix(seq, i, i+klen));
		map[meh]++;
	}

	for(int i = 0; i < length(kmers); i++)
	{
		string meh;
		assign(meh,kmers[i]);
		counts[i] = (int)map[meh];
	}

	//cout << "arg " << length(kmers) << endl;
}

void markov(String<Dna5String> kmers, Dna5String seq, int klen, int counts[], int markovOrder, double kmerProb[])
{
	//count the kmers using markov order
	int newmarkov = markovOrder + 1;
	String<Dna5String> markovKmers = defineKmers(newmarkov);
	
	//count the occurances of the markov kmers
	int markovKcounts [length(markovKmers)];
	count(markovKmers, seq, newmarkov, markovKcounts);

	//sum the total number of kmers
	int sumCounts = 0;
	for(int i = 0; i < length(markovKmers); i++)
		sumCounts = sumCounts + markovKcounts[i];

	//calculate frequency of kmers
	double markovFreq [length(markovKmers)];
	for(int i = 0; i < length(markovKmers); i++)
	{
		markovFreq[i] = (double)markovKcounts[i] / (double)sumCounts;
	}

	double total = 0.0;

	//for each kmer in kmers
	for(int i = 0; i < length(kmers); i++)
	{
		double prob = 1.0;
		Dna5String kmer = kmers[i];
		
		for(int l = 0; l < length(kmer); l++)
		{
			int j = l + markovOrder + 1;
			Infix<Dna5String>::Type inf = infix(kmer,l,j);

			for(int m = 0; m < length(markovKmers); m++)
			{
				if(inf == markovKmers[m])
				{
					prob = prob * markovFreq[m];
				} else {
				}
			}

			if(j == length(kmer))
			{
				break;
			}
		}
		kmerProb[i] = prob;
		total = total + prob;
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

double d2(array_type ref, int qry[], int nokmers, int val)
{
        double sumqCrC = 0.0;
        double sumqC2 = 0.0;
        double sumrC2 = 0.0;

        for(int i = 0; i < nokmers; i++){
                sumqCrC = sumqCrC + (qry[i] * ref[val][i]);
                sumqC2 = sumqC2 + (qry[i] * qry[i]);
                sumrC2 = sumrC2 + (ref[val][i] * ref[val][i]);
        }

        double score = sumqCrC / (sqrt(sumqC2) * sqrt(sumrC2));
        return 0.5*(1-score);
}

double d2(int ref[], int qry[], int nokmers)
{
        double sumqCrC = 0.0;
        double sumqC2 = 0.0;
        double sumrC2 = 0.0;

        for(int i = 0; i < nokmers; i++){
                sumqCrC = sumqCrC + (qry[i] * ref[i]);
                sumqC2 = sumqC2 + (qry[i] * qry[i]);
                sumrC2 = sumrC2 + (ref[i] * ref[i]);
        }

        double score = sumqCrC / (sqrt(sumqC2) * sqrt(sumrC2));
        return 0.5*(1-score);
}


double euler(array_type ref, int qry[], int nokmers, int val)
{
        double score = 0.0;
        double rN = 0.0;
        double qN = 0.0;

        for(int i = 0; i < nokmers; i++){
                rN = rN + ref[val][i];
                qN = qN + qry[i];
        }

        for(int i = 0; i < nokmers; i++){
                double rF = ref[val][i] / rN;
                double qF = qry[i] / qN;
                score = score + (pow((rF - qF), 2));
        }

        return pow(score, 0.5);
}

double d2s(array_type refCounts, int qryCounts[], int nokmers, array_type2 refProbs, double qryProbs[], int val)
{
        double score = 0.0;
        double D2S = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        double rN = 0.0;
        double qN = 0.0;
        double rtot = 0.0;
        double qtot = 0.0;

        for(int i = 0; i < nokmers; i++)
        {
                qN = qN + qryCounts[i];
                qtot = qtot + qryProbs[i];
                rN = rN + refCounts[val][i];
                rtot = rtot + refProbs[val][i];
        }

        for(int i = 0; i < nokmers; i++)
        {
                double qC = qryCounts[i];
                double rC = refCounts[val][i];
                double qP = qryProbs[i];
                double rP = refProbs[val][i];

                double qCt = qC - (qN * qP);
                double rCt = rC - (rN * rP);
                double dist = sqrt(( qCt * qCt )+( rCt * rCt ));

                if(dist == 0.0)
                {
                        std::cout << "Div by zero" << std::endl;
                }

                D2S = D2S + ( ( qCt * rCt ) / dist );
                sum1 = sum1 + ( ( qCt * qCt ) / dist );
                sum2 = sum2 + ( ( rCt * rCt ) / dist );
        }

        score = 0.5 * (1 - ( (D2S) / (sqrt(sum1)*sqrt(sum2)) ) );

        return score;

}

double d2star(array_type refCounts, int qryCounts[], int nokmers, array_type2 refProbs, double qryProbs[], int val)
{
        double score = 0.0;
        double D2Star = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        double qN = 0.0;
        double rN = 0.0;

        for(int i = 0; i < nokmers; i++)
        {
                qN = qN + qryCounts[i];
                rN = rN + refCounts[val][i];
        }

        for(int i = 0; i < nokmers; i++)
        {
                double qC = qryCounts[i];
                double rC = refCounts[val][i];
                double qP = qryProbs[i];
                double rP = refProbs[val][i];

                double q_np = qN * qP;
                double r_np = rN * rP;
                double qCt = qC - q_np;
                double rCt = rC - r_np;

                double dist = sqrt(q_np) * sqrt(r_np);

                if(dist == 0.0)
                {
                        std::cout << "Div by zero" << std::endl;
                }

                D2Star = D2Star + ( ( qCt * rCt ) / dist );
                sum1 = sum1 + ( ( qCt * qCt ) / q_np );
                sum2 = sum2 + ( ( rCt * rCt ) / r_np );
        }
        score = 0.5 * (1 - ( (D2Star) / (sqrt(sum1)*sqrt(sum2)) ) );
        return score;
}

Dna getRevCompl(Dna const & nucleotide)
{
	if (nucleotide == (Dna)'A')
		return (Dna)'T';
	if (nucleotide == (Dna)'T')
		return (Dna)'A';
	if (nucleotide == (Dna)'C')
		return (Dna)'G';
        if (nucleotide == (Dna)'G')
                return (Dna)'C';	
	return (Dna)'N';
}

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

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("alfsc");
	addOption(parser, seqan::ArgParseOption("k", "klen", "Kmer Length.", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "klen", "3");
	addOption(parser, seqan::ArgParseOption("q", "query-file", "Path to the query file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	addOption(parser, seqan::ArgParseOption("r", "reference-file", "Path to the reference file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "query-file");
	addOption(parser, seqan::ArgParseOption("m", "markov-order", "Markov Order", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "markov-order", "1");
	addOption(parser, seqan::ArgParseOption("n", "num-hits", "Number of top hits to return", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "num-hits", "10");
	addOption(parser, seqan::ArgParseOption("t", "distance-type", "The method of calculating the distance between two sequences.", seqan::ArgParseArgument::STRING, "STR"));
	setValidValues(parser, "distance-type", "d2 kmer d2s d2star");
	setDefaultValue(parser, "distance-type", "d2");
	addOption(parser, seqan::ArgParseOption("nr", "no-reverse", "Do not use reverse compliment."));
	setShortDescription(parser, "Alignment-free sequence comparison.");
	setVersion(parser, "0.1");
	setDate(parser, "September 2015");
	addUsageLine(parser, "-q query.fasta -r reference.fasta [\\fIOPTIONS\\fP] ");
	addDescription(parser, "Perform Alignment-free k-tuple frequency comparisons from two fasta files.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

        getOptionValue(options.klen, parser, "klen");
        getOptionValue(options.nohits, parser, "num-hits");
        getOptionValue(options.markovOrder, parser, "markov-order");
        getOptionValue(options.type, parser, "distance-type");
        options.noreverse = isSet(parser, "no-reverse");
	getOptionValue(options.queryFileName, parser, "query-file");
	getOptionValue(options.referenceFileName, parser, "reference-file");

	return seqan::ArgumentParser::PARSE_OK;

}

//this function does the work
void workfunction(CharString input1, CharString input2, ModifyStringOptions & options)
{

	//make definitions
	String<Dna5String> kmers = defineKmers(options.klen);
	bool firsttime = true;
	StringSet<CharString> queryids;
        StringSet<Dna5String> queryseqs;
        StringSet<CharString> refids;
        StringSet<Dna5String> refseqs;

	SeqFileIn queryFileIn(toCString(input1));
	SeqFileIn refFileIn(toCString(input2));

	//reads query file into RAM
	readRecords(queryids, queryseqs, queryFileIn);

	//reads reference into RAM
	readRecords(refids, refseqs, refFileIn);

        //for each sequence in the query
        for(int q = 0; q < length(queryids); q++)
        {
		CharString queryid = queryids[q];
		Dna5String queryseq = queryseqs[q];

		if(options.noreverse != true)
		{
			queryseq = doRevCompl(queryseq);
		}

		//define counts array to be the same size as the kmers
		int querycounts [length(kmers)];
		double qryMarkovProbs [length(kmers)];

		//return a kmer occurance histogram for the query sequence.
		count(kmers, queryseq, options.klen, querycounts);

		cout << "Query: ";
		for (int i = 0; i < length(kmers); i++){
			cout << querycounts[i] << " " << endl;
		}
		cout << endl;
		

		//if markov, do markov probability of the query sequence
		if(options.type == "d2s" || options.type == "d2star")
		{
			markov(kmers, queryseq, options.klen, querycounts, options.markovOrder, qryMarkovProbs);
		}      		

		typedef array_type::index index;
		array_type allrefcounts(boost::extents[length(refids)][length(kmers)]);

		//start querying reference
                for(int r = 0; r < length(refids); r++)
		{

			int refcounts [length(kmers)];
			double refProbs [length(kmers)];

			//if we do reverse compliment, so it here;
			Dna5String refseq = refseqs[r];
			if(options.noreverse != true && firsttime == true)
			{
				refseq = doRevCompl(refseqs[r]);
			}

			//now count 
			if(options.type == "d2")
			{
                                if(firsttime){
					count(kmers, refseq, options.klen, refcounts);
                                        for(int m = 0; m < length(kmers); m++)
                                        {
                                                allrefcounts[r][m] = refcounts[m];
                                        }
                                }

				cout << "Reference: ";
				for (int i = 0; i < length(kmers); i++){
					cout << refcounts[i] << " " << endl;
				}
				cout << endl;
                        }

		}

		firsttime = false;
	}
}

array_type calculateReference(String<Dna5String> kmers, ModifyStringOptions & options)
{
        StringSet<CharString> refids;
        StringSet<Dna5String> refseqs;
        SeqFileIn refFileIn(toCString(options.referenceFileName));
        array_type allrefcounts(boost::extents[length(refids)][length(kmers)]);

        //reads reference into RAM
        readRecords(refids, refseqs, refFileIn);
	
	//start querying reference
	for(int r = 0; r < length(refids); r++)
	{
		cout << r << endl;
		int refcounts [length(kmers)];
		double refProbs [length(kmers)];

		//if we do reverse compliment, so it here;
		Dna5String refseq = refseqs[r];
		if(options.noreverse != true)
		{
			refseq = doRevCompl(refseqs[r]);
		}

		//now count 
		if(options.type == "d2")
		{
			count(kmers, refseq, options.klen, refcounts);
			for(int m = 0; m < length(kmers); m++)
			{
				allrefcounts[r][m] = refcounts[m];
			}
		}
	}

	cout << "Finished calculating reference" << endl;

	return allrefcounts;
}

////i think the way this needs to work is as follows;
//the function actually reads the next line since the file is the producer!!!!
void worker(String<Dna5String> kmers, array_type reference_counts)
{
	while(1)
	{
		Dna5String queryseq;
		CharString queryid;
		m.lock();
		if(!atEnd(queryFileIn))
		{
			readRecord(queryid, queryseq, queryFileIn);
		} else 
		{
			m.unlock();
			return;
		}
		m.unlock();
		//workfunction(CharString input1, CharString input2, ModifyStringOptions & options);
		
		cout << "Thread ID " <<  this_thread::get_id() << " " << length(queryseq) << endl;
	}
}

int main(int argc, char const ** argv)
{

	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
/*
	//decide if we're running pairwise or not
        if(options.referenceFileName == "")
        {
		workfunction(options.queryFileName, options.queryFileName, options);
		cout << "Doing pairwise" << endl;
        }
        else
        {
		workfunction(options.queryFileName, options.referenceFileName, options);
		cout << "Standard" << endl;
        }
*/
	cout << "Calculating K-mers." << endl;
	//calculate kmers
	String<Dna5String> kmers = defineKmers(options.klen);

	cout << "Creating reference counts" << endl;

	//precalculate reference!!
	typedef array_type::index index;
	array_type reference_counts = calculateReference(kmers, options);

	cout << "Processing" << endl;

	open(queryFileIn, (toCString(options.queryFileName)));
	num_threads = 4;
	thread workers[num_threads];
	for(int w = 0; w < num_threads; w++)
	{
		workers[w] = thread(worker, kmers, reference_counts);
		cout << "thread " << w << " started." << endl;
	}

	for(int w = 0; w < num_threads; w++)
	{
		workers[w].join();
	}
/*
        //make definitions
        String<Dna5String> kmers = defineKmers(options.klen);

	//read in query sequence file
	SeqFileIn seqFileIn(toCString(options.queryFileName));
	
	StringSet<CharString> qids;
	StringSet<Dna5String> qseqs;
	readRecords(qids, qseqs, seqFileIn);

	bool firsttime = true;

	//read in reference sequence file
	//lets now go through all of the references
	StringSet<CharString> ids;
	StringSet<Dna5String> seqs;

	SeqFileIn seqFileIn2(toCString(options.referenceFileName));
	readRecords(ids, seqs, seqFileIn2);

	//here, if the number of seq's in the reference file is < the number of top
	//hits requested, set top hits to the number of seq's in the reference
	if(length(ids) < options.nohits)
		options.nohits = length(ids);

	typedef array_type::index index;
	typedef array_type2::index index;
	array_type allrefcounts_new(boost::extents[length(ids)][length(kmers)]);
	array_type2 allrefprobs_new(boost::extents[length(ids)][length(kmers)]);

	//define datastructure to store all counts+markovs
	int allrefcounts [0][0];
	double allrefprobs [0][0];

	double matrix [length(qids)][length(ids)];

	cout << "Information " << endl;
	cout << "Query fasta contains " << length(qids) << " seqs " << endl;
	cout << "Reference fasta contains " << length(ids) << " seqs " << endl;
	cout << "With a Kmer size of " << options.klen << " there are " << length(kmers) << " kmers." << endl;

	//for each sequence in the query
	for(int q = 0; q < length(qids); q++)
	{
		clock_t begin = clock();
		CharString id = qids[q];
		Dna5String seq = qseqs[q];

	        if(options.noreverse != true)
	        {
	                seq = doRevCompl(seq);
	        }

		//define counts array to be the same size as the kmers
		int counts [length(kmers)];

		//return the number of occurances of a kmer in query
		count(kmers, seq, options.klen, counts);

		double qryMarkovProbs [length(kmers)];

		//if markov, do markov
		if(options.type == "d2s" || options.type == "d2star")
		{
			markov(kmers, seq, options.klen, counts, options.markovOrder, qryMarkovProbs);
		}	

		//define arrays to record hits and zero them
		double hits [options.nohits];
		int hitpositions [options.nohits];
		for(int i = 0; i < options.nohits; i++)
		{
			hits[i] = 1.0;
			hitpositions[i] = 0;
		}

		//define arrays to store data
		int refcounts [length(kmers)];
		double refMarkovProbs [length(kmers)];

		//start querying reference
		for(int i = 0; i < length(ids); i++){

			//if we do reverse compliment, so it here;
			Dna5String allseq = seqs[i];
		        if(options.noreverse != true && firsttime == true)
		        {
		                allseq = doRevCompl(seqs[i]);
		        }

			if(options.type == "d2")
			{
				if(firsttime){
					count(kmers, allseq, options.klen, refcounts);
					for(int m = 0; m < length(kmers); m++)
					{
						allrefcounts_new[i][m] = refcounts[m];
					}
				}
				double dist = d2(allrefcounts_new, counts, length(kmers), i);
				//recordall(nohits, hits, dist, i, hitpositions);
				matrix[q][i] = dist;
			}
			else if (options.type == "kmer")
			{
                                if(firsttime){
                                        count(kmers, allseq, options.klen, refcounts);
                                        for(int m = 0; m < length(kmers); m++)
                                        {
						allrefcounts_new[i][m] = refcounts[m];
                                        }
                                }
				double dist = euler(allrefcounts_new, counts, length(kmers), i);
				//recordall(nohits, hits, dist, i, hitpositions);
				matrix[q][i] = dist;
			}
			else if (options.type == "d2s")
			{
                                if(firsttime){
                                        count(kmers, allseq, options.klen, refcounts);
					markov(kmers, allseq, options.klen, refcounts, options.markovOrder, refMarkovProbs);
                                        for(int m = 0; m < length(kmers); m++)
                                        {
						allrefcounts_new[i][m] = refcounts[m];
						allrefprobs_new[i][m] = refMarkovProbs[m];
                                        }
                                }
				double dist = d2s(allrefcounts_new, counts, length(kmers), allrefprobs_new, qryMarkovProbs, i);
				//recordall(nohits, hits, dist, i, hitpositions);
				matrix[q][i] = dist;
			}
			else if (options.type == "d2star")
			{
                                if(firsttime){
                                        count(kmers, allseq, options.klen, refcounts);
                                        markov(kmers, allseq, options.klen, refcounts, options.markovOrder, refMarkovProbs);
                                        for(int m = 0; m < length(kmers); m++)
                                        {
						allrefcounts_new[i][m] = refcounts[m];
						allrefprobs_new[i][m] = refMarkovProbs[m];
                                        }
                                }
				double dist = d2star(allrefcounts_new, counts, length(kmers), allrefprobs_new, qryMarkovProbs, i);
				//recordall(nohits, hits, dist, i, hitpositions);	
				matrix[q][i] = dist;
			}
			else 
			{
				std::cout << "It's all gone wrong!" << std::endl;
			}
		}

		firsttime = false;

		std::cout << "Top Hits for " << id << std::endl;
		std::cout << "------------ " << std::endl;

		//print top hits
        	for(int i = 0; i < options.nohits; i++)
		{
                        std::cout << "      " << i << " " << hits[i] << " " << ids[hitpositions[i]] << " " << hitpositions[i] << std::endl;
		}

		clock_t end = clock();
		cout << "Query " << q << " took " << double(end - begin) / CLOCKS_PER_SEC << " secs" << endl;

	}

	ofstream myfile;
	myfile.open ("matrix.txt");
	myfile << "Names ";
	for(int i = 0; i < length(qids); i++)
	{
		myfile << "Contig_" << i << ",";
	}
	myfile << endl;
	for(int i = 0; i < length(qids); i++)
	{

		myfile << ids[i] << ",";

		for(int j = 0; j < length(ids); j++)
		{
			myfile << matrix[i][j] << ",";
		}
		myfile << endl;
	}
	myfile.close();
*/
	return 0;
}
