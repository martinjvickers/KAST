/*
 *	ALFSC - Alignment-free Sequence Comparion
 *	Version 0.1
 *	Written by Dr. Martin Vickers (mjv08@aber.ac.uk)
 *
 *	Todo:
 *	1) Ensure d2s and d2star are working as desired
 *	2) Create matrices outputs
 *	3) Create taxonomy outputs
 *	4) Improve performance of kmer/markovKmer lookup
 *		e.g. use hash tables
 *	5) Store the query counts and markovProbs from 
 *	   the first run into memory or precalcule and 
 *	   store as reference index file in order to 
 *	   improve performance.
 *	6) Multithread the program, probably best bet is
 *	   to simply run each query input as a thread.
 *	7) Create nice output files for the results
 *	   	Would be nice if they were compatible 
 *	   	with other software. e.g. blastoutputs
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
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

using namespace seqan;
using namespace std;

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

void count(String<Dna5String> kmers, Dna5String seq, int klen, int counts[])
{
	//zero the array
	for(int i = 0; i < length(kmers); i++)
	{
		counts[i] = 0;
	}

	//count the number of occurances of each kmer
	for(int i = 0; i < length(seq)-klen; i++)
	{
		for(int j = 0; j < length(kmers); j++)
		{
			if(infix(seq, i, i+klen) == kmers[j])
			{
				counts[j] = counts[j] + 1;
			}
		}
	}
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

//I want to record the top hits, seq's and names
void recordall(int nohits, double hits[], double value, int seqcurrpos, int hitpos[])
{
        if(value > hits[0]){
                hits[0] = value;
		hitpos[0] = seqcurrpos; //copy position move
                int j;
                double temp;
		int temppos;

                for(int i = 0; i < nohits; i++)
                {
                        j = i;
                        while (j > 0 && hits[j] < hits[j-1])
                        {
                                temp = hits[j];
				temppos = hitpos[j]; //copy position move
                                hits[j] = hits[j-1];
				hitpos[j] = hitpos[j-1]; //copy position move
                                hits[j-1] = temp;
				hitpos[j-1] = temppos; //copy position move
                                j--;
                        }
                }
        }
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

double euler(int ref[], int qry[], int nokmers)
{
        double score = 0.0;
        double rN = 0.0;
        double qN = 0.0;

        for(int i = 0; i < nokmers; i++){
		rN = rN + ref[i];
		qN = qN + qry[i];
        }

	for(int i = 0; i < nokmers; i++){
		double rF = ref[i] / rN;
		double qF = qry[i] / qN;
		score = score + (pow((rF - qF), 2));
	}

        return pow(score, 0.5);
}

double d2s(int refCounts[], int qryCounts[], int nokmers, double refProbs[], double qryProbs[])
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
                rN = rN + refCounts[i];
                rtot = rtot + refProbs[i];
	}

	for(int i = 0; i < nokmers; i++)
	{
		double qC = qryCounts[i];
		double rC = refCounts[i];
		double qP = qryProbs[i];
		double rP = refProbs[i];

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

double d2star(int refCounts[], int qryCounts[], int nokmers, double refProbs[], double qryProbs[])
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
                rN = rN + refCounts[i];
        }

        for(int i = 0; i < nokmers; i++)
        {
                double qC = qryCounts[i];
                double rC = refCounts[i];
                double qP = qryProbs[i];
                double rP = refProbs[i];

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
	return (Dna)'C';
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

int main(int argc, char const ** argv)
{
	// Setup ArgumentParser.
	seqan::ArgumentParser parser("alfsc");
	addOption(parser, seqan::ArgParseOption("k", "klen", "Kmer Length.", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "klen", "3");
	addOption(parser, seqan::ArgParseOption("q", "query-file", "Path to the query file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "query-file");
	addOption(parser, seqan::ArgParseOption("r", "reference-file", "Path to the reference file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "reference-file");
        addOption(parser, seqan::ArgParseOption("m", "markov-order", "Markov Order", seqan::ArgParseArgument::INTEGER, "INT"));
        setDefaultValue(parser, "markov-order", "1");
        addOption(parser, seqan::ArgParseOption("n", "num-hits", "Number of top hits to return", seqan::ArgParseArgument::INTEGER, "INT"));
        setDefaultValue(parser, "num-hits", "10");	
	addOption(parser, seqan::ArgParseOption("t", "distance-type", "The method of calculating the distance between two sequences.", seqan::ArgParseArgument::STRING, "STR"));
	setValidValues(parser, "distance-type", "d2 kmer d2s d2star");
	setDefaultValue(parser, "distance-type", "d2");
	addOption(parser, seqan::ArgParseOption("nr", "no-reverse", "Do not use reverse compliment."));
	setShortDescription(parser, "Alignment-free comparison software.");
	setVersion(parser, "0.1");
	setDate(parser, "September 2015");
	addUsageLine(parser, "-q query.fasta -r reference.fasta [\\fIOPTIONS\\fP] ");
	addDescription(parser, "Perform Alignment-free k-tuple frequency comparisons from two fasta files.");

	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	unsigned klen;
	getOptionValue(klen, parser, "klen");

	int nohits;
	getOptionValue(nohits, parser, "num-hits");

        int markovOrder;
        getOptionValue(markovOrder, parser, "markov-order");

	CharString type;
	getOptionValue(type, parser, "distance-type");

	bool noreverse = false;
	noreverse = isSet(parser, "no-reverse");

	// define the kmers to be used
	String<Dna5String> kmers = defineKmers(klen);

	//read in sequence file
	CharString queryFileName;
	getOptionValue(queryFileName, parser, "query-file");
	SeqFileIn seqFileIn(toCString(queryFileName));
	StringSet<CharString> qids;
	StringSet<Dna5String> qseqs;
	readRecords(qids, qseqs, seqFileIn);

	//!!!!!!!!!!!It's here that we make our loop to process query file.
	for(int q = 0; q < length(qids); q++)
	{

		CharString id = qids[q];
		Dna5String seq = qseqs[q];

	        if(noreverse != true)
	        {
	                seq = doRevCompl(seq);
	        }

		//define counts array to be the same size as the kmers
		int counts [length(kmers)];

		clock_t begin = clock();

		//return the number of occurances of a kmer in query
		count(kmers, seq, klen, counts);

		clock_t end = clock();
		cout << "Time Taken to count the kmers " << double(end - begin) / CLOCKS_PER_SEC << " secs" << endl;

		double qryMarkovProbs [length(kmers)];

		//if markov, do markov
		if(type == "d2s" || type == "d2star")
		{
			clock_t begin = clock();
			markov(kmers, seq, klen, counts, markovOrder, qryMarkovProbs);
			clock_t end = clock();
			cout << "Time Taken to calc Markov the kmers " << double(end - begin) / CLOCKS_PER_SEC << " secs " << endl;
		}	

		//lets now go through all of the references
		StringSet<CharString> ids;
		StringSet<Dna5String> seqs;
	
		//read in reference file
	        CharString referenceFileName;
	        getOptionValue(referenceFileName, parser, "reference-file");
	        SeqFileIn seqFileIn2(toCString(referenceFileName));
		readRecords(ids, seqs, seqFileIn2);

		//define arrays to record and zero them
		double hits [nohits];
		int hitpositions [nohits];

		for(int i = 0; i < nohits; i++)
		{
			hits[i] = 0;
			hitpositions[i] = 0;
		}

		//start querying reference
		for(int i = 0; i < length(ids); i++){
			int refcounts [length(kmers)];
			double refMarkovProbs [length(kmers)];

			//if we do reverse compliment, so it here;
			Dna5String allseq = seqs[i];
		        if(noreverse != true)
		        {
		                allseq = doRevCompl(seqs[i]);
		        }

			if(type == "d2")
			{
				clock_t begin = clock();
				count(kmers, allseq, klen, refcounts);
				recordall(nohits, hits, d2(refcounts, counts, length(kmers)), i, hitpositions);
				clock_t end = clock();
	                        cout << "Time Taken to calc reference values and compare to query " << double(end - begin) / CLOCKS_PER_SEC << " secs " << endl;
			}
			else if (type == "kmer")
			{
				count(kmers, allseq, klen, refcounts);
	                        recordall(nohits, hits, euler(refcounts, counts, length(kmers)), i, hitpositions);
			}
			else if (type == "d2s")
			{
				count(kmers, allseq, klen, refcounts);
				markov(kmers, allseq, klen, refcounts, markovOrder, refMarkovProbs);
				recordall(nohits, hits, d2s(refcounts, counts, length(kmers), refMarkovProbs, qryMarkovProbs), i, hitpositions);
			}
			else if (type == "d2star")
			{
				count(kmers, allseq, klen, refcounts);
        	                markov(kmers, allseq, klen, refcounts, markovOrder, refMarkovProbs);
        	                recordall(nohits, hits, d2star(refcounts, counts, length(kmers), refMarkovProbs, qryMarkovProbs), i, hitpositions);
			}
			else 
			{
				std::cout << "It's all gone wrong!" << std::endl;
			}
		}

		std::cout << "Top Hits for " << id << std::endl;
		std::cout << "------------ " << std::endl;

		//print top hits
        	for(int i = nohits-1; i >= 0; i--)
        	        std::cout << "	" << ids[hitpositions[i]] << " " << hitpositions[i] << " " << hits[i] << std::endl;

	}

	return 0;
}
