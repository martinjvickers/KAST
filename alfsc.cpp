/*
 *	ALFSC - Alignment-free Sequence Comparison
 *	Version 0.1
 *	Written by Dr. Martin Vickers (mjv08@aber.ac.uk)
 *
 *	Todo:
 *	1) Write unit tests
 *	6) Create nice output files for the results
 *	   	Would be nice if they were compatible 
 *	   	with other software. e.g. blastoutputs
 *	*) Memory management, need to ensure decent use of memory
 *	*) sort out input to allow a single input fasta to work as a pairwise comparison.
 */

#include "common.h"
#include "distances.h"
#include "utils.h"

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array<double, 2> array_type2;

mutex m;
mutex n;
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
	int num_threads;
};

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
	addOption(parser, seqan::ArgParseOption("c", "num-cores", "Number of Cores.", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "num-cores", "1");
	setShortDescription(parser, "Alignment-free sequence comparison.");
	setVersion(parser, "0.0.1");
	setDate(parser, "July 2016");
	addUsageLine(parser, "-q query.fasta -r reference.fasta [\\fIOPTIONS\\fP] ");
	addDescription(parser, "Perform Alignment-free k-tuple frequency comparisons from two fasta files.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// Only extract  options if the program will continue after parseCommandLine()
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

        getOptionValue(options.klen, parser, "klen");
        getOptionValue(options.nohits, parser, "num-hits");
        getOptionValue(options.markovOrder, parser, "markov-order");
        getOptionValue(options.type, parser, "distance-type");
        options.noreverse = isSet(parser, "no-reverse");
	getOptionValue(options.queryFileName, parser, "query-file");
	getOptionValue(options.referenceFileName, parser, "reference-file");
	getOptionValue(options.num_threads, parser, "num-cores");

	return seqan::ArgumentParser::PARSE_OK;

}

/*
 * This is the main body of work. Pop off the next sequence from the 
 * query file and compare it to each sequence in the reference.
 */
void worker(ModifyStringOptions options)
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
	
                if(options.noreverse != true)
                {
                	queryseq = doRevCompl(queryseq);
		}
		
		unordered_map<string, long long int> countmap;
		count(queryseq, options.klen, countmap);

		StringSet<CharString> refids;
		StringSet<Dna5String> refseqs;
		SeqFileIn refFileIn(toCString(options.referenceFileName));

		unordered_map<string,thingy> markovthingy;
	
                //if markov, do markov
                if(options.type == "d2s" || options.type == "d2star")
                {
			markov(queryseq, options.klen, options.markovOrder, markovthingy);
                }       

		//reads reference into RAM
		readRecords(refids, refseqs, refFileIn);

		//to store the top hits from the reference
		double hits [options.nohits];
		int hitpositions [options.nohits];

		for(int i = 0; i < options.nohits; i++)
                {
                        hits[i] = 1.0;
                        hitpositions[i] = 0;
                }

		for(int r = 0; r < length(refids); r++)
		{
			Dna5String referenceseq = refseqs[r];
                        
			if(options.noreverse != true)
                        {
                                referenceseq = doRevCompl(refseqs[r]);
                        }

			double dist;

			if (options.type == "d2")
			{
				unordered_map<string, long long int> refmap;
				count(referenceseq, options.klen, refmap);
				dist = d2(refmap, countmap);
			} 
			else if(options.type == "kmer")
			{
                                unordered_map<string, long long int> refmap;
                                count(referenceseq, options.klen, refmap);
                                dist = euler(refmap, countmap);
			}
			
			else if (options.type == "d2s")
			{
				unordered_map<string,thingy> refmarkovthingy;
				markov(referenceseq, options.klen, options.markovOrder, refmarkovthingy);
				dist = d2s(markovthingy, refmarkovthingy);
			}
			else if (options.type == "d2star")
			{
				unordered_map<string,thingy> refmarkovthingy;
				markov(referenceseq, options.klen, options.markovOrder, refmarkovthingy);
				dist = d2star(markovthingy, refmarkovthingy);
			}
			recordall(options.nohits, hits, dist, r, hitpositions);
		}

		//print out the top hits
		n.lock();
                std::cout << "Top Hits for " << queryid << std::endl;
                std::cout << "------------ " << std::endl;
                for(int i = 0; i < options.nohits; i++)
                {
			cout << "      " << i << " " << hits[i] << " " << refids[hitpositions[i]] << " " << hitpositions[i] << endl;
                }
		n.unlock();
	}
}

int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	
	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	//open file and launch threads
	open(queryFileIn, (toCString(options.queryFileName)));
	thread workers[options.num_threads];
	for(int w = 0; w < options.num_threads; w++)
	{
		workers[w] = thread(worker, options);
	}

	//do not exit until all the threads have finished
	for(int w = 0; w < options.num_threads; w++)
	{
		workers[w].join();
	}
	
	return 0;
}
