/*
ALFSC - Alignment-free Sequence Comparison
Version 0.0.1
Written by Dr. Martin Vickers (mjv08@aber.ac.uk)

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
#include "distances.h"
#include "utils.h"

//mutex n;
mutex m;
SeqFileIn queryFileIn;
vector<unordered_map<string,long long int>> reference_counts_vec;
vector<unordered_map<string,markov_dat>> reference_markov_vec;

/*
Overload the SeqFileBuffer_ so that it uses Iupac String. In this way 
the input file is checked against Iupac and any non-A/C/G/T is silently 
converted into a N.
Idea put forward by h-2 in issue https://github.com/seqan/seqan/issues/1196
*/
/*
namespace seqan {
	template <typename TString, typename TSSetSpec, typename TSpec>
	struct SeqFileBuffer_<StringSet<TString, TSSetSpec>, TSpec>
	{
		typedef String<Iupac> Type;
	};
}
*/

/*
Parse our commandline options
*/
seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("alfsc");
	addOption(parser, seqan::ArgParseOption("k", "klen", "Kmer Length.", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "klen", "3");
	addOption(parser, seqan::ArgParseOption("d", "debug", "Debug Messages."));
	addOption(parser, seqan::ArgParseOption("q", "query-file", "Path to the file containing your query sequence data.\n", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "query-file", toCString(concat(getFileExtensions(SeqFileIn()), ' ')));
	addOption(parser, seqan::ArgParseOption("r", "reference-file", "Path to the file containing your reference sequence data.", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "reference-file", toCString(concat(getFileExtensions(SeqFileIn()), ' ')));
	//setRequired(parser, "query-file");
	addOption(parser, seqan::ArgParseOption("p", "pairwise-file", "Path to the file containing your sequence data which you will perform pairwise comparison on.", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "pairwise-file", toCString(concat(getFileExtensions(SeqFileIn()), ' ')));
	addOption(parser, seqan::ArgParseOption("m", "markov-order", "Markov Order", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "markov-order", "1");
	addOption(parser, seqan::ArgParseOption("n", "num-hits", "Number of top hits to return", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "num-hits", "10");
	addOption(parser, seqan::ArgParseOption("t", "distance-type", "The method of calculating the distance between two sequences.", seqan::ArgParseArgument::STRING, "STR"));
	setValidValues(parser, "distance-type", "d2 kmer d2s d2star manhattan chebyshev hao dai");
	setDefaultValue(parser, "distance-type", "d2");
	addOption(parser, seqan::ArgParseOption("nr", "no-reverse", "Do not use reverse compliment."));
	addOption(parser, seqan::ArgParseOption("c", "num-cores", "Number of Cores.", seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(parser, seqan::ArgParseOption("u", "use-ram", "Use RAM to store reference counts once computed. Very fast but will use a lot of RAM if you have a large reference and/or large kmer size."));
	setDefaultValue(parser, "num-cores", "1");
	setShortDescription(parser, "Alignment-free sequence comparison.");
	setVersion(parser, "0.0.2");
	setDate(parser, "September 2016");
	addUsageLine(parser, "-q query.fasta -r reference.fasta [\\fIOPTIONS\\fP] ");
	addUsageLine(parser, "-p mydata.fasta [\\fIOPTIONS\\fP] ");
	addDescription(parser, "Perform Alignment-free k-tuple frequency comparisons from two fasta files.");

	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// Only extract options if the program will continue after parseCommandLine()
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	//begin extracting options
        getOptionValue(options.klen, parser, "klen");
        getOptionValue(options.nohits, parser, "num-hits");
        getOptionValue(options.markovOrder, parser, "markov-order");
        getOptionValue(options.type, parser, "distance-type");
        options.noreverse = isSet(parser, "no-reverse");
	options.debug = isSet(parser, "debug");
	options.useram = isSet(parser, "use-ram");
	getOptionValue(options.queryFileName, parser, "query-file");
	getOptionValue(options.referenceFileName, parser, "reference-file");
	getOptionValue(options.pairwiseFileName, parser, "pairwise-file");
	getOptionValue(options.num_threads, parser, "num-cores");

	if(isSet(parser, "pairwise-file")){
		if(isSet(parser, "reference-file") == true || isSet(parser, "query-file") == true)
		{
			cerr << "If you are performing a pairwise comparison, you do not need to specify a query (-q) and a reference (-r) file. If you are performing a reference/query based search you do not need to specify a pairwise-file (-p). See alfsc -h for details." << endl;
			return seqan::ArgumentParser::PARSE_ERROR;
		}
	}

	if(isSet(parser, "reference-file") == true && isSet(parser, "query-file") == false)
	{
		cerr << "You have specified a reference (-r) file but not a query (-q) file. See alfsc -h for details." << endl;
		return seqan::ArgumentParser::PARSE_ERROR;
	}

	if(isSet(parser, "reference-file") == false && isSet(parser, "query-file") == true)
        {
                cerr << "You have specified a query (-q) file but not a reference (-r) file. See alfsc -h for details." << endl;
                return seqan::ArgumentParser::PARSE_ERROR;
        }

	if(isSet(parser, "reference-file") == false && isSet(parser, "query-file") == false && isSet(parser, "pairwise-file") == false)
	{
		cerr << "You have not specifed any input file. See alfsc -h for details." << endl;
                return seqan::ArgumentParser::PARSE_ERROR;
	}

	return seqan::ArgumentParser::PARSE_OK;

}

void pairwise(ModifyStringOptions options)
{
	//if we do this with markov
	if(options.type == "d2s" || options.type == "d2star" || options.type == "hao" || options.type == "dai")
	{
		//iterate through reference_markov_vec
		for(auto const& p: reference_markov_vec)
        	{
			for(auto const& q: reference_markov_vec)
        		{
				if(options.type == "d2s")
					cout << d2s(p,q) << " ";
				else if(options.type == "d2star")
					cout << d2star(p,q) << " ";
				else if(options.type == "hao")
					cout << hao(p,q) << " ";
				else if(options.type == "dai")
					cout << dAI(p,q) << " ";
				else
					cerr << "Error: distance type not defined!" << endl;
			}

			cout << endl;

		}
	} 
	else if(options.type == "d2" || options.type == "kmer" || options.type == "manhatten" || options.type == "chebyshev")
	{
                //iterate through reference_markov_vec
                for(auto const& p: reference_counts_vec)
                {
                        for(auto const& q: reference_counts_vec)
                        {
				if(options.type == "d2")
					cout << d2(p,q) << " ";
				else if(options.type == "kmer")
					cout << euler(p,q) << " ";
				else if(options.type == "manhattan")
                                        cout << manhattan(p,q) << " ";
				else if(options.type == "chebyshev")
                                        cout << chebyshev(p,q) << " ";
				else
					cerr << "Error: distance type not defined!" << endl;
                        }
			cout << endl;
                }
	}
	else {
		cerr << "Error: How'd that happen?" << endl;
	}
}


/*
 * This is the main body of work. Pop off the next sequence from the 
 * query file and compare it to each sequence in the reference.
 */
void worker(ModifyStringOptions options)
{
	while(1)
	{
		//gets the next queryseq off from the file
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
		
		//stores our query counts information
		unordered_map<string, long long int> query_countmap;

                //stores our query markov information
                unordered_map<string,markov_dat> query_markovmap;
	
                //if markov, do markov
                if(options.type == "d2s" || options.type == "d2star" || options.type == "hao" || options.type == "dai")
                {
			markov(queryseq, options.klen, options.markovOrder, query_markovmap);
			if(options.useram == true)
			{
				gettophits(options, query_markovmap, queryid, reference_markov_vec);
			} else {
				gettophits(options, query_markovmap, queryid);
			}
                } else if(options.type == "d2" || options.type == "kmer" || options.type == "manhattan" || options.type == "chebyshev")
		{
			count(queryseq, options.klen, query_countmap);
			if(options.useram == true)
			{
				gettophits(options, query_countmap, queryid, reference_counts_vec);
			} else 
			{
				gettophits(options, query_countmap, queryid);
			}
		} else 
		{
			cout << "WARNING: distance measure specified is not implemented" << endl;
		}
	}
}

void precompute(ModifyStringOptions options, CharString reference)
{
	//begin to read in the file
	StringSet<CharString> refids;
	StringSet<Dna5String> refseqs;
	SeqFileIn refFileIn(toCString(reference));
	readRecords(refids, refseqs, refFileIn);

	for(int r = 0; r < length(refids); r++)
	{
		Dna5String referenceseq = refseqs[r];
		if(options.noreverse != true)
                {
                	referenceseq = doRevCompl(refseqs[r]);
                }

		//if markov, do markov
                if(options.type == "d2s" || options.type == "d2star" || options.type == "hao" || options.type == "dai")
                {
			unordered_map<string, markov_dat> refmap; //store our current count
			markov(referenceseq, options.klen, options.markovOrder, refmap); //count!
			reference_markov_vec.push_back(refmap); //insert results to global store

                } else if(options.type == "d2" || options.type == "kmer" || options.type == "manhattan" || options.type == "chebyshev")
                {
			unordered_map<string, long long int> refmap; //store our current count
			count(referenceseq, options.klen, refmap); //count!
			reference_counts_vec.push_back(refmap); //insert results to global store
                } else {
			cout << "WARNING: distance measure specified is not implemented" << endl;
		}
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

	//if we are doing a pairwise comparison
	if(options.pairwiseFileName != NULL)
	{
		cout << "Performing a pairwise comparison." << endl;
		precompute(options, options.pairwiseFileName);
		//can't really multithread pairwise at this point (maybe we can do this in the future)
		pairwise(options);
	} 
	else if (options.referenceFileName != NULL && options.queryFileName != NULL)
	{
		cout << "Performing a query/reference based search" << endl;

		//precompute the reference
		if(options.useram == true)
		{
			precompute(options, options.referenceFileName);
		}
cout << "precomputed" << endl;
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
	}

	return 0;
}
