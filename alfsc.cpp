/*
ALFSC - Alignment-free Sequence Comparison
Version 0.0.5
Written by Dr. Martin Vickers (martin.vickers@jic.ac.uk)

MIT License

Copyright (c) 2017 Martin James Vickers

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
#include "seq.cpp"
#include "utils.cpp"
#include "distances.cpp"
#include <limits>
#include <ctime>
mutex m; //read query mutex
mutex n; //write to console mutex
SeqFileIn queryFileIn;
map<string,bool> kmermap;
vector<Seq> v;

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
	addOption(parser, seqan::ArgParseOption("p", "pairwise-file", "Path to the file containing your sequence data which you will perform pairwise comparison on.", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "pairwise-file", toCString(concat(getFileExtensions(SeqFileIn()), ' ')));
	addOption(parser, seqan::ArgParseOption("m", "markov-order", "Markov Order", seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Output file.", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");
	setDefaultValue(parser, "markov-order", "1");
	addOption(parser, seqan::ArgParseOption("n", "num-hits", "Number of top hits to return", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "num-hits", "10");
	addOption(parser, seqan::ArgParseOption("t", "distance-type", "The method of calculating the distance between two sequences.", seqan::ArgParseArgument::STRING, "STR"));
	setValidValues(parser, "distance-type", "d2 kmer d2s d2s-opt d2star manhattan chebyshev hao dai");
	setDefaultValue(parser, "distance-type", "d2");
	addOption(parser, seqan::ArgParseOption("f", "output-format", ".", seqan::ArgParseArgument::STRING, "STR"));
	setValidValues(parser, "output-format", "tabular");
        setDefaultValue(parser, "output-format", "tabular");
	addOption(parser, seqan::ArgParseOption("nr", "no-reverse", "Do not use reverse compliment."));
	addOption(parser, seqan::ArgParseOption("c", "num-cores", "Number of Cores.", seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(parser, seqan::ArgParseOption("u", "use-ram", "Use RAM to store reference counts once computed. Very fast but will use a lot of RAM if you have a large reference and/or large kmer size."));
	setDefaultValue(parser, "num-cores", "1");
	setShortDescription(parser, "Alignment-free sequence comparison.");
	setVersion(parser, "0.0.5");
	setDate(parser, "January 2017");
	addUsageLine(parser, "-q query.fasta -r reference.fasta -o results.txt [\\fIOPTIONS\\fP] ");
	addUsageLine(parser, "-p mydata.fasta -o results.txt [\\fIOPTIONS\\fP] ");
	addDescription(parser, "Perform Alignment-free k-tuple frequency comparisons from sequences. This can be in the form of two input files (e.g. a reference and a query) or a single file for pairwise comparisons to be made.");

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
	getOptionValue(options.outputFileName, parser, "output-file");
	getOptionValue(options.num_threads, parser, "num-cores");
	getOptionValue(options.output_format, parser, "output-format");

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
		printHelp(parser);
		return seqan::ArgumentParser::PARSE_ERROR;
	}
	if(isSet(parser, "reference-file") == false && isSet(parser, "query-file") == true)
        {
                cerr << "You have specified a query (-q) file but not a reference (-r) file. See alfsc -h for details." << endl;
		printHelp(parser);
                return seqan::ArgumentParser::PARSE_ERROR;
        }
	if(isSet(parser, "reference-file") == false && isSet(parser, "query-file") == false && isSet(parser, "pairwise-file") == false)
	{
		cerr << "You have not specifed any input file. See alfsc -h for details." << endl;
		printHelp(parser);
                return seqan::ArgumentParser::PARSE_ERROR;
	}
	return seqan::ArgumentParser::PARSE_OK;
}

int pairwise(ModifyStringOptions options)
{
	SeqFileIn pairwiseFileIn;
	StringSet<IupacString> pairwiseseq;
        StringSet<CharString> pairwiseid;
	
	if(!open(pairwiseFileIn, (toCString(options.pairwiseFileName))))
	{
		cerr << "Error: could not open file " << toCString(options.pairwiseFileName) << endl;
		return 1;
	}

	readRecords(pairwiseid, pairwiseseq, pairwiseFileIn);	

	const int size = length(pairwiseid);
	double array[size][size];
	int last = 1;

	for(int i = 0; i < length(pairwiseid); i++)
	{
		Seq refseqobj(pairwiseseq[i], pairwiseid[i], options.noreverse, options.klen, options.markovOrder);

		for(int j = 0; j < last; j++)
        	{
			Seq qryseqobj(pairwiseseq[j], pairwiseid[j], options.noreverse, options.klen, options.markovOrder);

			double dist;
			if(options.type == "kmer")
				dist = euler(refseqobj, qryseqobj, options);
			else if(options.type == "d2")
				dist = d2(refseqobj, qryseqobj, options);
			else if(options.type == "d2s")
				dist = d2sopt(refseqobj, qryseqobj, options, kmermap);
			else if(options.type == "d2s-opt")
				dist = d2sopt(refseqobj, qryseqobj, options, kmermap);
			array[i][j] = dist;
			array[j][i] = dist;
		}
		last++;
	}

	for(int i = 0; i < length(pairwiseid); i++)
	{
		for(int j = 0; j < length(pairwiseid); j++)
		{
			cout << array[i][j] << " ";
		}
		cout << endl;
	}

	return 0;
}

//this is where we do stuff
int mainloop(ModifyStringOptions options)
{
	while(1)
	{
		//get the next query record
		IupacString queryseq;
		CharString queryid;

		m.lock();
		if(!atEnd(queryFileIn))
		{
			readRecord(queryid, queryseq, queryFileIn);
		}
		else
		{
			m.unlock();
			return 0;
		}
		m.unlock();

		Seq qryseqobj(queryseq, queryid, options.noreverse, options.klen, options.markovOrder);

		clock_t begin = clock();

		map<double,string> results;
		map<string,double> results2;

		if(options.useram == true)
		{
			for(int i = 0; i < v.size(); i++)
			{
				double dist;
				if(options.type == "kmer")
                                        dist = euler(v[i], qryseqobj, options);
                                else if(options.type == "d2")
                                        dist = d2(v[i], qryseqobj, options);
                                else if(options.type == "d2s")
                                        dist = d2sopt(v[i], qryseqobj, options, kmermap);
                                else if(options.type == "d2s-opt")
                                        dist = d2sopt(v[i], qryseqobj, options, kmermap);

				//record
				results[dist] = toCString(v[i].getID());
				results2[toCString(v[i].getID())] = dist;
				//cout << toCString(v[i].getID()) << endl;
			}
		} else {
			//now go through each reference seq
        	        CharString refid;
                	IupacString refseq;
	                SeqFileIn refFileIn;
			if(!open(refFileIn, (toCString(options.referenceFileName))))
			{
				cerr << "Error: could not open file " << toCString(options.referenceFileName) << endl;
				return 1;
			}

			while(!atEnd(refFileIn))
			{
				readRecord(refid, refseq, refFileIn);
				Seq refseqobj(refseq, refid, options.noreverse, options.klen, options.markovOrder);

				double dist;
	
				if(options.type == "kmer")
					dist = euler(refseqobj, qryseqobj, options);
				else if(options.type == "d2")
					dist = d2(refseqobj, qryseqobj, options);
				else if(options.type == "d2s")
					dist = d2sopt(refseqobj, qryseqobj, options, kmermap);
				else if(options.type == "d2s-opt")
					dist = d2sopt(refseqobj, qryseqobj, options, kmermap);
		
				//record
				results[dist] = toCString(refseqobj.getID());
				results2[toCString(refseqobj.getID())] = dist;
				
			}
		}

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		n.lock();
		cout << "############################ " << queryid << " took " << elapsed_secs << endl;
		int c = 0;
		for(auto const& p: results2)
		{
			cout << p.first << " " << p.second << endl;
			c++;
			if(c > options.nohits)
				break;
		}
		n.unlock();
	}

	return 0;
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
		//can't really multithread pairwise at this point (maybe we can do this in the future)
		pairwise(options);
	}
	else if (options.referenceFileName != NULL && options.queryFileName != NULL)
	{

		//only do this if we want to do it
		if(options.useram == true)
		{
			cout << " reading in reference objects " << endl;

        	        //now go through each reference seq
        	        CharString refid;
        	        IupacString refseq;
        	        SeqFileIn refFileIn;
			if(!open(refFileIn, (toCString(options.referenceFileName))))
			{
				cerr << "Error: could not open file " << toCString(options.referenceFileName) << endl;
				return 1;
			}
			while(!atEnd(refFileIn))
        	        {
				readRecord(refid, refseq, refFileIn);
                        	Seq refseqobj(refseq, refid, options.noreverse, options.klen, options.markovOrder);
				refseqobj.getCounts(options.klen);
				refseqobj.getTotalMarkov(options.klen, options.markovOrder);
				v.push_back(refseqobj);
				cout << sizeof(refseqobj) << endl;
			}
			cout << v.size() << endl;
		}

		if(options.type == "d2s-opt")
		{
			kmermap = makeall(options);
		} else if(options.type == "d2s")
		{
			kmermap = makecomplete(options);
		}

		if(!open(queryFileIn, (toCString(options.queryFileName))))
		{
			cerr << "Error: could not open file " << toCString(options.queryFileName) << endl;
			return 1;
		}

		thread workers[options.num_threads];
		for(int w = 0; w < options.num_threads; w++)
		{
			workers[w] = thread(mainloop, options);
		}

		//do not exit until all the threads have finished
		for(int w = 0; w < options.num_threads; w++)
		{
			workers[w].join();
		}
	}

	return 0;
}
