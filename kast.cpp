/*
KAST - Kmer Alignment-free Search Tool
Version 0.0.9
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

#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <string>
#include <thread>
#include <map>
#include <vector>

#include "common.h"

#include "distance.h"
#include "utils.h"
using namespace seqan;
using namespace std;

map<string, bool> kmer_count_map;
SeqFileIn queryFileIn;
mutex m, n, r;
vector<map<string, unsigned int>> ref_counts_vec;
vector<map<string, double>> ref_markov_vec;
StringSet<CharString> referenceids;
StringSet<Dna5String> referenceseqs;
ofstream outfile; //output file
int current_row = 0;
vector< vector<double> > array_threaded;

//main threaded loop
int mainloop(ModifyStringOptions options)
{
	while(1)
	{
		String<AminoAcid> queryseq;
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

		map<double, int> results;

		String<AminoAcid> seq;
		if(options.noreverse == false)
			seq = doRevCompl(queryseq);
		else
			seq = queryseq;

		map<string, unsigned int> querycounts = count(seq, options.klen);
		map<string, double> querymarkov;
		if(options.type == "d2s" || options.type == "hao" || options.type == "d2star")
		{
			querymarkov = markov(options.klen, seq, options.markovOrder, kmer_count_map); // only do this if one of the markov distance methods
		}

		if(options.lowram != true)
                {

			for(int j = 0; j < length(referenceids); j++)
                	{
				double dist;

				if(options.type == "kmer")
					dist = euler(options, ref_counts_vec[j], querycounts);
				else if(options.type == "d2")
					dist = d2(options, ref_counts_vec[j], querycounts);
				else if(options.type == "d2s" || options.type == "d2s-opt")
					dist = d2s(options, kmer_count_map, ref_counts_vec[j], ref_markov_vec[j], querycounts, querymarkov);
				else if(options.type == "d2star")
					dist = d2star(options, kmer_count_map, ref_counts_vec[j], ref_markov_vec[j], querycounts, querymarkov);
				else if(options.type == "hao")
					dist = hao(options, kmer_count_map, ref_counts_vec[j], ref_markov_vec[j], querycounts, querymarkov);
				else if(options.type == "manhattan")
					dist = manhattan(options, ref_counts_vec[j], querycounts);
				else if(options.type == "chebyshev")
					dist = chebyshev(options, ref_counts_vec[j], querycounts);

				results.insert(pair<double, int> (dist, j));
				if(results.size() > options.nohits)
				{
					map<double, int>::iterator it = results.end();
					results.erase(--it);
				}
			}

		} else {
			CharString refid;
			String<AminoAcid> refseq;
			SeqFileIn refFileIn;
			int counter = 0;
			if(!open(refFileIn, (toCString(options.referenceFileName))))
                       	{
				cerr << "Error: could not open file " << toCString(options.referenceFileName) << endl;
				return 1;
			}

			while(!atEnd(refFileIn))
                        {
				readRecord(refid, refseq, refFileIn);

				String<AminoAcid> rseq;
	                        if(options.noreverse == false)
					rseq = doRevCompl(refseq);
				else
					rseq = refseq;

				map<string, unsigned int> refcounts = count(rseq, options.klen);
				map<string, double> refmarkov;
				if(options.type == "d2s" || options.type == "hao" || options.type == "d2star")
				{
					refmarkov = markov(options.klen, rseq, options.markovOrder, kmer_count_map); // only do this if one of the markov distance methods
				}

				double dist;
				if(options.type == "kmer")
					dist = euler(options, refcounts, querycounts);
				else if(options.type == "d2")
					dist = d2(options, refcounts, querycounts);
				else if(options.type == "d2s")
                                        dist = d2s(options, kmer_count_map, refcounts, refmarkov, querycounts, querymarkov);
				else if(options.type == "d2star")
                                        dist = d2star(options, kmer_count_map, refcounts, refmarkov, querycounts, querymarkov);
				else if(options.type == "hao")
                                        dist = hao(options, kmer_count_map, refcounts, refmarkov, querycounts, querymarkov);
				else if(options.type == "manhattan")
                                        dist = manhattan(options, refcounts, querycounts);
				else if(options.type == "chebyshev")
                                        dist = chebyshev(options, refcounts, querycounts);

				results.insert(pair<double, int> (dist, counter));
                                if(results.size() > options.nohits)
                                {
                                        map<double, int>::iterator it = results.end();
                                        results.erase(--it);
                                }

				counter++;
			}
		}

		if(options.tabout == true)
		{
			n.lock();

			StringSet<CharString> split;
			strSplit(split, queryid);
			CharString qName = split[0];
			outfile << "############################ " << length(queryseq) << "\t" << gc_ratio(queryseq) << "\t" << qName << endl;

			for(pair<double, int> p: results)
			{
				StringSet<CharString> split2;
				strSplit(split2, referenceids[p.second]);
				outfile << p.first << "\t" << length(referenceseqs[p.second]) << "\t" << gc_ratio(referenceseqs[p.second]) << "\t" << split2[0] << endl;
//				outfile << p.first << "\t" << length(referenceseqs[p.second]) << "\t" << referenceids[p.second] << endl;
				
			}
			n.unlock();
		} else if(options.blastlike == true) {

			n.lock();
			outfile << ">" << queryid << endl;
			outfile << "Length=" << length(queryseq) << endl;
			outfile << endl;

			for(pair<double, int> p: results)
			{
				outfile << "Query\tLength=" << length(queryseq) << "\tGC-Ratio=" << gc_ratio(queryseq) << endl;
				outfile << "Ref\tLength=" << length(referenceseqs[p.second]) << "\tGC-Ratio=" << gc_ratio(referenceseqs[p.second]) << endl;
			}
			n.unlock();

		} else {
			n.lock();
                	outfile << "############################ " << queryid << endl;
                	for(pair<double, int> p: results)
                	{
                		outfile << referenceids[p.second] << " " << p.first << endl;
                	}
			n.unlock();
		}
	}

	return 0;
}

int pwthread(ModifyStringOptions options, StringSet<CharString> pairwiseid, StringSet<String<AminoAcid>> pairwiseseq)
{
        //not sure about this loop condition
        while(current_row < length(pairwiseid))
        {
                //get row locally and increment it
                int i;
                r.lock();
                i = current_row;
                if(i >= length(pairwiseid))
                        return 0;
                current_row++;
                r.unlock();

		String<AminoAcid> seq;
		if(options.noreverse == false)
			seq = doRevCompl(pairwiseseq[i]);

		map<string, unsigned int> querycounts = count(seq, options.klen);
		map<string, double> querymarkov;
		if(options.type == "d2s" || options.type == "hao" || options.type == "d2star")
		{
			querymarkov = markov(options.klen, seq, options.markovOrder, kmer_count_map);
		}

		for(int j = 0; j < i; j++)
                {
			double dist;
			String<AminoAcid> refseq;
			if(options.noreverse == false)
				refseq = doRevCompl(pairwiseseq[j]);

			map<string, unsigned int> refcounts = count(refseq, options.klen);
			map<string, double> refmarkov;
			if(options.type == "d2s" || options.type == "hao" || options.type == "d2star")
			{
				refmarkov = markov(options.klen, refseq, options.markovOrder, kmer_count_map);
			}
			
			if(options.type == "kmer")
                        	dist = euler(options, refcounts, querycounts);
                        else if(options.type == "d2")
                        	dist = d2(options, refcounts, querycounts);
                        else if(options.type == "d2s")
                        	dist = d2s(options, kmer_count_map, refcounts, refmarkov, querycounts, querymarkov);
                        else if(options.type == "d2star")
                        	dist = d2star(options, kmer_count_map, refcounts, refmarkov, querycounts, querymarkov);
                        else if(options.type == "hao")
                        	dist = hao(options, kmer_count_map, refcounts, refmarkov, querycounts, querymarkov);
                        else if(options.type == "manhattan")
                        	dist = manhattan(options, refcounts, querycounts);
                        else if(options.type == "chebyshev")
                        	dist = chebyshev(options, refcounts, querycounts);
			
			array_threaded[i][j] = dist;
			array_threaded[j][i] = dist;
		}
        }

        return 0;
}

int threaded_pw(ModifyStringOptions options)
{
        //so I'd need a pairwise set of records
        SeqFileIn pairwiseFileIn;
        StringSet<String<AminoAcid>> pairwiseseq;
        StringSet<CharString> pairwiseid;

        if(!open(pairwiseFileIn, (toCString(options.pairwiseFileName))))
        {
                cerr << "Error: could not open file " << toCString(options.pairwiseFileName) << endl;
                return 1;
        }

        readRecords(pairwiseid, pairwiseseq, pairwiseFileIn);

        //sort out our global matrix
        const int size = length(pairwiseid);
        int last = 1;
        array_threaded.resize(size);
        for(int i = 0; i < size; i++)
                array_threaded[i].resize(size);

        //run our threads, this is where we do the work
        //thread workers[options.num_threads];
        int arraySize = options.num_threads;
        thread * workers = new thread[arraySize];
        for(int w = 0; w < options.num_threads; w++)
        {
                workers[w] = thread(pwthread, options, pairwiseid, pairwiseseq);
        }

        //do not exit until all the threads have finished
        for(int w = 0; w < options.num_threads; w++)
        {
                workers[w].join();
        }

	if(options.phylyp = true)
	{
		outfile << length(pairwiseid) << endl;
		for(int i = 0; i < length(pairwiseid); i++)
                {
                        StringSet<CharString> split;
                        strSplit(split, pairwiseid[i]);
                        CharString qName = split[0];

			outfile << qName << "\t";
			for(int j = 0; j < length(pairwiseid); j++)
			{
				outfile << array_threaded[i][j] << "\t";
			}
			outfile << endl;
		}
	} else {
	        //write out pairwise information to file
	        for(int i = 0; i < length(pairwiseid); i++)
	        {
			outfile << pairwiseid[i] << " ";
	                for(int j = 0; j < length(pairwiseid); j++)
	                {
	                        outfile << array_threaded[i][j] << " ";
	                }
	                outfile << endl;
	        }
	}

        return 0;
}

int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	try {
		outfile.open(toCString(options.outputFileName), std::ios_base::out);
	} catch (const ifstream::failure& e) {
		cout << "Error: could not open output file " << toCString(options.outputFileName) << endl;
		return 1; //if you can't open the output file then why bother trying to run the rest
	}

	if(options.pairwiseFileName != NULL)
        {
		//do pairwise
		threaded_pw(options);
	}
	else if (options.referenceFileName != NULL && options.queryFileName != NULL)
        {
		//read in reference
		SeqFileIn seqRefFileIn(toCString(options.referenceFileName));
		readRecords(referenceids, referenceseqs, seqRefFileIn);

		if(options.type == "d2s" || options.type == "hao" || options.type == "d2star")
		{
			kmer_count_map = makecomplete(options);
			cout << "d2s kmers" << length(kmer_count_map) << endl;
		}
		else if(options.type == "d2s-opt")
		{
			kmer_count_map = makequick(options, referenceseqs);
			cout << "d2s-opt kmers" << length(kmer_count_map) << endl;
		}

		//read reference and counts into memory
		for(int i = 0; i < length(referenceids); i++)
		{
			String<AminoAcid> seq;
			if(options.noreverse == false)
				seq = doRevCompl(referenceseqs[i]);
			else
				seq = referenceseqs[i];

			ref_counts_vec.push_back(count(seq, options.klen));
			if(options.type == "d2s" || options.type == "hao" || options.type == "d2star" || options.type == "d2s-opt")
				ref_markov_vec.push_back(markov(options.klen, seq, options.markovOrder, kmer_count_map));
		}

		if(!open(queryFileIn, (toCString(options.queryFileName))))
		{
			cerr << "Error: could not open file " << toCString(options.queryFileName) << endl;
			return 1;
		}

		int arraySize = options.num_threads;
		thread * workers = new thread[arraySize];

		for(int w = 0; w < options.num_threads; w++)
		{
			workers[w] = thread(mainloop, options);
		}

		for(int w = 0; w < options.num_threads; w++)
		{
			workers[w].join();
		}
	}
	return 0;
}
