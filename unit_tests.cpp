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
#include "seq.cpp"
#include "distances.cpp"
#include "utils.cpp"

/*
void testd2starN()
{
        cout << "d2star test"<<endl;
        int klen = 3;
        int markovorder = 1;
        Dna5String qryseq = doRevCompl("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
        Dna5String refseq = doRevCompl("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");

        unordered_map<string, markov_dat> refmap;
        markov(refseq, klen, markovorder, refmap);

        cout << "Comparing seq : " << refseq << endl;
        for(pair<string, markov_dat> p: refmap)
        {
                cout << p.first << " " << p.second.count << " " << p.second.prob << endl;
        }

        unordered_map<string, markov_dat> querymap;
        markov(qryseq, klen, markovorder, querymap);

        cout << "Comparing seq : " << qryseq << endl;
        for(pair<string, markov_dat> p: querymap)
        {
                cout << p.first << " " << p.second.count << " " << p.second.prob << endl;
        }

        double result = d2star(refmap, querymap);

        cout << result << endl;

}


void testd2star()
{
	cout << "d2star test"<<endl;
        int klen = 3;
	int markovorder = 1;
        Dna5String qryseq = doRevCompl("gattGCCTCTCATTTTCTCTCCCATATTATAGGGTGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGATCTTCGATTTTTTGGCAACCCAAAATGGAGGCGGACGAACGAGATGATAATGATAAGATGATTCAAAAAGACAATGCACGACAGAGAGAGCAGAAAAGATAA");
        Dna5String refseq = doRevCompl("CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC");

        unordered_map<string, markov_dat> refmap;
        markov(refseq, klen, markovorder, refmap);

        cout << "Comparing seq : " << refseq << endl;
        for(pair<string, markov_dat> p: refmap)
        {
                cout << p.first << " " << p.second.count << " " << p.second.prob << endl;
        }

        unordered_map<string, markov_dat> querymap;
        markov(qryseq, klen, markovorder, querymap);

        cout << "Comparing seq : " << qryseq << endl;
        for(pair<string, markov_dat> p: querymap)
        {
                cout << p.first << " " << p.second.count << " " << p.second.prob << endl;
        }

        double result = d2star(refmap, querymap);

        cout << result << endl;

}

void testd2s(){

	int klen = 3;
	Dna5String qryseq = doRevCompl("GATTGCCTCTCATTTTCTCTCCCATATTATAGGGTGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGATCTTCGATTTTTTGGCAACCCAAAATGGAGGCGGACGAACGAGATGATAATGATAAGATGATTCAAAAAGACAATGCACGACAGAGAGAGCAGAAAAGATAA");
	Dna5String refseq = doRevCompl("CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC");
	
	unordered_map<string, long long int> refmap;
	count(refseq, klen, refmap);

	cout << "Comparing seq : " << refseq << endl;
	for(pair<string, long long int> p: refmap)
	{
		cout << p.first << " " << p.second << endl;
	}

	unordered_map<string, long long int> querymap;
	count(qryseq, klen, querymap);	

	cout << "Comparing seq : " << qryseq << endl;
	for(pair<string, long long int> p: querymap)
        {
		cout << p.first << " " << p.second << endl;
        }

	double result = d2(refmap, querymap);

	cout << result << endl;

}
*/

/*
Same as testCount_1() but with an addition of an N at the beginning and end of the sequence. The results should be identical as N's are dropped automagically
*/
/*void testCount_2()
{
        Dna5String seq = "NTGACTGACTGACTGACTGACTGACTGACTGACN";
        int klen = 3;
        unordered_map<string, long long int> map;
        count(seq, klen, map);
        for(pair<string, long long int> p: map)
        {
                if(p.first == "TGA" && p.second != 8)
                        cout << "TGA should occur 8 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
                else if(p.first == "GAC" && p.second != 8)
                        cout << "GAC should occur 8 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
                else if(p.first == "ACT" && p.second != 7)
                        cout << "ACT should occur 7 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
                else if(p.first == "CTG" && p.second != 7)
                        cout << "CTG should occur 7 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
                else if(p.first != "GAC" && p.first != "TGA" && p.first != "ACT" && p.first != "CTG")
                        cout << "We have an unknown 3mer in this test: " << p.first << endl;
        }
}
*/


/*
Very simple test of counting to ensure that the correct results are being returned for a known result
*/
/*void testCount_1()
{
	Dna5String seq = "TGACTGACTGACTGACTGACTGACTGACTGACKJ";
	int klen = 3;
	unordered_map<string, long long int> map;
	count(seq, klen, map);
	for(pair<string, long long int> p: map)
	{
		if(p.first == "TGA" && p.second != 8)
			cout << "TGA should occur 8 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
		else if(p.first == "GAC" && p.second != 8)
			cout << "GAC should occur 8 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
		else if(p.first == "ACT" && p.second != 7)
			cout << "ACT should occur 7 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
		else if(p.first == "CTG" && p.second != 7)
			cout << "CTG should occur 7 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
		else if(p.first != "GAC" && p.first != "TGA" && p.first != "ACT" && p.first != "CTG")
			cout << "We have an unknown 3mer in this test: " << p.first << endl;
	}
}*/

/*Simple test of reverse compliment*/
int testRevCompl_1()
{
        Dna5String input_dnaSeq = "TGAC";
	Dna5String expected_output_dnaSeq = "TGACNNNGTCA";
        Dna5String recieved_output_dnaSeq = doRevCompl(input_dnaSeq);
	if(expected_output_dnaSeq == recieved_output_dnaSeq)
	{
		cout << "Test reverse compliment PASSED" << endl;
		return 0;
	} else {
		cout << "Test reverse compliment FAILED: Sent " << input_dnaSeq << " expected " << expected_output_dnaSeq << " recieved " << recieved_output_dnaSeq << endl;
		return 1;
	}
}

//check that when you give it something, it returns what you want
int testRevCompl_2()
{
	if(getRevCompl('A') != 'T')
	{
		cout << "ERROR submitted A but got " << getRevCompl('A') << endl;
		return 1;
	}

	if(getRevCompl('T') != 'A')
        {
                cout << "ERROR submitted T but got " << getRevCompl('T') << endl;
		return 1;
        }

	if(getRevCompl('C') != 'G')
        {
                cout << "ERROR submitted C but got " << getRevCompl('C') << endl;
		return 1;
        }

	if(getRevCompl('G') != 'C')
        {
                cout << "ERROR submitted G but got " << getRevCompl('G') << endl;
		return 1;
        }

	if(getRevCompl('N') != 'N')
        {
                cout << "ERROR submitted N but got " << getRevCompl('N') << endl;
		return 1;
        }
	if(getRevCompl('M') != 'N')
	{
		cout << "ERROR submitted M but got " << getRevCompl('M') << endl;
		return 1;
	}
	if(getRevCompl('X') != 'N')
        {
                cout << "ERROR submitted X but got " << getRevCompl('X') << endl;
		return 1;
        }

	cout << "Test Single Base getRevCompl PASSED" << endl;
	return 0 ;

}

int main(int argc, char const ** argv)
{
	int returncode = 0;

	if(!testRevCompl_1())
		returncode = 1;

	if(!testRevCompl_2())
		returncode = 1;

	//tests of the counting function
//	testCount_1();
//	testCount_2();

//	testd2s();
//	testd2star();

//	testd2starN();
	return returncode;
}
