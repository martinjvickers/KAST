/*
KAST - 
Version 0.0.7
Written by Dr. Martin Vickers

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

#include "distance.h"
#include "utils.h"

int testeuler(){

	int klen = 3;
	IupacString qryseq = doRevCompl("GATTGCCTCTCATTTTCTCTCCCATATTATAGGGTGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGATCTTCGATTTTTTGGCAACCCAAAATGGAGGCGGACGAACGAGATGATAATGATAAGATGATTCAAAAAGACAATGCACGACAGAGAGAGCAGAAAAGATAA");
	IupacString refseq = doRevCompl("CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC");
	
	ModifyStringOptions options;
	map<string, unsigned int> refcounts = count(refseq, klen);
	map<string, unsigned int> querycounts = count(qryseq, klen);
	
	double dist = euler(options, refcounts, querycounts);

	cout << "Euler " << dist << endl;

	return 0;
}

int testd2(){

        int klen = 3;
        IupacString qryseq = doRevCompl("GATTGCCTCTCATTTTCTCTCCCATATTATAGGGTGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGATCTTCGATTTTTTGGCAACCCAAAATGGAGGCGGACGAACGAGATGATAATGATAAGATGATTCAAAAAGACAATGCACGACAGAGAGAGCAGAAAAGATAA");
        IupacString refseq = doRevCompl("CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC");

        ModifyStringOptions options;
        map<string, unsigned int> refcounts = count(refseq, klen);
        map<string, unsigned int> querycounts = count(qryseq, klen);

        double dist = d2(options, refcounts, querycounts);

        cout << "d2 " << dist << endl;

        return 0;
}

int testmanhattan(){

        int klen = 3;
        IupacString qryseq = doRevCompl("GATTGCCTCTCATTTTCTCTCCCATATTATAGGGTGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGATCTTCGATTTTTTGGCAACCCAAAATGGAGGCGGACGAACGAGATGATAATGATAAGATGATTCAAAAAGACAATGCACGACAGAGAGAGCAGAAAAGATAA");
        IupacString refseq = doRevCompl("CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC");

        ModifyStringOptions options;
        map<string, unsigned int> refcounts = count(refseq, klen);
        map<string, unsigned int> querycounts = count(qryseq, klen);

        double dist = manhattan(options, refcounts, querycounts);

        cout << "Manhattan " << dist << endl;

        return 0;
}


int testCount_1()
{
        IupacString seq = "NTGACTGACTGACTGACTGACTGACTGACTGACN";
        int klen = 3;
        map<string, unsigned int> counts = count(seq, klen);

        for(pair<string, unsigned int> p: counts)
        {
                if(p.first == "TGA" && p.second != 8)
                {
                        cout << "TGA should occur 8 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
                        return 1;
                }
                else if(p.first == "GAC" && p.second != 8)
                {
                        cout << "GAC should occur 8 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
                        return 1;
                }
                else if(p.first == "ACT" && p.second != 7)
                {
                        cout << "ACT should occur 7 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
                        return 1;
                }
                else if(p.first == "CTG" && p.second != 7)
                {
                        cout << "CTG should occur 7 times in str: " << seq << " but it's being counted " << p.second << " times." << endl;
                        return 1;
                }
                else if(p.first != "GAC" && p.first != "TGA" && p.first != "ACT" && p.first != "CTG")
                {
                        cout << "We have an unknown 3mer in this test: " << p.first << endl;
			return 1;
                }
        }

        return 0;
}


/*Simple test of reverse compliment*/
int testRevCompl_1()
{
        Dna5String input_dnaSeq = "TGAC";
	Dna5String expected_output_dnaSeq = "TGACNNNGTCA";
        Dna5String recieved_output_dnaSeq = doRevCompl(input_dnaSeq);
	if(expected_output_dnaSeq == recieved_output_dnaSeq)
	{
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

	return 0 ;

}

int main(int argc, char const ** argv)
{
	int returncode = 0;

	//testing reverse compliments
	if(testRevCompl_1() != 0)
	{
		returncode = 1;
	} 
	else 
	{
		cout << "Test reverse compliment PASSED" << endl;
	}

	if(testRevCompl_2() != 0)
	{
                returncode = 1;
        }
        else
        {
		cout << "Test Single Base getRevCompl PASSED" << endl;
        }

	//test counting
	if(testCount_1() != 0)
	{
                returncode = 1;
        }
        else
        {
		cout << "Test Count PASSED" << endl;
        }

	//test distances
	if(testeuler() != 0)
        {
                returncode = 1;
        }
        else
        {
                cout << "Test Euler PASSED" << endl;
        }

        if(testd2() != 0)
        {
                returncode = 1;
        }
        else
        {
                cout << "Test Euler PASSED" << endl;
        }

        if(testmanhattan() != 0)
        {
                returncode = 1;
        }
        else
        {
                cout << "Test Euler PASSED" << endl;
        }

	return returncode;
}
