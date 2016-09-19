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
#include "distances.h"

/*
	HANG ON, there are no N's in the test files. WHY ARE THE RESULTS DIFFERENT?


	This way works, and works really well (better if i could tidy that nasty bit about having to reiterate a bunch of times). 

*/
double d2(unordered_map<string, long long int> refmap, unordered_map<string, long long int> querymap)
{
        double sumqCrC = 0.0;
        double sumqC2 = 0.0;
        double sumrC2 = 0.0;

	//i need to iterate over all unique keys from both refmap and querymap..

	//I'm not sure if i can do this faster. Here I create a map, add all the key to it
	//which should crush any duplicate keys. Since I don't care about the values here
	//it doesn't matter
	unordered_map<string, long long int> ourkmers;
	ourkmers = refmap;
	ourkmers.insert(querymap.begin(),querymap.end());

	//now we know all the kmers we have, iterate over them
	for(pair<string, long long int> p: ourkmers)
        {
                sumqCrC = sumqCrC + (refmap[p.first] * querymap[p.first]);
                sumqC2 = sumqC2 + (querymap[p.first] * querymap[p.first]);
                sumrC2 = sumrC2 + (refmap[p.first] * refmap[p.first]);
        }

        double score = sumqCrC / (sqrt(sumqC2) * sqrt(sumrC2));
        return 0.5*(1-score);
}

/*

*/
double euler(unordered_map<string, long long int> refmap, unordered_map<string, long long int> querymap)
{
        double score = 0.0;
        double rN = 0.0;
        double qN = 0.0;

        //i need to iterate over all unique keys from both refmap and querymap..
	unordered_map<string, long long int> ourkmers;
	ourkmers = refmap;
	ourkmers.insert(querymap.begin(),querymap.end());

	for(pair<string, long long int> p: ourkmers)
	{
		rN = rN + refmap[p.first];
		qN = qN + querymap[p.first];
	}

	for(pair<string, long long int> p: ourkmers)
        {
                double rF = refmap[p.first] / rN;
		double qF = querymap[p.first] / qN;
		score = score + (pow((rF - qF), 2));
        }

        return pow(score, 0.5);
}

/*okay, this does actually work. The errors that I was chasing are actually down to there being lots of
  N's in the yeast reference which when comparing the result with the old method of doing it means you
  have a wildly different number.

  When comparing sm1.fasta/sm2.fasta results you get a very similar score as there are very few N's.
*/
double d2s(unordered_map<string,markov_dat> qrymarkovthingy, unordered_map<string,markov_dat> refmarkovthingy)
{
        double score = 0.0;
        double D2S = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        double rN = 0.0;
        double qN = 0.0;
        double rtot = 0.0;
        double qtot = 0.0;

	for(pair<string, markov_dat> p: qrymarkovthingy)
	{
		qN = qN + p.second.count;
		qtot = qtot + p.second.prob;
	}
        for(pair<string, markov_dat> p: refmarkovthingy)
        {
		rN = rN + p.second.count;
		rtot = rtot + p.second.prob;
        }

        unordered_map<string, markov_dat> ourkmers;
        ourkmers = refmarkovthingy;
        ourkmers.insert(qrymarkovthingy.begin(),qrymarkovthingy.end());

	//this loop is the difficult one because we don't know about kmers that were never in our data
	//but this just shouldn't matter. Why is it different? 
	for(pair<string, markov_dat> p: ourkmers)
	{

                double qC = qrymarkovthingy[p.first].count;
                double rC = refmarkovthingy[p.first].count;
                double qP = qrymarkovthingy[p.first].prob;
                double rP = refmarkovthingy[p.first].prob;

		//i think this is wrong...
                //double qCt = qC - (qN * qP);
                //double rCt = rC - (rN * rP);
		//it should be
		double qCt = qC - (qtot * qP);
                double rCt = rC - (rtot * rP);

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

double d2star(unordered_map<string,markov_dat> qrymarkovthingy, unordered_map<string,markov_dat> refmarkovthingy)
{

        double score = 0.0;
        double D2Star = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        double qN = 0.0;
        double rN = 0.0;
	
        for(pair<string, markov_dat> p: qrymarkovthingy)
        {
                qN = qN + p.second.count;
        }
        for(pair<string, markov_dat> p: refmarkovthingy)
        {
                rN = rN + p.second.count;
        }

	//make a superset of unordered maps so we know all of our kmers
	unordered_map<string, markov_dat> ourkmers;
	ourkmers = refmarkovthingy;
	ourkmers.insert(qrymarkovthingy.begin(),qrymarkovthingy.end());

	//now iterate through it
        for(pair<string, markov_dat> p: ourkmers)
        {

                double qC = qrymarkovthingy[p.first].count;
                double rC = refmarkovthingy[p.first].count;
                double qP = qrymarkovthingy[p.first].prob;
                double rP = refmarkovthingy[p.first].prob;

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
