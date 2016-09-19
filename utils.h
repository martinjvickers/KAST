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

Iupac getRevCompl(Iupac const & nucleotide);
Dna5String doRevCompl(Dna5String seq);
void count(Dna5String seq, int klen, unordered_map<string, long long int> & map);
void markov(Dna5String seq, int klen, int markovOrder, unordered_map<string,markov_dat> & markovmap);
void recordall(int nohits, double hits[], double value, int seqcurrpos, int hitpos[]);
void gettophits(ModifyStringOptions options, unordered_map<string, long long int> query_countmap, CharString queryid);
void gettophits(ModifyStringOptions options, unordered_map<string, markov_dat> query_markovmap, CharString queryid);

void gettophits(ModifyStringOptions options, unordered_map<string, markov_dat> query_markovmap, CharString queryid, vector<unordered_map<string,markov_dat>> reference_markov_vec);
void gettophits(ModifyStringOptions options, unordered_map<string, long long int> query_countsmap, CharString queryid, vector<unordered_map<string,long long int>> reference_counts_vec);
