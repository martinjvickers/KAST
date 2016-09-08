Dna getRevCompl(Dna const & nucleotide);
Dna5String doRevCompl(Dna5String seq);
String<Dna5String> defineKmers(int kmerlength);
void count(Dna5String seq, int klen, unordered_map<string, long long int> & map);
void markov(Dna5String seq, int klen, int markovOrder, unordered_map<string,markov_dat> & markovmap);
void recordall(int nohits, double hits[], double value, int seqcurrpos, int hitpos[]);
