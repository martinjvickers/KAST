Dna getRevCompl(Dna const & nucleotide);
Dna5String doRevCompl(Dna5String seq);
String<Dna5String> defineKmers(int kmerlength);
void count(Dna5String seq, int klen, unordered_map<string, long long int> & map);
void count(String<Dna5String> kmers, Dna5String seq, int klen, int counts[]);
void count(String<Dna5String> kmers, Dna5String seq, int klen, long long int counts[]);
void markov(Dna5String seq, int klen, int markovOrder, unordered_map<string,thingy> & markovthingy);
void markov(String<Dna5String> kmers, Dna5String seq, int klen, long long int counts[], int markovOrder, double kmerProb[]);
void recordall(int nohits, double hits[], double value, int seqcurrpos, int hitpos[]);
