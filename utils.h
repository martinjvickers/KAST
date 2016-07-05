Dna getRevCompl(Dna const & nucleotide);
Dna5String doRevCompl(Dna5String seq);
String<Dna5String> defineKmers(int kmerlength);
void count(String<Dna5String> kmers, Dna5String seq, int klen, int counts[]);
void count(String<Dna5String> kmers, Dna5String seq, int klen, long long int counts[]);
void markov(String<Dna5String> kmers, Dna5String seq, int klen, long long int counts[], int markovOrder, double kmerProb[]);
void recordall(int nohits, double hits[], double value, int seqcurrpos, int hitpos[]);
