#ifndef COMMON_H
#define COMMON_H

#include <seqan/sequence.h>
#include <seqan/reduced_aminoacid.h>

using namespace seqan2;
using namespace std;

typedef SimpleType<unsigned char, ReducedAminoAcid_<Murphy10> > ReducedAminoAcidMurphy10;

/*
User defined options struct
*/
struct ModifyStringOptions
{
   int klen;
   int nohits;
   int markovOrder;
   CharString type;
   CharString sequenceType;
   CharString output_format;
   bool noreverse;
   bool calcgc;
   bool noheader;
   CharString queryFileName = NULL;
   CharString referenceFileName = NULL;
   CharString pairwiseFileName = NULL;
   CharString interleavedFileName = NULL;
   int num_threads;
   bool debug;
   bool lowram;
   bool phylyp = true;
   bool tabout;
   bool blastlike;
   CharString outputFileName = NULL;
   vector<CharString> mask;
   int effectiveLength;
   double score_cutoff;
   double filter_percent = 0;
   int64_t filter_bp = 0;
};

#endif
