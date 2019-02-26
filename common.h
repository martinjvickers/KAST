/*
MIT License

Copyright (c) 2019 Martin James Vickers

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

#ifndef COMMON_H
#define COMMON_H

#include <seqan/sequence.h>
#include <seqan/reduced_aminoacid.h>

using namespace seqan;
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
   int num_threads;
   bool debug;
   bool lowram;
   bool phylyp = true;
   bool tabout;
   bool blastlike;
   CharString outputFileName = NULL;
   vector<CharString> mask;
   int effectiveLength;
   int score_cutoff;
};

#endif
