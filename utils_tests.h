#include "utils.h"

using namespace seqan2;
using namespace std;

SEQAN_DEFINE_TEST(count_dna)
{

}

SEQAN_DEFINE_TEST(mask_count_dna)
{
   String<Dna5> seq = "AGGCAGCGTACGAACCTACTGGAGTTGCGGTATGGGACCAGGCGACCTCTGATGCAGAGATACAGGAGCGCCGCGCCGGGTCTTCCTTGTAGAAGTCCTG";
   String<unsigned> counts;
   int klen = 5;
   int effectiveLength = 1;
   vector<CharString> mask;
   mask.push_back("10000");
   countKmersNew(counts, seq, klen, effectiveLength, mask);

   /*
   Result should be;
   A = 21
   C = 24
   G = 33
   T = 18
   */

   /*Set up a results kmer counts with the known result*/
   String<unsigned> kmerCounts;
   Shape<Dna> myShape;
   resize(myShape, effectiveLength);
   int kmerNumber = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(myShape));
   seqan2::clear(kmerCounts);
   seqan2::resize(kmerCounts, kmerNumber, 0);
   kmerCounts[0] = 21;
   kmerCounts[1] = 24;
   kmerCounts[2] = 33;
   kmerCounts[3] = 18;

   for(int i = 0; i < length(counts); i++)
   {
      SEQAN_ASSERT_EQ(kmerCounts[i], counts[i]);
   }
}

SEQAN_DEFINE_TEST(mask_count_raa)
{
   String<ReducedAminoAcidMurphy10> seq = "MVLTIYPDELVQIVS";
   String<unsigned> counts;
   int klen = 5;
   int effectiveLength = 1;
   vector<CharString> mask;
   mask.push_back("10000");
   countKmersNew(counts, seq, klen, effectiveLength, mask);

   String<unsigned> kmerCounts;
   Shape<ReducedAminoAcidMurphy10> myShape;
   resize(myShape, effectiveLength);
   int kmerNumber = _intPow((unsigned)ValueSize<ReducedAminoAcidMurphy10>::VALUE, weight(myShape));
   seqan2::clear(kmerCounts);
   seqan2::resize(kmerCounts, kmerNumber, 0);

   kmerCounts[0] = 0;
   kmerCounts[1] = 2;
   kmerCounts[2] = 0;
   kmerCounts[3] = 1;
   kmerCounts[4] = 0;
   kmerCounts[5] = 0;
   kmerCounts[6] = 6;
   kmerCounts[7] = 0;
   kmerCounts[8] = 1;
   kmerCounts[9] = 1;

   for(unsigned int i = 0; i < length(counts); i++)
   {
      SEQAN_ASSERT_EQ(kmerCounts[i], counts[i]);
//      String<ReducedAminoAcidMurphy10> orig;
//      unhash(orig, i, effectiveLength);
//      cout << orig << "\t" << counts[i] << endl;
//     cout << "kmerCounts[" << i << "] = " << counts[i] << ";" << endl;
   }
}

SEQAN_DEFINE_TEST(mask_count_aa)
{
   String<AminoAcid> seq = "MVLTIYPDELVQIVS";
   String<unsigned> counts;
   int klen = 5;
   int effectiveLength = 1;
   vector<CharString> mask;
   mask.push_back("10000");
   countKmersNew(counts, seq, klen, effectiveLength, mask);

   String<unsigned> kmerCounts;
   Shape<ReducedAminoAcidMurphy10> myShape;
   resize(myShape, effectiveLength);
   int kmerNumber = _intPow((unsigned)ValueSize<AminoAcid>::VALUE, weight(myShape));
   seqan2::clear(kmerCounts);
   seqan2::resize(kmerCounts, kmerNumber, 0);

   kmerCounts[0] = 0;
   kmerCounts[1] = 0;
   kmerCounts[2] = 0;
   kmerCounts[3] = 1;
   kmerCounts[4] = 1;
   kmerCounts[5] = 0;
   kmerCounts[6] = 0;
   kmerCounts[7] = 0;
   kmerCounts[8] = 1;
   kmerCounts[9] = 0;
   kmerCounts[10] = 0;
   kmerCounts[11] = 2;
   kmerCounts[12] = 1;
   kmerCounts[13] = 0;
   kmerCounts[14] = 0;
   kmerCounts[15] = 1;
   kmerCounts[16] = 0;
   kmerCounts[17] = 0;
   kmerCounts[18] = 0;
   kmerCounts[19] = 1;
   kmerCounts[20] = 0;
   kmerCounts[21] = 2;
   kmerCounts[22] = 0;
   kmerCounts[23] = 1;
   kmerCounts[24] = 0;
   kmerCounts[25] = 0;
   kmerCounts[26] = 0;

   for(int i = 0; i < length(counts); i++)
   {
      SEQAN_ASSERT_EQ(kmerCounts[i], counts[i]);
   }
}
