#include <seqan/basic.h>
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include "distance.h"
#include "utils.h"

using namespace seqan;
using namespace std;

/*
   Prep the query and reference sequences and perform counts
*/
void prep(String<unsigned> & qrycounts, String<unsigned> & refcounts, unsigned int k)
{
   String<Dna5> qryseq = "AGGCAGCGTACGAACCTACTGGAGTTGCGGTATGGGACCAGGCGACCTCTGATGCAGAGATACAGGAGCGCCGCGCCGGGTCTTCCTTGTAGAAGTCCTG";
   String<Dna5> refseq = "CGGAGACCTCCGTGGACGGGGAAGTCCTGCGCGGGTCAGACGTACGCCCCGATTAGTTGCCCGGACGCCCGGTTGGCAGAAGTGACGGCGACTGCCCTCA";

   String<Dna5> qryseqrc = qryseq;
   reverseComplement(qryseqrc);
   append(qryseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"); // this should probably the same size as options.klen
   append(qryseq, qryseqrc);

   String<Dna5> refseqrc = refseq;
   reverseComplement(refseqrc);
   append(refseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"); // this should probably the same size as options.klen
   append(refseq, refseqrc);

   countKmersNew(qrycounts, qryseq, k);
   countKmersNew(refcounts, refseq, k);
}

void prep(String<unsigned> & qrycounts, String<unsigned> & refcounts, 
          String<double> & qrymarkov, String<double> & refmarkov,
          unsigned int k, unsigned int m)
{
   String<Dna5> qryseq = "AGGCAGCGTACGAACCTACTGGAGTTGCGGTATGGGACCAGGCGACCTCTGATGCAGAGATACAGGAGCGCCGCGCCGGGTCTTCCTTGTAGAAGTCCTG";
   String<Dna5> refseq = "CGGAGACCTCCGTGGACGGGGAAGTCCTGCGCGGGTCAGACGTACGCCCCGATTAGTTGCCCGGACGCCCGGTTGGCAGAAGTGACGGCGACTGCCCTCA";

   String<Dna5> qryseqrc = qryseq;
   reverseComplement(qryseqrc);
   append(qryseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"); // this should probably the same size as options.klen
   append(qryseq, qryseqrc);

   String<Dna5> refseqrc = refseq;
   reverseComplement(refseqrc);
   append(refseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"); // this should probably the same size as options.klen
   append(refseq, refseqrc);

   countKmersNew(qrycounts, qryseq, k);
   countKmersNew(refcounts, refseq, k);

   markov(qrymarkov, qrycounts, qryseq, k, m);
   markov(refmarkov, refcounts, refseq, k, m);
}

void prep_aa(String<unsigned> & qrycounts, String<unsigned> & refcounts,
             unsigned int k)
{
   String<AminoAcid> qryseq = "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHKIPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTLMGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL";
   String<AminoAcid> refseq = "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQIATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH";

   countKmersNew(qrycounts, qryseq, k);
   countKmersNew(refcounts, refseq, k);
}

void prep_aa(String<unsigned> & qrycounts, String<unsigned> & refcounts,
             String<double> & qrymarkov, String<double> & refmarkov,
             unsigned int k, unsigned int m)
{
   String<AminoAcid> qryseq = "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHKIPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTLMGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL";
   String<AminoAcid> refseq = "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQIATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH";

   countKmersNew(qrycounts, qryseq, k);
   countKmersNew(refcounts, refseq, k);

   markov(qrymarkov, qrycounts, qryseq, k, m);
   markov(refmarkov, refcounts, refseq, k, m);
}

/*
   Begin running tests
*/


SEQAN_DEFINE_TEST(d2_dna)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.10619));
   expected_results.push_back(make_pair(5, 0.3537));
   expected_results.push_back(make_pair(7, 0.47872));
   expected_results.push_back(make_pair(9, 0.49457));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(d2(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(d2_aa)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.47198));
   expected_results.push_back(make_pair(5, 0.5));
   expected_results.push_back(make_pair(6, 0.5));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep_aa(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(d2(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(euler_dna)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.10306));
   expected_results.push_back(make_pair(5, 0.09317));
   expected_results.push_back(make_pair(7, 0.10092));
   expected_results.push_back(make_pair(9, 0.10369));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(euler(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(euler_aa)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.11299));
   expected_results.push_back(make_pair(5, 0.113633));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep_aa(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(euler(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(manhattan_dna)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.63265));
   expected_results.push_back(make_pair(5, 1.45833));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(manhattan(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(manhattan_aa)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 1.91489));
   expected_results.push_back(make_pair(5, 2));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep_aa(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(manhattan(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(bc_dna)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.316326530612));
   expected_results.push_back(make_pair(5, 0.729166666667));
   expected_results.push_back(make_pair(7, 0.957446808511));
   expected_results.push_back(make_pair(9, 0.989130434783));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(bray_curtis_distance(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(bc_aa)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.943343));
   expected_results.push_back(make_pair(5, 1));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep_aa(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(bray_curtis_distance(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(ngd_dna)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.316326530612));
   expected_results.push_back(make_pair(5, 0.729166666667));
   expected_results.push_back(make_pair(7, 0.957446808511));
   expected_results.push_back(make_pair(9, 0.989130434783));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(normalised_google_distance(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(ngd_aa)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.957447));
   expected_results.push_back(make_pair(5, 1));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep_aa(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(normalised_google_distance(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(chebyshev_dna)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.03061));
   expected_results.push_back(make_pair(5, 0.01042));
   expected_results.push_back(make_pair(7, 0.00532));
   expected_results.push_back(make_pair(9, 0.00543));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(chebyshev(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(chebyshev_aa)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.0169492));
   expected_results.push_back(make_pair(5, 0.00862069));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      prep_aa(qrycounts, refcounts, result.first);
      SEQAN_ASSERT_IN_DELTA(chebyshev(refcounts, qrycounts), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(d2s_dna)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.43246389442301436));
   expected_results.push_back(make_pair(5, 0.3808187760303444));
   expected_results.push_back(make_pair(7, 0.3111486571852624));
   expected_results.push_back(make_pair(9, 0.3105504265761406));

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      String<double> qrymarkov, refmarkov;
      prep(qrycounts, refcounts, qrymarkov, refmarkov, result.first, 1);
      SEQAN_ASSERT_IN_DELTA(d2s(refcounts, qrycounts, refmarkov, qrymarkov), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(d2s_aa)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.377069));
   expected_results.push_back(make_pair(4, 0.379908));
   expected_results.push_back(make_pair(5, 0.391587));

   unsigned markovOrder = 1;

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      String<double> qrymarkov, refmarkov;
      prep_aa(qrycounts, refcounts, qrymarkov, refmarkov, result.first, markovOrder);
      SEQAN_ASSERT_IN_DELTA(d2s(refcounts, qrycounts, refmarkov, qrymarkov), result.second, 0.0001);
   }
}


SEQAN_DEFINE_TEST(d2star_dna)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.4027100011247771));
   expected_results.push_back(make_pair(5, 0.4333069320392635));
   expected_results.push_back(make_pair(7, 0.4811460538701716));
   expected_results.push_back(make_pair(9, 0.4938153388000316));

   unsigned markovOrder = 1;

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      String<double> qrymarkov, refmarkov;
      prep(qrycounts, refcounts, qrymarkov, refmarkov, result.first, markovOrder);
      SEQAN_ASSERT_IN_DELTA(d2star(refcounts, qrycounts, refmarkov, qrymarkov), result.second, 0.0001);
   }
}

SEQAN_DEFINE_TEST(d2star_aa)
{
   vector<pair<unsigned int, double>> expected_results;
   expected_results.push_back(make_pair(3, 0.497807));
   expected_results.push_back(make_pair(4, 0.500178));
   expected_results.push_back(make_pair(5, 0.500023));

   unsigned markovOrder = 1;

   for(pair<unsigned int, double> result : expected_results)
   {
      String<unsigned> qrycounts, refcounts;
      String<double> qrymarkov, refmarkov;
      prep_aa(qrycounts, refcounts, qrymarkov, refmarkov, result.first, markovOrder);
      SEQAN_ASSERT_IN_DELTA(d2star(refcounts, qrycounts, refmarkov, qrymarkov), result.second, 0.0001);
   }
}

