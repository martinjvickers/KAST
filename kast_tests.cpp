#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include "distance_tests.h"
#include "utils_tests.h"

using namespace seqan;
using namespace std;

SEQAN_BEGIN_TESTSUITE(KAST_tests)
{
   // Call roundUp() tests
   SEQAN_CALL_TEST(d2_dna);
   SEQAN_CALL_TEST(euler_dna);
   SEQAN_CALL_TEST(manhattan_dna);
   SEQAN_CALL_TEST(bc_dna);
   SEQAN_CALL_TEST(ngd_dna);
   SEQAN_CALL_TEST(chebyshev_dna);
   SEQAN_CALL_TEST(d2s_dna);
   SEQAN_CALL_TEST(d2star_dna);

   SEQAN_CALL_TEST(d2_aa);
   SEQAN_CALL_TEST(euler_aa);
   SEQAN_CALL_TEST(manhattan_aa);
   SEQAN_CALL_TEST(bc_aa);
   SEQAN_CALL_TEST(ngd_aa);
   SEQAN_CALL_TEST(chebyshev_aa);
   SEQAN_CALL_TEST(d2s_aa);
   SEQAN_CALL_TEST(d2star_aa);

   SEQAN_CALL_TEST(d2_raa);
   SEQAN_CALL_TEST(euler_raa);
   SEQAN_CALL_TEST(manhattan_raa);
   SEQAN_CALL_TEST(bc_raa);
   SEQAN_CALL_TEST(ngd_raa);
   SEQAN_CALL_TEST(chebyshev_raa);
   SEQAN_CALL_TEST(d2s_raa);
   SEQAN_CALL_TEST(d2star_raa);

   //SEQAN_CALL_TEST(count_dna);

   SEQAN_CALL_TEST(mask_count_dna);
   SEQAN_CALL_TEST(mask_count_raa);
   SEQAN_CALL_TEST(mask_count_aa);

   //SEQAN_CALL_TEST(rounding_test_2);
   //SEQAN_CALL_TEST(rounding_test_3);

   // Call checkSorted() tests
   //SEQAN_CALL_TEST(check_sorted_test_1);
   //SEQAN_CALL_TEST(check_sorted_test_2);
}
SEQAN_END_TESTSUITE


