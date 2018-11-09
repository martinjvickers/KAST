#!/bin/sh

#d2, euclid, d2s, D2S, d2s-opt, d2star, D2Star, manhattan, chebyshev, hao, dai, bc, ngd

echo "[INFO] Running Search App Tests"

echo "[TEST] Manhattan"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t manhattan -o test_manhattan_k3
diff test_manhattan_k3 results/search/manhattan_k3
rm test_manhattan_k3

echo "[TEST] d2"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t d2 -o test_d2_k3
diff test_d2_k3 results/search/d2_k3
rm test_d2_k3

echo "[TEST] D2S"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t D2S -o test_D2S_k3
diff test_D2S_k3 results/search/D2S_k3
rm test_D2S_k3

echo "[TEST] D2Star"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t D2Star -o test_D2Star_k3
diff test_D2Star_k3 results/search/D2Star_k3
rm test_D2Star_k3

echo "[TEST] euclid"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t euclid -o test_euclid_k3
diff test_euclid_k3 results/search/euclid_k3
rm test_euclid_k3

echo "[TEST] chebyshev"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t chebyshev -o test_chebyshev_k3
diff test_chebyshev_k3 results/search/chebyshev_k3
rm test_chebyshev_k3

echo "[TEST] bc"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t bc -o test_bc_k3
diff test_bc_k3 results/search/bc_k3
rm test_bc_k3

echo "[TEST] ngd"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t ngd -o test_ngd_k3
diff test_ngd_k3 results/search/ngd_k3
rm test_ngd_k3

echo "[TEST] dai"
./kast -r example_data/seq1.fa -q example_data/yeast.fasta -k 3 -c 1 -t dai -o test_dai_k3
diff test_dai_k3 results/search/dai_k3
rm test_dai_k3

echo "[INFO] Running Pairwise App Tests"

echo "[TEST] ngd"
./kast -p example_data/yeast.fasta -k 3 -c 1 -t ngd -o test_ngd_k3
diff test_ngd_k3 results/pairwise/ngd_k3
rm test_ngd_k3

echo "[TEST] d2"
./kast -p example_data/yeast.fasta -k 3 -c 1 -t d2 -o test_d2_k3
diff test_d2_k3 results/pairwise/d2_k3
rm test_d2_k3

echo "[TEST] d2s"
./kast -p example_data/yeast.fasta -k 3 -c 1 -t D2S -o test_d2s_k3
diff test_d2s_k3 results/pairwise/d2s_k3
rm test_d2s_k3

echo "[TEST] d2star"
./kast -p example_data/yeast.fasta -k 3 -c 1 -t D2Star -o test_d2star_k3
diff test_d2star_k3 results/pairwise/d2star_k3
rm test_d2star_k3

echo "[TEST] euclid"
./kast -p example_data/yeast.fasta -k 3 -c 1 -t euclid -o test_euclid_k3
diff test_euclid_k3 results/pairwise/euclid_k3
rm test_euclid_k3

echo "[TEST] manhattan"
./kast -p example_data/yeast.fasta -k 3 -c 1 -t manhattan -o test_manhattan_k3
diff test_manhattan_k3 results/pairwise/manhattan_k3
rm test_manhattan_k3

echo "[TEST] dai"
./kast -p example_data/yeast.fasta -k 3 -c 1 -t dai -o test_dai_k3
diff test_dai_k3 results/pairwise/dai_k3
rm test_dai_k3

echo "[TEST] chebyshev"
./kast -p example_data/yeast.fasta -k 3 -c 1 -t chebyshev -o test_chebyshev_k3
diff test_chebyshev_k3 results/pairwise/chebyshev_k3
rm test_chebyshev_k3

echo "[INFO] Running Pairwise App Tests AminoAcid"

echo "[TEST] ngd"
./kast -p example_data/protein.fasta -k 3 -c 1 -t ngd -o test_ngd_k3_aa -s aa
rm test_ngd_k3_aa

echo "[TEST] d2"
./kast -p example_data/protein.fasta -k 3 -c 1 -t d2 -o test_d2_k3_aa -s aa
rm test_d2_k3_aa

echo "[TEST] euclid"
./kast -p example_data/protein.fasta -k 3 -c 1 -t euclid -o test_euclid_k3_aa -s aa
rm test_euclid_k3_aa

echo "[TEST] manhattan"
./kast -p example_data/protein.fasta -k 3 -c 1 -t manhattan -o test_manhattan_k3_aa -s aa
rm test_manhattan_k3_aa
