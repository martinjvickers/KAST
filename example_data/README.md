## Using d2-tools

NOTE: Only works on python 2.X

* Prepare files

```
mvickers@n108379:~/Desktop/d2-tools$ cat seq1.fa 
>seq1
AGGCAGCGTACGAACCTACTGGAGTTGCGGTATGGGACCAGGCGACCTCTGATGCAGAGATACAGGAGCGCCGCGCCGGGTCTTCCTTGTAGAAGTCCTG
mvickers@n108379:~/Desktop/d2-tools$ cat seq2.fa 
>seq2
CGGAGACCTCCGTGGACGGGGAAGTCCTGCGCGGGTCAGACGTACGCCCCGATTAGTTGCCCGGACGCCCGGTTGGCAGAAGTGACGGCGACTGCCCTCA
mvickers@n108379:~/Desktop/d2-tools$ cat meh.txt 
seq1.fa
seq2.fa
```

* Do the count

```
./Tuple_Count.py -l meh.txt -k 1
./Tuple_Count.py -l meh.txt -k 2
./Tuple_Count.py -l meh.txt -k 3
./Tuple_Count.py -l meh.txt -k 4
./Tuple_Count.py -l meh.txt -k 5
./Tuple_Count.py -l meh.txt -k 6
./Tuple_Count.py -l meh.txt -k 7
./Tuple_Count.py -l meh.txt -k 8
./Tuple_Count.py -l meh.txt -k 9
./Tuple_Count.py -l meh.txt -k 10
```

* Do the markov

```
./Markov_Probability_ZeroToThree.py -l meh.txt -k 3
```

* Calculate distance

```
./calculate_dissimiliraty.py -k 3 -l meh.txt -o res -d d2
```

## Using ALFSC

```
./ALFSC-contam-v002.py -q example_data/seq1.fa -r example_data/seq2.fa -t d2s -k5 -m 1
```

Get our d2tools results


```
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2 -o d2_k3_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2 -o d2_k5_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2 -o d2_k7_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2 -o d2_k9_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2S -o d2S_k3_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2S -o d2S_k5_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2S -o d2S_k7_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2S -o d2S_k9_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Eu -o Eu_k3_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Eu -o Eu_k5_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Eu -o Eu_k7_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Eu -o Eu_k9_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Ch -o Ch_k3_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Ch -o Ch_k5_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Ch -o Ch_k7_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Ch -o Ch_k9_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2Star -o d2Star_k9_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2Star -o d2Star_k7_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2Star -o d2Star_k5_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2Star -o d2Star_k3_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Ma -o Ma_k9_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Ma -o Ma_k7_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Ma -o Ma_k5_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Ma -o Ma_k3_m0 -r 0
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2 -o d2_k3_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2 -o d2_k5_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2 -o d2_k7_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2 -o d2_k9_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2S -o d2S_k3_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2S -o d2S_k5_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2S -o d2S_k7_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2S -o d2S_k9_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Eu -o Eu_k3_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Eu -o Eu_k5_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Eu -o Eu_k7_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Eu -o Eu_k9_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Ch -o Ch_k3_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Ch -o Ch_k5_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Ch -o Ch_k7_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Ch -o Ch_k9_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2Star -o d2Star_k9_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2Star -o d2Star_k7_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2Star -o d2Star_k5_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2Star -o d2Star_k3_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Ma -o Ma_k9_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Ma -o Ma_k7_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Ma -o Ma_k5_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Ma -o Ma_k3_m1 -r 1
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2 -o d2_k3_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2 -o d2_k5_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2 -o d2_k7_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2 -o d2_k9_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2S -o d2S_k3_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2S -o d2S_k5_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2S -o d2S_k7_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2S -o d2S_k9_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Eu -o Eu_k3_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Eu -o Eu_k5_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Eu -o Eu_k7_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Eu -o Eu_k9_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Ch -o Ch_k3_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Ch -o Ch_k5_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Ch -o Ch_k7_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Ch -o Ch_k9_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2Star -o d2Star_k9_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2Star -o d2Star_k7_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2Star -o d2Star_k5_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2Star -o d2Star_k3_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Ma -o Ma_k9_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Ma -o Ma_k7_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Ma -o Ma_k5_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Ma -o Ma_k3_m3 -r 2
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2 -o d2_k3_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2 -o d2_k5_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2 -o d2_k7_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2 -o d2_k9_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2S -o d2S_k3_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2S -o d2S_k5_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2S -o d2S_k7_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2S -o d2S_k9_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Eu -o Eu_k3_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Eu -o Eu_k5_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Eu -o Eu_k7_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Eu -o Eu_k9_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Ch -o Ch_k3_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Ch -o Ch_k5_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Ch -o Ch_k7_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Ch -o Ch_k9_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d d2Star -o d2Star_k9_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d d2Star -o d2Star_k7_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d d2Star -o d2Star_k5_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d d2Star -o d2Star_k3_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 9 -l meh.txt -d Ma -o Ma_k9_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 7 -l meh.txt -d Ma -o Ma_k7_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 5 -l meh.txt -d Ma -o Ma_k5_m3 -r 3
/usr/bin/python ./calculate_dissimiliraty.py -k 3 -l meh.txt -d Ma -o Ma_k3_m3 -r 3


```
