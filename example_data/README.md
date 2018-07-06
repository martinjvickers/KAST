## Using d2-tools

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
