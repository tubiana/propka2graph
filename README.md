# propka2graph
This python script take propka output and then do some stat on residues pka.

This programme can be used to get statistic from propka output (by residues).
Create a graph per residues with the pka calculed.

usage: propka2graph.py [-h] -f [FILES [FILES ...]] [-c CHAINS]
                       [-p COMPARAISON] [-b BASIC]

Parse and compare propka files.

optional arguments:
  -h, --help            show this help message and exit
  -f [FILES [FILES ...]], --files [FILES [FILES ...]]
                        propka file(s)
  -c CHAINS, --chains CHAINS
                        chain to select (separated by a coma)
  -p COMPARAISON, --comparaison COMPARAISON
                        Compare 2 propka file (Y/N)
  -b BASIC, --basic BASIC
                        basic graph? (Y/N)

Example : 

python propka2graph.py -f 1ihm.pka -c A -b Y

CAUTION:
It works only with standart amino acids.

