## SYNOPSIS
All programs permit to analyze and visualize the taxonomic abundance of a set of sequences compared against a sequence database with BLAST.

1. `taxoptimizer.py` parses a blast output report generated with option `-m8` and add NCBI Taxonomy database information for each HSP.

2. `rankoptimizer.py` analyzes taxonomy abundance of a set of sequences, firstly pre-processed by `taxoptimizer.py` program, and format result with [Krona](https://github.com/marbl/Krona),
an interactive metagenomic visualization in a Web browser. 

3. `kronaextract.py` extracts sub-list of Query ID from a set of sequences matching a given taxon name and/or their offset number from `rankoptimizer.py` output.

## INSTALL

In order to work properly, some prerequisites are needed:

* **Install [bsddb3](https://pypi.python.org/pypi/bsddb3) python library**
```
$ pip install bsddb3
```
* **Install [Berkeley DB library](http://www.oracle.com)** 
  * *Mac OSX*
```
brew install berkeley-db4
```
  * *Ubuntu/Debian*
```
sudo apt-get install libdb-dev
```
  * *CentOS*
```
sudo yum install libdb-devel
```

* **Install [taxodb_ncbi](https://github.com/C3BI-pasteur-fr/taxodb_ncbi) pack**

See `README.md` to install it

* **Install [taxo_rrna](https://github.com/C3BI-pasteur-fr/taxo_rrna) pack**

See `README.md` to install it

* **Install [golden](https://github.com/C3BI-pasteur-fr/golden) suite**

See `README.md` to install it


* **Install `taxo_pack`**
```
 $ tar xfvz taxo_pack-<version>.tar.gz
 $ cd taxo_pack-<version>
 $ python setup.py install 
 ```
 
## PREPARE DATA:

* Create [golden](https://github.com/C3BI-pasteur-fr/golden) indexes using `goldin`
```
$ cd golden_indexes
$ mkdir uniprot
$ tar xfvz ~/taxo_pack-<version>/test/part_of_uniprot.db.tar.gz
$ mv part_of_uniprot.db uniprot
$ goldin uniprot uniprot
```
[golden](https://github.com/C3BI-pasteur-fr/golden) needs to know where its indexes are located.
To do this set a environment variable which point to directory where are located indexes.
```
$ export GOLDENDATA='/path/to/golden_indexes'
```

* Create NCBI Taxonomy Database using [taxodb_ncbi](https://github.com/C3BI-pasteur-fr/taxodb_ncbi) and/or [taxo_rrna](https://github.com/C3BI-pasteur-fr/taxo_rrna)
```
 $ taxodb_ncbi.py -n names.dmp -d nodes.dmp -b ncbi_taxodb.bdb
```

## USAGE

1. `taxoptimizer.py`
```

```
Command line example:
```
$ taxoptimizer.py -i ~/taxo_pack-2.0/test/sequence_test.blast.m8 -o sequence_test.taxo -b ~/taxo_pack-2.0/test/ncbi_taxodb.bdb -t ncbi
```

2. `rankoptimizer.py`
```
usage: rankoptimizer [options] -i <FILE>

rankoptimizer analyzes the taxonomy abundance of a set of sequences, pre-
processed by taxoptimizer program, and format result with Krona, an
interactive metagenomic visualization in a Web browser. By default, only the
best HSP of each sequence is reported.

optional arguments:
  -h, --help            show this help message and exit

Options:
  -i File, --in File    Tabulated input file. Blast report with additional
                        NCBI Taxonomy database informations from taxoptimizer
                        program. (default: None)
  -s File, --krona_js File
                        Krona javascript library. Official distribution: https
                        ://github.com/marbl/Krona/blob/master/KronaTools/src/k
                        rona-2.0.js (default: None)
  -c TAXCOLUMN, --taxcolumn TAXCOLUMN
                        First column with NCBI Taxonomy informations.
                        (default: 14)
  -C SCORECOLUMN, --scorecolumn SCORECOLUMN
                        Column's number with HSP scores (default: 12)

Output options:
  -k File, --krona File
                        xml output file with Krona Specification. (default:
                        None)
  -t File, --text File  taxonomy abundance with dendrogram tree representation
                        (default: None)
  -v File, --htmlx File
                        html output with Krona specification and Krona
                        javascript library. xml style (default: None)
  -V File, --htmlj File
                        html output with Krona specification and Krona
                        javascript library. json style (default: None)
  -j File, --json File  json output with Krona specification (default: None)
  -p File, --dump File  Python dump output with Krona specification. (default:
                        None)
  -a, --lca             Report lowest common ancestor of the taxonomic
                        abundance (default: False)

Specific options:
  -d DELTA, --delta DELTA
                        Accepted percentage variation of HSP score value.
                        (default: 0)
  -l, --clean           Remove 'cellular organism' from taxonomy. (default:
                        False)
  -R, --rank            Report rank information in xml krona (-k) and html
                        krona (-v) formats. By default, "superkingdom, phylum,
                        ..." information are not analyzed. Note: If two
                        queries have the same taxonomy but only one has rank
                        information, queries are in differents nodes into the
                        tree. (default: False)
  -U, --identical       Reintroduce abundance of identical reads in the
                        analyze. Identical reads are treat as unique sequence.
                        The number of identical reads is put at the end of the
                        query name (separated by '_'). (default: False)

Krona 2.1, an interactive metagenomic visualization tool in a Web browser.
(https://github.com/marbl/Krona/): Ondov BD, Bergman NH, and Phillippy AM.
Interactive metagenomic visualization in a Web browser. BMC Bioinformatics.
2011 Sep 30; 12(1):385.
```
Command line example:
```
$ wget https://github.com/marbl/Krona/blob/master/KronaTools/src/krona-2.0.js
$ rankoptimizer.py -i ~/taxo_pack-2.0/test/sequence_test.tr.bl8.taxo -k R_k.xml -t R_t.txt -v R_v.html -V R_Vj.html -j R_j.json -p R_p.dmp -a -s krona-2.0.js
```
3. `kronaextract.py`
```
usage: kronaextract [options] -i <FILE>  -n <STRING>

optional arguments:
  -h, --help            show this help message and exit

Options:
  -i file, --in file    Xml input file with Krona 2.1 Specification. A
                        rankoptimizer ouptut file is recommended) (default:
                        None)
  -n STRING, --taxo_name STRING
                        Taxonomic name. (default: None)
  -o file, --out file   Output file. (default: None)
  -s str, --split_prefix str
                        Split output file into two files with the given prefix
                        name (default: None)

kronaextract extracts sub-list of Query ID from a set of sequences matching a
given taxon name and/or their offset number in the blast report. The output
file could be split into two files: the 'prefix.seq' file contains reads names
and the 'prefix.offset' file contains the corresponding taxoptimizer's line
offset.
```
Command line example
```
$ kronaextract.py -i R_k.xml -n 'Retroviridae'  -o Retroviridae.out -s Retroviridae
```


