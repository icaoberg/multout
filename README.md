multout
=======
A program for finding a robust estimate of shape, location and Mahalanobis distances in high dimension.

This software was used and referenced in the following articles

* Rocke and Woodruff. (1996) [Idenfication of Outliers in Multivariate 
Data](http://www.jstor.org/discover/10.2307/2291724?sid=21105809106893&uid=2&uid=4&uid=3739256&uid=3739864). 
Journal of the American Statistical Association 91:1047-1061.

and

* Markey, Boland and Murphy. (1999) [Towards Objective Selection of Representative 
Microscope Images](http://murphylab.web.cmu.edu/publications/72-markey1999.pdf). Biophysical Journal 
76:2230-2237.
 
### Development branch status
[![Build 
Status](https://travis-ci.org/icaoberg/multout.svg?branch=dev)](https://travis-ci.org/icaoberg/multout)

### Master branch status
[![Build 
Status](https://travis-ci.org/icaoberg/multout.svg?branch=master)](https://travis-ci.org/icaoberg/multout)

Installation
------------
To build multout type

```
make
```

This should make two files: multout and mulcross.

Usage
-----
```
[icaoberg@lanec1 multout]$ ./multout  
multout version 3.03.1
Copyright 1992,93,94,95,96 by David L. Woodruff and David M. Rocke
Usage: multout infile outfile [iterations [parmsfile]]

 where infile contains: p n
                        data record one (p data elements)
                        data record two
                            . . . 
and iterations is an integer controlling the time spent.
For light contamination use 1 iteration; for heavy, try (n)(p) or more.
If it is not specified, (n)(p) is used.

Example: multout suspect.dat suspect.out 10000
```

Example
-------
This repository contains an example input file named MULCROSS.DAT. To use this document as an example
type

```
./multout MULCROSS.DAT MULCROSS.OUT 100
```

You should see output similar to 

```
[icaoberg@lanec1 multout]$ ./multout MULCROSS.DAT MULCROSS.OUT 100
multout version 3.03.1
Copyright 1992,93,94,95,96 by David L. Woodruff and David M. Rocke
Begin Partition Cell 1
Begin Partition Cell 2
Begin Partition Cell 3
Begin Partition Cell 4
Analysis report written to MULCROSS.OUT.
Beginning outlier detection.
Done.
Final report written to MULCROSS.OUT
```

Bugs and Questions
------------------
Feel free to contact the original author, David Woodruff at `dlwoodruff AT ucdavis DOT edu`

To submit bugs about the source code in this repository please visit

https://github.com/icaoberg/multout

For any other inquiries visit those links as well.

TODO
----
- [ ] Port `multout` to python
- [ ] Make test for `multout` for python that mimics the C version
