multout
=======
A program for finding a robust estimate of shape, location and Mahalanobis distances in high dimension.

This software was used as a tool in the following article

* Markey, M.K., Boland, M.V., and Murphy, R.F. (1999) [Towards Objective Selection of Representative 
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
multout version 3.03
Copyright 1992,93,94,95,96 by David L. Woodruff and David M. Rocke
Usage: multout infile outfile [iterations [parmsfile]]

 where infile contains: p n
                        data record one (p data elements)
                        data record two
                            . . . 
and iterations is an integer controlling the time spent.
For light contamination use 1 iteration; for heavy, try (n)(p) or more.
If it is not specified, (n)(p) is used.
```
Example: multout suspect.dat suspect.out 10000
