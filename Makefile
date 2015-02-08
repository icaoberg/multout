INCLUDE=/usr/include/
MEX=/usr/local/bin/mex

all:
	gcc -o mulcross mulcross.c -L $(INCLUDE) -lm
	gcc -o multout multout.c -L $(INCLUDE) -lm

mex:
	$(MEX) ml_multout.c
