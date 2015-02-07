INCLUDE=/usr/include/

all:
	gcc -o mulcross mulcross.c -L $(INCLUDE) -lm
	gcc -o multout multout.c -L $(INCLUDE) -lm

test:
	./multout MULCROSS.DAT MULCROSS.OUT 1000
