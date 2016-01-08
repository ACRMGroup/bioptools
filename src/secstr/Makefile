CC = gcc 
LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include
#COPT = -O3 -ansi -pedantic -Wall -I$(INCDIR) -L$(LIBDIR)
COPT = -g -ansi -pedantic -Wall -I$(INCDIR) -L$(LIBDIR)
LIBS = -lbiop -lgen -lm -lxml2

SOBJS  = main.o secstr.o

secstr : $(SOBJS)
	$(CC) $(COPT) -o $@ $(SOBJS) $(LIBS)

.c.o : 
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f $(SOBJS)

distclean : clean
	\rm -f sstruc

