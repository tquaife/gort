CC = gcc
FF = gfortran
CFLAGS = -Wall -g #-O3 
LIBS =  -lm 
INCLS = -I./ -I/usr/include -I./include

vpath %.f90 ./PROSPECT-D

OBJECTS = gortt_pn_kopen.o gortt_brdf.o gortt_albedo.o gortt_lidar.o
FOBJECTS =  dataSpec_PDB.o tav_abs.o prospect_DB.o

		
gortt:	${OBJECTS} gortt.o makefile ${FOBJECTS}
		${FF} ${CFLAGS} -o $@ ${FOBJECTS} gortt.o ${OBJECTS} ${INCLS} ${LIBS}
		#------------gortt compilation complete----------------


.c.o: $<
		$(CC) ${INCLS} $(CFLAGS) -c $<
		
#.f90.o: $<
%.o : %.f90
		$(FF) $(CFLAGS) -c $<
	    

#Regression suite: compares model output against tests/golden.
#See tests/run.sh. "make bless" re-captures the goldens and should
#only be run on a commit that is MEANT to change the numbers.

test:	gortt
		@cd tests && ./run.sh

bless:	gortt
		@cd tests && ./bless.sh


clean:
		\rm -f *.o *~ *% *.mod

#n.b. does not remove the gortt binary
.PHONY:	test bless clean

