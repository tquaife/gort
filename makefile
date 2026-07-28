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


#Dependencies. Without these, editing include/gortt.h rebuilds nothing
#and you end up testing a stale binary - which is especially dangerous
#now that "make test" reports a confident PASS.

${OBJECTS} gortt.o:		include/gortt.h makefile
gortt.o:				include/soil_rho.h

#tav_abs and prospect_DB both "use dataSpec_PDB", so need its .mod first.
#Only the left-to-right order of FOBJECTS was keeping this working, which
#would fail under make -j.
tav_abs.o prospect_DB.o:	dataSpec_PDB.o


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

