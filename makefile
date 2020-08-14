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
	    

clean:
		\rm -f *.o *~ *% *.mod

