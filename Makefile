CC = g++

# Define any options for the C compiler here
# COPT =  -ansi -Wall -pedantic -g -DGUNZIP_SUPPORT

# Define the archive/library command here
AR = ar r

# If ranlib required, set here; otherwise set to echo
RANLIB = ranlib
#RANLIB = echo

# Destination in which to place the library
LIBDEST = $(HOME)/lib

###########################################################
#          DON'T CHANGE ANYTHING BELOW THIS POINT         #
###########################################################

# Files for libgen.a
OFILES = algorthm.o contdist.o correlat.o datamanp.o descript.o \
         deviate.o discdist.o eigensys.o kernels.o linalg.o mathx.o \
	 matutils.o noncdist.o normdist.o numerror.o rankdist.o \
	 ranking.o residual.o rootfind.o sort.o vecutils.o

all : libnumerics.a

libnumerics.a : $(OFILES)
	$(AR) libnumerics.a $? 
	$(RANLIB) libnumerics.a

.c.o :
	$(CC) $(COPT) -o $@ -c $<

clean :
	rm -f $(OFILES)

install :
	cp libnumerics.a $(LIBDEST)


