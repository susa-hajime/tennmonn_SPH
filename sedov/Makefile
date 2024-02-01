#F90 = ifort
#CC = cc
#LD = ifort
F90 = gfortran
CC = gcc
LD = gfortran

FFLAGS = -O -ffree-line-length-512
#FFLAGS =  -O -mcmodel=large
#FFLAGS =  -CB -check all
#FFLAGS =  -O -openmp -mcmodel=large -shared-intel


CFLAGS = -O4 -DDEBUG

LIBSDIR = 
CRT0DIR =

FLIB = 

MYLIB=

CRT0 =
###
EXECFILE = run.exe
MODULES1   = define_hoge.o define_tree.o
MODULES2 =
OBJS1   = main.o sph.o resize.o tree_subs.o maketree.o  nb_tree_walk.o
OBJS2   =

###
SRCSM	= ${MODULES1:.o=.f90}
SRCSM2	= ${MODULES2:.o=.f}
SRCS1   = ${OBJS1:.o=.f90}
SRCS2   = ${OBJS2:.o=.f}

$(EXECFILE): $(MODULES1) $(MODULES2) $(OBJS1) $(OBJS2)
	rm -f $(EXECFILE)
	$(LD) $(LIBSDIR) $(FFLAGS) $(CRT0) $(MODULES1) $(MODULES2) $(OBJS1) $(OBJS2) $(FLIB) $(LIBC) $(MYLIB) $(ERRNOLIB) -o $@

$(MODULES1):$(SRCSM)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCSM)
$(MODULES2):$(SRCSM2)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCSM2)
$(OBJS1):$(SRCS1)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCS1)
$(OBJS2):$(SRCS2)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCS2)


clean:
	rm -f *.o *.mod $(EXECFILE)
###
$(MODULES1): $(SRCSM)
$(MODULES2): $(SRCSM2)
$(OBJS1): $(SRCS1) $(SRCSM) define_tree.f90 define_hoge.f90
$(OBJS2): $(SRCS2) define_tree.f90 define_hoge.f90














