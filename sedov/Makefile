F90 = ifort
CC = cc
LD = ifort
#F90 = gfortran
#CC = gcc
#LD = gfortran

#FFLAGS =  -O -mcmodel=large
#FFLAGS =  -CB -check all
FFLAGS =  -O -openmp -mcmodel=large -shared-intel

CFLAGS = -O4 -DDEBUG

LIBSDIR = 
CRT0DIR =

FLIB = 

MYLIB=

CRT0 =
###
EXECFILE = run.exe
MODULES   = define_hoge.o define_tree.o
MODULES2 =
MODULES3 = distance.o
OBJS1   = main.o sph.o resize.o tree_subs.o maketree.o  nb_tree_walk.o
OBJS2   =
OBJS3   = 
OBJS4   = 
###
SRCSM	= ${MODULES:.o=.f90}
SRCSM2	= ${MODULES2:.o=.f}
SRCSM3	= ${MODULES3:.o=.f90}
SRCS1   = ${OBJS1:.o=.f90}
SRCS2   = ${OBJS2:.o=.f}
SRCS3   = ${OBJS3:.o=.f}
SRCS4   = ${OBJS4:.o=.f90}

$(EXECFILE): $(MODULES) $(MODULES2) $(MODULES3) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(G6CPPLIB)
	rm -f $(EXECFILE)
	$(LD) $(LIBSDIR) $(FFLAGS) $(CRT0) $(MODULES) $(MODULES2) $(MODULES3) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(G6CPPLIB) $(FLIB) $(LIBC) $(MYLIB) $(ERRNOLIB) -o $@

$(MODULES):$(SRCSM)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCSM)
$(MODULES2):$(SRCSM2)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCSM2)
$(MODULES3):$(SRCSM3)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCSM3)
$(OBJS1):$(SRCS1)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCS1)
$(OBJS2):$(SRCS2)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCS2)
$(OBJS3):$(SRCS3)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCS3)
$(OBJS4):$(SRCS4)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCS4)

clean:
	rm -f *.o *.mod $(EXEFILE)
###
$(MODULES): $(SRCSM)
$(MODULES2): $(SRCSM2)
$(OBJS1): $(SRCS1) $(SRCSM) define_tree.f90
$(OBJS2): $(SRCS2) define_tree.f90
$(OBJS3): $(SRCS3)
$(OBJS4): $(SRCS4)














