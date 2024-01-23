#F90 = ifort
#CC = cc
#LD = ifort
F90 = gfortran
CC = gcc
LD = gfortran
FFLAGS = -O -ffree-line-length-none

#FFLAGS = -O -ffree-line-length-500
#FFLAGS =  -O -mcmodel=large
#FFLAGS =  -CB -check all -traceback

CFLAGS = -O4 -DDEBUG

LIBSDIR = 
CRT0DIR =

FLIB = 

MYLIB=

CRT0 =
###
EXECFILE = run.exe
MODULES   = define_hoge.o define_phys.o
OBJS1   = main.o sph.o time_step.o fld_Genergy_Renergy_update.o
OBJS2   =
###
SRCSM	= ${MODULES:.o=.f90}
SRCS1   = ${OBJS1:.o=.f90}
SRCS2   = ${OBJS2:.o=.f}

$(EXECFILE): $(MODULES) $(OBJS1) $(OBJS2) 
	rm -f $(EXECFILE)
	$(LD) $(LIBSDIR) $(FFLAGS) $(CRT0) $(MODULES) $(MODULES2) $(MODULES3) $(OBJS1) $(OBJS2) $(FLIB) $(LIBC) $(MYLIB) $(ERRNOLIB) -o $@

$(MODULES):$(SRCSM)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCSM)
$(OBJS1):$(SRCS1)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCS1)
$(OBJS2):$(SRCS2)
	$(F90) $(FFLAGS) $(LIBSDIR) -c $(SRCS2)

clean:
	rm -f *.o *.mod $(EXECFILE)
###
$(MODULES): $(SRCSM)
$(OBJS1): $(SRCS1) $(SRCSM) define_hoge.f90
$(OBJS2): $(SRCS2) define_hoge.f90















