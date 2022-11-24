FLAGS           =
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       =

### This can be defined from outside:
# make COMP=gfortran
UNAME = ${shell uname}
ifeq ($(UNAME),Linux)
	STATIC=yes
	COMP=ifort
	MKL=yes
else
	COMP=gfortran
	STATIC=no
	MKL=no
endif

$(info $(COMP) MKL $(MKL) STATIC $(STATIC) )

MODEL0=harm
MODEL=harm
CMODEL=model_$(MODEL).o
ifneq (,$(findstring NH3,$(MODEL)))
	CMODEL=model_$(MODEL).o nh3code2.o
endif
$(info $(CMODEL) )


DEBUG=no


F77C= $(COMP)
F90C= $(COMP)
F90OPTS=-O2 -I.
F77OPTS=-O2 -I.
ifeq ($(COMP),ifort)
	F90OPTS=-O2 -traceback -I.
	F77OPTS=-O2 -traceback -I.
	MP=-qopenmp
	LAPACK=-mkl 
        LSTATIC=-static -static-intel -static-libgcc -qopenmp-lib=compat -static-libstdc++
else
	F90OPTS=-O2 -fbacktrace -I.
	F77OPTS=-O2 -fbacktrace -I.
	MP=-fopenmp
	LSTATIC=-static -static-libgfortran
	LAPACK=-llapack -lblas
endif
ifeq ($(COMP),gfortran)
ifeq ($(MKL),yes)
	LAPACK=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
endif
endif

LINKOPTS=-I. $(LAPACK) $(MP)
ifeq ($(STATIC),yes)
	LINKOPTS=-I. $(MP) $(LSTATIC) $(LAPACK)
endif

### Debugging gfortran
ifeq ($(COMP),gfortran)
	ifeq ($(DEBUG),yes)
		MP=
		F90C= gfortran
		F77C= gfortran
		F90OPTS=-O0 -g -Wall -fbacktrace -fcheck=all -Wextra -I. -pg
		F77OPTS=-O0 -g -Wall -fbacktrace -fcheck=all -Wextra -I. -pg
		LINKOPTS= -I. -O0 -g -Wall -fbacktrace -Wextra -pg $(LAPACK)
	endif	
endif

### Debugging ifort
ifeq ($(COMP),ifort)
	ifeq ($(DEBUG),yes)
		MP=
		F90C=ifort
		F77C=ifort
		F90OPTS=-O0 -check all -traceback -p -I. -g -pg
		F77OPTS=-O0 -check all -traceback -p -I. -g -pg
		LINKOPTS= -I. -pg $(LAPACK) 
	endif	
endif

BIN=./

OBJS= diag-lapack.o cspline.o
OBJS1=getdata.o gauss.o
OBJS2=diab-proj.o laser.o evaluation_xyz.o splining.o cspline.o rk4.o \
	transform.o diagonal.o DirectDynamics.o propag.o \
	writing.o diag-lapack.o \
	jacobi.o

ifeq ($(MODEL),$(MODEL0))
	EXEMODEL=""
else
	EXEMODEL=_$(MODEL)
endif

all:
	make init
	make init0
	make init1
	make init2
	make init-geoms
	make init_rv
	make init_v
	make spectra
	make dyn
	make internal
	make run_model
	make xyz2

init: $(OBJS) sharc_init.o
	$(F90C) $(F90OPTS) \
	-o sharc_init.exe sharc_init.o $(OBJS) \
	$(LINKOPTS)

init0: $(OBJS) sharc_init0.o sharc_init-aux.o
	$(F90C) $(F90OPTS) \
	-o sharc_init0.exe sharc_init0.o sharc_init-aux.o $(OBJS) \
	$(LINKOPTS)

init1: $(OBJS) sharc_init1.o sharc_init-aux.o
	$(F90C) $(F90OPTS) -o sharc_init1.exe sharc_init1.o sharc_init-aux.o $(OBJS) \
	$(LINKOPTS)

init2: $(OBJS) sharc_init2.o sharc_init-aux.o
	$(F90C) $(F90OPTS) -o sharc_init2.exe sharc_init2.o sharc_init-aux.o $(OBJS) \
	$(LINKOPTS)

init-geoms: $(OBJS) sharc_init-geoms.o sharc_init-aux.o
	$(F90C) $(F90OPTS) -o sharc_init-geoms.exe sharc_init-geoms.o sharc_init-aux.o $(OBJS) \
	$(LINKOPTS)


init_rv:	sharc_init_rv.o
	$(F90C) $(F90OPTS) -o sharc_init_rv.exe sharc_init_rv.o \
        $(LINKOPTS)

init_v:	sharc_init_v.o
	$(F90C) $(F90OPTS) -o sharc_init_v.exe sharc_init_v.o $(OBJS) \
        $(LINKOPTS)

spectra: $(OBJS) ../models/$(CMODEL) $(OBJS1) sharc_spectra.o
	$(F90C) $(F90OPTS) \
	-o sharc_spectra$(EXEMODEL).exe sharc_spectra.o \
        $(OBJS) $(CMODEL) $(OBJS1) \
	$(LINKOPTS)

dyn: $(OBJS2) ../models/$(CMODEL) sharc_dyn.o inp_dyn.o SH.o
	$(F90C) $(F90OPTS) -o sharc$(EXEMODEL).exe \
	inp_dyn.o sharc_dyn.o SH.o $(OBJS2) $(CMODEL) \
	$(LINKOPTS)

internal:	$(OBJS) sharc_internal.o 
	$(F90C) $(F90OPTS) -o sharc_internal.exe sharc_internal.o $(OBJS) \
	$(LINKOPTS)
	
run_model:	$(OBJS) ../models/$(CMODEL) run_model.o
	$(F90C) $(F90OPTS) -o model$(EXEMODEL).exe \
	run_model.o $(CMODEL) $(OBJS) \
	$(LINKOPTS)

xyz2:	
	$(F90C) $(F90OPTS) -o xyz2.exe xyz2.f90 $(LINKOPTS)

clean:
	rm -f *.o *.mod

.SUFFIXES : .f90 .f

.f90.o:
	$(F90C) $(F90OPTS) $(MP) -c $<

.f.o:
	$(F77C) $(F77OPTS) $(MP) -c $<

