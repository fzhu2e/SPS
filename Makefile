PROG =	sp_dynamic

SRCS =	sp_dynamic.f90 sp_module_boundary.f90 sp_module_constant.f90 \
	sp_module_debug.f90 sp_module_initiate.f90 sp_module_integrate.f90 \
	sp_module_interpolate.f90 sp_module_math.f90 sp_module_model.f90 \
	sp_module_output.f90 sp_module_tendency.f90 sp_module_timeschemes.f90

OBJS =	sp_dynamic.o sp_module_boundary.o sp_module_constant.o \
	sp_module_debug.o sp_module_initiate.o sp_module_integrate.o \
	sp_module_interpolate.o sp_module_math.o sp_module_model.o \
	sp_module_output.o sp_module_tendency.o sp_module_timeschemes.o

LIBS =	

F90 = ifort
#F90FLAGS = -O0 -g -w -DDEBUG
#F90FLAGS = -O3 -w 
F90FLAGS = -O3 -w -fp-model precise
LDFLAGS =
CONFIG =  
#CONFIG = -openmp

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) $(CONFIG) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) $(CONFIG) -c $<

sp_dynamic.o: sp_module_boundary.o sp_module_constant.o sp_module_debug.o \
	sp_module_initiate.o sp_module_integrate.o sp_module_model.o \
	sp_module_output.o
sp_module_boundary.o: sp_module_constant.o sp_module_model.o
sp_module_debug.o: sp_module_constant.o sp_module_model.o
sp_module_initiate.o: sp_module_constant.o sp_module_debug.o \
	sp_module_interpolate.o sp_module_model.o
sp_module_integrate.o: sp_module_boundary.o sp_module_constant.o \
	sp_module_debug.o sp_module_model.o sp_module_timeschemes.o
sp_module_interpolate.o: sp_module_constant.o sp_module_debug.o \
	sp_module_model.o
sp_module_math.o: sp_module_constant.o
sp_module_model.o: sp_module_constant.o
sp_module_output.o: sp_module_constant.o sp_module_model.o
sp_module_tendency.o: sp_module_constant.o sp_module_debug.o \
	sp_module_interpolate.o sp_module_model.o
sp_module_timeschemes.o: sp_module_boundary.o sp_module_constant.o \
	sp_module_debug.o sp_module_interpolate.o sp_module_math.o \
	sp_module_model.o sp_module_tendency.o
