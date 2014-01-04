PROG =	sp_dynamic

SRCS =	sp_dynamic.f90 sp_module_advection.f90 sp_module_boundary.f90 \
	sp_module_cldfra.f90 sp_module_constant.f90 sp_module_debug.f90 \
	sp_module_gridvar.f90 sp_module_initiate.f90 sp_module_integrate.f90 \
	sp_module_interpolate.f90 sp_module_model.f90 sp_module_output.f90 \
	sp_module_subgrid.f90 sp_module_tendency.f90 \
	sp_module_timeschemes.f90 sp_module_wsm6.f90

OBJS =	sp_dynamic.o sp_module_advection.o sp_module_boundary.o \
	sp_module_cldfra.o sp_module_constant.o sp_module_debug.o \
	sp_module_gridvar.o sp_module_initiate.o sp_module_integrate.o \
	sp_module_interpolate.o sp_module_model.o sp_module_output.o \
	sp_module_subgrid.o sp_module_tendency.o sp_module_timeschemes.o \
	sp_module_wsm6.o

LIBS =	

F90 = ifort
F90FLAGS = -O3 -w -fp-model precise -openmp
LDFLAGS = -openmp

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

sp_dynamic.o: sp_module_boundary.o sp_module_constant.o sp_module_debug.o \
	sp_module_gridvar.o sp_module_initiate.o sp_module_integrate.o \
	sp_module_model.o sp_module_output.o
sp_module_advection.o: sp_module_constant.o sp_module_debug.o \
	sp_module_gridvar.o sp_module_model.o
sp_module_boundary.o: sp_module_constant.o sp_module_debug.o \
	sp_module_gridvar.o sp_module_model.o
sp_module_cldfra.o: sp_module_constant.o sp_module_debug.o \
	sp_module_gridvar.o sp_module_model.o
sp_module_debug.o: sp_module_constant.o sp_module_model.o
sp_module_gridvar.o: sp_module_constant.o sp_module_debug.o sp_module_model.o
sp_module_initiate.o: sp_module_constant.o sp_module_debug.o \
	sp_module_gridvar.o sp_module_model.o
sp_module_integrate.o: sp_module_boundary.o sp_module_constant.o \
	sp_module_debug.o sp_module_gridvar.o sp_module_model.o \
	sp_module_subgrid.o sp_module_timeschemes.o sp_module_wsm6.o
sp_module_interpolate.o: sp_module_constant.o sp_module_debug.o \
	sp_module_gridvar.o sp_module_model.o
sp_module_model.o: sp_module_constant.o
sp_module_output.o: sp_module_constant.o sp_module_model.o
sp_module_subgrid.o: sp_module_constant.o sp_module_debug.o \
	sp_module_gridvar.o sp_module_model.o
sp_module_tendency.o: sp_module_advection.o sp_module_constant.o \
	sp_module_debug.o sp_module_gridvar.o sp_module_model.o
sp_module_timeschemes.o: sp_module_boundary.o sp_module_constant.o \
	sp_module_debug.o sp_module_gridvar.o sp_module_interpolate.o \
	sp_module_model.o sp_module_tendency.o
sp_module_wsm6.o: sp_module_constant.o sp_module_debug.o sp_module_gridvar.o \
	sp_module_model.o
