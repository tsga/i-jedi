SHELL=	/bin/ksh
FC	=$(FCMP)
LIBS	=$(LIBSM)

OBJS	= apply_incr_noahmp_soil.o \
	  soil_increments.o	   \
	  set_soilveg_snippet.o    \
	  sflx_snippet.o

CMD	= apply_soil_incr.x

$(CMD): $(OBJS) 
	${FC} $(FFLAGS) $(OBJS) $(LIBS) -o $(CMD)

apply_incr_noahmp_soil.o: soil_increments.o  apply_incr_noahmp_soil.f90 
	$(FC) $(FFLAGS) -c apply_incr_noahmp_soil.f90

soil_increments.o: set_soilveg_snippet.o sflx_snippet.o soil_increments.f90
	$(FC) $(FFLAGS) -c soil_increments.f90

set_soilveg_snippet.o: set_soilveg_snippet.f90
	$(FC) $(FFLAGS) -c set_soilveg_snippet.f90

sflx_snippet.o: sflx_snippet.f90
	$(FC) $(FFLAGS) -c sflx_snippet.f90


install:
	-cp -p $(CMD) ../../../install/bin/
	-cp -p $(CMD) ../../../build/bin/
clean:
	-rm -f $(OBJS) *.mod #$(CMD)

