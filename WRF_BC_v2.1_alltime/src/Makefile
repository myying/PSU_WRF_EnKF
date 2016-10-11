#

include ../compiler_macros

OBJS		= \
		update_wrf_bc.o			\
		module_netcdf_interface.o	\
		module_couple_uv.o

update_wrf_bc.exe:	$(OBJS)
		$(FC) -o $@ $(FCFLAGS) $(OBJS) $(LIBS)

update_wrf_bc.o: \
		module_netcdf_interface.o	\
		module_couple_uv.o		\
		update_wrf_bc.F
		$(FC) -c $(FCFLAGS) $(FPPFLAG) \
		update_wrf_bc.F

module_netcdf_interface.o: \
		module_netcdf_interface.F
		$(FC) -c $(FCFLAGS) $(FPPFLAG) \
		module_netcdf_interface.F

module_couple_uv.o: \
		module_couple_uv.F
		$(FC) -c $(FCFLAGS) $(FPPFLAG) \
		module_couple_uv.F

clean:	
	$(RM) $(RM_LIST)

