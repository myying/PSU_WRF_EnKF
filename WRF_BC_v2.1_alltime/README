WRF_BC:

Purpose: update wrf bc file.

source code list (in src directory):
Makefile
module_couple_uv.F
module_netcdf_interface.F
update_wrf_bc.F


Files needed:

wrf_3dvar_output        #3DVAR analysis
wrfbdy_d01              #Original WRF BC
wrfinput_d01		#Original WRF input (only for cycling)

------------------------------------------------------------------
NOTE:

Keep an extra copy of wrfbdy_d01, it will be overwritten.
no wrfinput_d01 copy needed, it won't be changed.

After run "WRF_BC", rename wrf_3dvar_output as wrfinput_d01

------------------------------------------------------------------

Compile:

	make

	(make clean to remove objs and execs)

	(Before compile, make sure you have proper paths to NETCDF.)

Run:

	new_bc.csh

	(make sure you have linked (named) all files mentioned above)

The meaning of namelist file parame.in:
 wrf_3dvar_output_file : WRF3DVAR analysis (output of wrf3dvar).
 wrf_bdy_file          : wrfbdy_d01 (output of SI).
 wrf_input_from_si     : wrfinput_d01 (output of SI).

 cycling = .false.     : false, cold start; true, cycling.
 debug   = .false.     : set true to print out some message.
 low_bdy_only = .false.: if make long-runs, e.g. SST, surface parameters
                         change in the middle of runs, set it to true.

Note:

 wrf_bdy_file will be changed, keep an original SI copy.
 wrf_input_from_si will not be changed.
 wrf_3dvar_output_file will change, only when low_bdy_only is true.
------------------------------------------------------------------

Once again, remember to rename wrf_3dvar_output as wrfinput_d01 before
you run WRF.

