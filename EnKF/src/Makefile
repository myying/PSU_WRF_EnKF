FC=mpif90
FCFLAG=-f90=ifort -i4 -w -ftz -align all -fno-alias -fp-model precise -assume realloc_lhs -convert big_endian

all: enkf.mpi replace_perturbationmean_by_initial.exe replace_geo_by_initial.exe replace_xfmean_with_gfs.exe enkf_gfs_hybrid_inflation.exe soairborne_to_3dvar.exe so_to_3dvar.exe ensemble_mean.exe replace_mean.exe replace_mean_outside_site.exe

clean:
	rm -f *.o *.mod *.mpi *.exe

mpi_module.o:
	$(FC) -c $(FCFLAG) -FR mpi_module.f
module_structure.o:
	$(FC) -c $(FCFLAG) -FR module_structure.f
module_netcdf.o:
	$(FC) -c $(FCFLAG) -I${NETCDF}/include module_netcdf.f
module_map_utils.o:
	$(FC) -c $(FCFLAG) module_map_utils.f
module_wrf_tools.o:
	$(FC) -c $(FCFLAG) -FR module_wrf_tools.f
module_radar.o:
	$(FC) -c $(FCFLAG) -FR module_radar.f
cal_roi.o:
	$(FC) -c $(FCFLAG) -FR cal_roi.f
hurricane_center.o:
	$(FC) -c $(FCFLAG) -FR hurricane_center.f
obs_io.o:
	$(FC) -c $(FCFLAG) -FR obs_io.f
xb.o:
	$(FC) -c $(FCFLAG) -FR xb.f
sub_enkf_util.o:
	$(FC) -c $(FCFLAG) -FR sub_enkf_util.f
enkf.o:
	$(FC) -c $(FCFLAG) -FR enkf.f
main.o:
	$(FC) -c $(FCFLAG) -FR main.f
ensemble_mean.o:
	$(FC) -c $(FCFLAG) -FR ensemble_mean.f
replace_mean.o:
	$(FC) -c $(FCFLAG) -FR replace_mean.f

enkf.mpi: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o enkf.o
	$(FC) -o enkf.mpi $(FCFLAG) -FR -I$(NETCDF)/include  module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o enkf.o -L$(NETCDF)/lib -lnetcdf -lnetcdff main.f

replace_perturbationmean_by_initial.exe: module_netcdf.o
	$(FC) -o replace_perturbationmean_by_initial.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o -L$(NETCDF)/lib -lnetcdf -lnetcdff replace_perturbationmean_by_initial.f

replace_geo_by_initial.exe: module_netcdf.o
	$(FC) -o replace_geo_by_initial.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o -L$(NETCDF)/lib -lnetcdf -lnetcdff replace_geo_by_initial.f

replace_mean_outside_site.exe: module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o
	$(FC) -o replace_mean_outside_site.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o -L$(NETCDF)/lib -lnetcdf -lnetcdff replace_mean_outside_site.f

replace_xfmean_with_gfs.exe: module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o
	$(FC) -o replace_xfmean_with_gfs.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o -L$(NETCDF)/lib -lnetcdf -lnetcdff replace_xfmean_with_gfs.f

enkf_gfs_hybrid_inflation.exe: module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o
	$(FC) -o enkf_gfs_hybrid_inflation.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o -L$(NETCDF)/lib -lnetcdf -lnetcdff enkf_gfs_hybrid_inflation.f

soairborne_to_3dvar.exe: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o
	$(FC) -o soairborne_to_3dvar.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o -L$(NETCDF)/lib -lnetcdf -lnetcdff soairborne_to_3dvar.f

so_to_3dvar.exe: module_structure.o mpi_module.o module_netcdf.o module_wrf_tools.o module_map_utils.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o
	$(FC) -o so_to_3dvar.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o -L$(NETCDF)/lib -lnetcdf -lnetcdff so_to_3dvar.f

ensemble_mean.exe: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o
	$(FC) -o ensemble_mean.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o -L$(NETCDF)/lib -lnetcdf -lnetcdff ensemble_mean.f

replace_mean.exe: module_netcdf.o
	$(FC) -o replace_mean.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o -L$(NETCDF)/lib -lnetcdf -lnetcdff replace_mean.f

wrf2obs3dvar: module_netcdf.o module_map_utils.o module_structure.o module_wrf_tools.o
	rm -f wrf2obs3dvar.exe
	$(FC) -o wrf2obs3dvar.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o module_map_utils.o module_structure.o module_wrf_tools.o -L$(NETCDF)/lib -lnetcdf -lnetcdff wrf2obs3dvar.f time_calculation.f
