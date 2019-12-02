##set in environment paths: NETCDF CRTM FFTW

FC=mpif90
FCFLAG=-f90=ifort -i4 -w -ftz -align all -fno-alias -fp-model precise -assume realloc_lhs -convert big_endian -traceback

all: enkf.mpi ensemble_mean.exe replace_mean.exe scale_decompose.exe alignment.exe

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
	$(FC) -c $(FCFLAG) -FR -I$(NETCDF)/include -I$(CRTM_DIR)/include xb.f
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
fft_pack.o:
	$(FC) -c $(FCFLAG) -FR -I$(FFTW)/include fft_pack.f
multiscale_utils.o:
	$(FC) -c $(FCFLAG) -FR multiscale_utils.f

enkf.mpi: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o enkf.o
	$(FC) -o enkf.mpi $(FCFLAG) -FR -I$(NETCDF)/include  module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o enkf.o -L$(CRTM_DIR)/lib -lCRTM -L$(NETCDF)/lib -lnetcdff main.f

ensemble_mean.exe: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o
	$(FC) -o ensemble_mean.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o -L$(CRTM_DIR)/lib -lCRTM -L$(NETCDF)/lib -lnetcdff ensemble_mean.f

scale_decompose.exe: module_structure.o mpi_module.o module_netcdf.o module_wrf_tools.o module_map_utils.o sub_enkf_util.o fft_pack.o multiscale_utils.o
	$(FC) -o scale_decompose.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_wrf_tools.o module_map_utils.o sub_enkf_util.o fft_pack.o multiscale_utils.o -L$(NETCDF)/lib -lnetcdff -L$(FFTW)/lib -lfftw3 scale_decompose.f

alignment.exe: module_structure.o mpi_module.o module_netcdf.o module_wrf_tools.o module_map_utils.o sub_enkf_util.o fft_pack.o multiscale_utils.o
	$(FC) -o alignment.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_wrf_tools.o module_map_utils.o sub_enkf_util.o fft_pack.o multiscale_utils.o -L$(NETCDF)/lib -lnetcdff -L$(FFTW)/lib -lfftw3 alignment.f

replace_mean.exe: module_netcdf.o mpi_module.o
	$(FC) -o replace_mean.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o mpi_module.o -L$(NETCDF)/lib -lnetcdff replace_mean.f

wrf2obs3dvar.exe: module_netcdf.o module_map_utils.o module_structure.o module_wrf_tools.o
	$(FC) -o wrf2obs3dvar.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o module_map_utils.o module_structure.o module_wrf_tools.o -L$(NETCDF)/lib -lnetcdff wrf2obs3dvar.f time_calculation.f

