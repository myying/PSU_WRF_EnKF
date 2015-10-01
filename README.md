
PSU WRF Ensemble-Variational Data Assimilation System
==============================================================================================================

Codes and Authors:
------------------
Control Scripts:        Michael Ying (yxy159@psu.edu), Jonathan Poterjoy, Yonghui Weng.
                        Last updated 2015/01/20

EnKF component code:    Fuqing Zhang, Yonghui Weng, Michael Ying, Jonathan Poterjoy.
                        Last updated 2014/08/01

4DVar compoent:         WRFDA package from NCAR.

WRF model:              from NCAR

WRFDA_4denvar:          WRFDA code modified for 4DEnVar: Jonathan Poterjoy.

MULTI_INC:              utility programs for running multi-incremental 4DVar

WRF_BC_v2.1_alltime:    WRF_BC_v2.1 modified by Zhiyong Meng to update/perturb lateral boundary condition
                        at a given time.

--------------------------------------------------------------------------------------------------------------

List of components

Top level files:
----------------
 config/katrina_enkf           configuration file(s), see README.config for details

 run_cycle.sh                  runs cycling DA from DATE_START (initialization step) to DATE_CYCLE_END

 run_forecast.sh               runs forecast from analyses at each cycle to DATE_END

 run_gen_be.sh                 runs the model and calculate be.dat from outputs using gen_be package

 gen_be/                       WRFDA gen_be package modified to allow parallelization

Namelists generators:
---------------------
 namelist_wps.sh               creates namelist.wps for WPS input

 namelist_wrf.sh               creates namelist.input for WRF input.
                               usage: namelist_wrf.sh use_for domain_id
                               use_for=wrf    - for WRF runs (all domains)
                               use_for=wrfvar - for da_wrfvar.exe (1 domain)
                               use_for=ndown  - for ndown.exe (2 domains, parent and child)

 namelist_wrfvar.sh            creates namelist.input for da_wrfvar.exe (the &wrfvar parts)

 namelist_obsproc.sh           creates namelist.obsproc for obsproc.exe 

Modules:
--------
 module_icbc.sh                runs geogrid.exe ungrib.exe and metgrid.exe from WPS and real.exe from WRF 
                               to prepare initial and boundary conditions.
                               input:  first guess ($FG_DIR) geog ($GEOG_DIR)
                               output: $WORK_DIR/rc/$DATE/wrfbdy_d01|wrfinput_d0?

 module_obsproc.sh             runs obsproc.exe to create LITTLE_R formatted inputs for EnKF and WRFDA
                               input: raw observations data (NCAR_LITTLE_R, MADIS, BUFR, ...)
                               output: $WORK_DIR/obs/$DATE/obs_gts_<time>.3DVAR|4DVAR
                               Note: if data is preprocessed, this module can be omitted.

 module_perturb_ic.sh          for generating ensemble perturbations.
                               input: be.dat, $WORK_DIR/rc/$DATE/wrfinput_d01
                               output: 100 random perturbations (rc/random_samples) and 
                                       the first NUM_ENS becomes the initial ensemble. 
                               If having nested domains, the perturbations are nested down to create 
                               wrfinput files for the nested domains.

 module_enkf.sh                runs the EnKF component.
                               input: observations, 
                                      $WORK_DIR/fc/$PREVDATE/wrfinput_<domain_id>_<DATE>_<member_id>
                               output: $WORK_DIR/fc/$DATE/wrfinput_<domain_id>_<member_id>
                               If coupled with 4DVar, the ensemble mean is replaced with 4DVar analysis
                               If REPLACE_MEAN=true, the ensemble mean is also replaced with specified
                               file.

 module_4dvar.sh               runs the 4DVar component.
                               input: observations, fg fg02, ep (if coupled with EnKF)
                               output: $WORK_DIR/fc/$DATE/wrfinput_<domain_id>_<DATE+OBS_WIN_MIN>

 module_wrf_ens.sh             runs ensemble forecast to next EnKF analysis time.
                               The boundary condition is perturbed using random_samples.
                               input: $WORK_DIR/fc/$DATE/wrfinput_<domain_id>_<member_id>
                                      $WORK_DIR/fc/wrfbdy_d01_<member_id>
                               output: $WORK_DIR/fc/$DATE/wrfinput_<domain_id>_<NEXTDATE>_<member_id>

 module_wrf_window.sh          runs wrf forecast from 4DVar analysis
                               input: $WORK_DIR/fc/$DATE/wrfinput_<domain_id>_<DATE+OBS_WIN_MIN>
                                      $WORK_DIR/fc/wrfbdy_d01
                               output: $WORK_DIR/fc/$DATE/wrfinput_<domain_id>_<DATE> (coupled with EnKF)
                                    or $WORK_DIR/fc/$DATE/wrfinput_<domain_id>_<NEXTDATE+OBS_WIN_MIN> 
                                    (4DVar only)

 module_wrf_window1.sh         runs wrf forecast across observation window
                               input: $WORK_DIR/fc/$PREVDATE/wrfinput_<domain_id>_<DATE+OBS_WIN_MIN>
                                      $WORK_DIR/fc/wrfbdy_d01_window
                               output: $WORK_DIR/fc/$PREVDATE/wrfinput_<domain_id>_<DATE+OBS_WIN_MAX>

 module_wrf_ens_window1.sh     runs ensemble forecast across observation window (preparing for ep*)
                               input: $WORK_DIR/fc/$PREVDATE/wrfinput_<domain_id>_<DATE+OBS_WIN_MIN>_<member_id>
                                      $WORK_DIR/fc/wrfbdy_d01_window_<member_id>
                               output: $WORK_DIR/fc/$PREVDATE/wrfinput_<domain_id>_<bin_date>_window1_<member_id>
                                       bin_date=DATE+[OBS_WIN_MIN:MINUTES_PER_SLOT:OBS_WIN_MAX]
                               
  
 Note: A schematic work flow diagram can be found at http://hfip.psu.edu/yxy159/PSU_DA_system_work_flow.jpg 

Utilities:
----------
 jstat                         provides a real-time summary of job status
                               usage: jstat $WORK_DIR

 job_submit.sh                 submit and run jobs (currently supported clusters: stampede and jet)
                               Two possible job submit modes:
                               If $JOB_SUBMIT_MODE=1: run_cycle.sh is submitted into queue, it requestes
                               resources and job_submit.sh executes program. The scheduling (choreograph)
                               of jobs are taken care of in each module.
                               If $JOB_SUBMIT_MODE=2: run_cycle.sh is executed directly, and job_submit.sh
                               creates separate run scripts and submit them. The scheduling is done by the
                               queue.
                                                             
 util.sh                       utility functions including math, time calculation, work flow control, etc.

 util_change_nc_att.ncl        ncl script that changes ncfile global attribute values

 util_linint_nc_time.ncl       ncl script that interpolates a ncfile in time

 calc_domain_moves.sh          calculates the preset move steps for nested inner domains that follows storm, 
                               move is defined by input tcvitals data ($TCVITALS_DIR)

 calc_ij_parent_start.sh       calculates initial nested domain locations in the parent domain according to 
                               tcvitals data

 multi_physics_draw.sh         if multi-physics ensemble is used, this randomly assign physics options in 
                               module_perturb_ic.sh. The following namelist_wrf.sh will read wrfinput files
                               to get options instead of from configuration file.

 multi_physics_reset.sh        for deterministic runs, this resets the physics options to that defined in 
                               the configuration file for wrfinput files.

--------------------------------------------------------------------------------------------------------------

Notes on some specific functionalities:

Boundary condition file
-----------------------
At `DATE_START`, a `wrfbdy_d01` file contains time steps form `DATE_START` to `DATE_END` is created. The update/
perturbation of boundary condition all happen to the same file, but at different time steps (`n_1` option in
WRF_BC_v2.1_alltime). A different copy is made for each member `wrfbdy_d01_<member_id>` and for different
purpose: `wrfbdy_d01` is for `module_wrf_window.sh`, `wrfbdy_d01_window` is for `module_wrf_window1.sh`.

Analysis time is different for EnKF/4DVar
-----------------------------------------
  In `run_cycle.sh` `DATE` starts from `DATE_START` and loops over all assimilation cycles (`DATE_CYCLE_START` to 
`DATE_CYCLE_END`). The analysis time for EnKF is `$DATE` and for 4DVar is `$DATE+$OBS_WIN_MIN`. [`OBS_WIN_MIN` 
`OBS_WIN_MAX`] defines the window in which observations are assimilated for the current cycle. For example,
if the window is [-3 3] and `DATE`=6, the EnKF component will assimilate data from 3 to 9 and assume they are
all valid at 6, while 4DVar will assimilate observations 3, 4, ..., 9 (7 bins) so that observations are
assimilated at the correct time. The bin size can be adjusted by `MINUTES_PER_SLOT(1)`. All numbers in hour 
units in this paragraph, but in minutes in configuration files.

  For hybrid method, extra wrf runs are necessary to match the different analysis times:
  1) previous `module_wrf_ens.sh` will save ensemble at `DATE`+`OBS_WIN_MIN` and calculate mean -> fg
  2) `module_wrf_window1.sh` integrates fg to `DATE+OBS_WIN_MAX` -> fg02
  3) `module_wrf_window.sh` integrates 4DVar analysis at `DATE`+`OBS_WIN_MIN` to `DATE` -> replace ensemble mean

E4DVar vs 4DEnVar
-----------------
  The difference between 4DEnVar and E4DVar is that 4DEnVar uses a set of ensemble perturbations to approximate
U(t)=M(t)U, therefore not using tangent linear and adjoint models. To generate ensemble perturbations, an 
extra ensemble forecast is run through the observation window (`module_wrf_ens_window1.sh`). E4DVar only need
ep at the beginning of observation window, while 4DEnVar needs ep for all observation bins.

Nested domains
--------------
  Currently the nested domains performs data assimilation separately, there could be discontinuity created
by assimilation different amount of observations. We hope that the model will damp out the discontinuity
just like the assimilation-induced imbalances. But there are still on-going development for solving this problem.

Moving (storm-following) nested domains
---------------------------------------
  Data assimilation is performed for each domain separately, for 4DVar the domains should not be moving, thus
fg and fg02 (and possible ep files) should come from fixed domain runs. For 4DEnVar, the eps should also come
from fixed domain runs. On the other hand, the ensemble forecast and deterministic runs from analysis should 
still be moving domain runs. In the schematic diagram (in 1), pink denotes moving domain runs and yellow denote 
fixed domain runs.

CPU choreography
----------------
  Ensemble forecast requires `NUM_ENS`\*`wrf_ntasks`, all cpus on a node are utilized: `wrf_ppn=HOSTPPN`. If 
`run_cycle.sh` header specifies `total_ntasks=NUM_ENS*wrf_ntasks`, all runs will complete at the same time. 
If `total_ntasks` is smaller, the runs will complete in more than one batchs.

  For EnKF/4DVar, due to larger memory demand we need to use smaller ppn, say `enkf_ppn=4` (`HOSTPPN=16`).
`enkf_ntasks=NMCPU*NICPU*NJCPU` and is usually set to (`NUM_ENS`+1)\*1\*1 (no domain decomposition). If `total_ntasks`
is set to `NUM_ENS*HOSTPPN` (`NUM_ENS` nodes), 3 EnKFs can be run simultaneously with each EnKF occupying 
(`NUM_ENS`+1) CPUs or (`NUM_ENS`+1)/(`HOSTPPN`/`enkf_ppn`) nodes.

  So far, no choreography designed for hybrid runs, since 4DVar and EnKF can be run simultaneously, it's hard 
to avoid conflicts over resources. So recommend using `JOB_SUBMIT_MODE=2` for hybrid runs.

--------------------------------------------------------------------------------------------------------------

Installation tips and test case:

1. Unpackage and compile all necessary code packages

2. Modify `run_cycle.sh` header to assign CPU resources. The header should conform to your system's requirements.
   Modify the `total_ntasks` part (the name of system variable may be different).
   If necessary, modify `job_submit.sh` or even add your own execution lines.

3. Modify configuration files to setup your experiment, change paths to the code packages, data files.
   Some model/DA options are given in the configuration file, if you want to personalize an option that's not
   in the configuration file, go modify the corresponding namelist generator/ module.
   The lowest part assigns CPU usage to different components (WRF runs, EnKF, 4DVar, etc)

4. Submit/execute `run_cycle.sh`

5. Check job status either by looking at standard error/output from `run_cycle.sh`, or use `jstat` to get a more
   dynamic look of each component's status

6. Once completed, you will have analysis for each assimilation cycle, then `run_forecast.sh` to run forecast 
   from these analyses.

Test case - Hurricane Katrina (2005)
------------------------------------
A sample configuration file has been setup for the Katrina case: `config/katrina_enkf`

Test case data set can be found here: http://hfip.psu.edu/yxy159/katrinaTestData.tar

Once you completed steps 1 and 2, set `DATA_DIR` to the test data directory, then just submit/execute `run_cycle.sh`.


