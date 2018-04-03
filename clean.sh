#!/bin/bash

for f in CAM_ABS_DATA CAM_AEROPT_DATA ETAMPNEW_DATA ETAMPNEW_DATA.expanded_rain ETAMPNEW_DATA.expanded_rain_DBL ETAMPNEW_DATA_DBL GENPARM.TBL LANDUSE.TBL MPTABLE.TBL README.namelist README.tslist RRTMG_LW_DATA RRTMG_LW_DATA_DBL RRTMG_SW_DATA RRTMG_SW_DATA_DBL RRTM_DATA RRTM_DATA_DBL SOILPARM.TBL URBPARM.TBL URBPARM_UZE.TBL VEGPARM.TBL co2_trans grib2map.tbl gribmap.txt ndown.exe nup.exe ozone.formatted ozone_lat.formatted ozone_plev.formatted real.exe tc.exe tr49t67 tr49t85 tr67t85 wrf.exe namelist.output; do
  rm $f
done

rm rsl.*
