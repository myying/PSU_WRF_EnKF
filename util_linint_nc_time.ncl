load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;linearly interpolate the input ncfile to given delta minute (dmin)
;ncl linint_time_ncfile.ncl dmin=* 'infile="*"'
begin
  system("rm -f tmp.nc")
  a=addfile(infile,"r")
  b=addfile("tmp.nc","c")
  varname=(/"SEAICE","VEGFRA","LAI","ALBBCK","SST"/)
  times=a->Times
  dim=dimsizes(times)
  nt=dim(0)
  mm=new(nt,"integer")
  do t=0,nt-1
    mm(t)=stringtointeger(systemfunc("date -u -d '"+charactertostring(times(t,0:9))+" "+charactertostring(times(t,11:15))+"' +%s"))
  end do
  mm=mm-mm(0)
  mm=mm/60
  mm1=ispan(mm(0),mm(nt-1),dmin)
  print(mm)
  print(mm1)
  times1=new((/dimsizes(mm1),19/),"character")
  do t=0,dimsizes(mm1)-1
    tstring=stringtocharacter(systemfunc("date -u -d '"+mm1(t)+" minute "+charactertostring(times(0,0:9))+" "+charactertostring(times(0,11:15))+"' +%Y-%m-%d_%H:%M:00"))
    times1(t,:)=tstring(0:18)
  end do
  times1!0="Time"
  times1!1="DateStrLen"
  b->Times=times1

  do v=0,dimsizes(varname)-1
    data=a->$varname(v)$
    time=new(dimsizes(data),"integer")
    do t=0,nt-1
      time(t,:,:)=mm(t)
    end do
    data1=linint1_n_Wrap(time,data,False,mm1,0,0)
    data1!0="Time"
    b->$varname(v)$=data1
    fileattdef(b,a)
    delete(data1)
  end do

end
