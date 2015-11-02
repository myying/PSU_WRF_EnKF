begin
a=addfile(ncfile,"w")
att=True
att@SIMULATION_START_DATE=time
fileattdef(a,att)
times=new((/1,20/),"character")
times(0,:)=stringtochar(time)
times!0="Time"
times!1="DateStrLen"
a->Times=times(:,0:18)
end

