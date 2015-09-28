;input: ncfile, attname, attvalue

begin
  a=addfile(ncfile,"w")
  att=True
  att@$attname$=attvalue
  fileattdef(a,att)
end
