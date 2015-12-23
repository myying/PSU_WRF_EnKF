!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program mean
  use constants
  use namelist_define
  use mpi_module
  use mapinfo_define
  use netcdf
  use map_utils
  use wrf_tools
  integer            :: i_unit=80010
  integer            :: unit
  integer            :: ie, ic, it, ix, jx, kx, len, i, j, k, n, ntot, ni,nj,nk
  integer            :: ii, jj, kk, m, nv, nm, gid, sid, iid,jid, g_comm,s_comm
  real               :: member_per_cpu
  character (len=10) :: wrf_file
  real, allocatable, dimension(:,:,:,:)   :: xm
  real, allocatable, dimension(:,:,:,:,:) :: x
  real, allocatable, dimension(:) :: randnum
  type(proj_info)    :: proj
  character (len=80) :: times
  character (len=24) :: finish_file_flg
!----------------------------------------------------------------
! Initialize parallel stuff
  call parallel_start()
! Get WRF model information from the 1th member
! ( reference to MM5_to_GrADS: module_wrf_to_grads_util.f)
  write(wrf_file,'(a5,i5.5)')'fort.',i_unit+1
  call get_wrf_info(wrf_file, ix, jx, kx, times, proj) ! saved wrf info to proj 
! Initilize and Read in namelist information
  call read_namelist(ix, jx, kx) 
! Decompose domain based on n_cpus (nprocs) and n_members (numbers_en) 
! figure out parallel strategy:
if(.not.manual_parallel) then
  nmcpu=1
  nicpu=1
  njcpu=1
  do !nmcpu as close to numbers_en as possible
    nmcpu=nmcpu*2
    if(nmcpu*2>(numbers_en+1) .or. nmcpu*2>nprocs) exit
  end do
  do
    if ( nprocs/(nmcpu*nicpu*njcpu*2) < 1 ) exit
    if ( nicpu <= njcpu ) then
      nicpu=nicpu*2
    else
      njcpu=njcpu*2
    endif
  end do
endif
! nicpu njcpu and nmcpu are from namelist if manual_parallel=true
member_per_cpu=real(numbers_en+1)/nmcpu
nm=ceiling(member_per_cpu)
if ( my_proc_id == 0 ) then
  write(*,*) '---------------------------------------------------'
  write(*,'(a,i4,a,i5,a)') ' PARALLEL SCHEME: ',numbers_en,' ensemble members + 1 mean, ',nprocs,' cpus'
  write(*,'(a,i3,a,i3,a,i2,a)') ' domain decomposed to (ni*nj)=(',nicpu,'*',njcpu,') slabs'
  write(*,'(a,i3,a,f5.2,a,i3)') ' members distributed to ',nmcpu, ' cpu groups, member_per_cpu =', member_per_cpu, ', nm=',nm
  if ((nmcpu*nicpu*njcpu)/=nprocs) stop '*****PARALLEL ERROR: nprocs!=nmcpu*nicpu*njcpu *****'
  if (member_per_cpu<1) stop '*****PARALLEL ERROR: member_per_cpu cannot be less than 1! *****'
endif

! determine which portion of grid my_proc_id is in charge of
gid=int(my_proc_id/(nicpu*njcpu))
sid=mod(my_proc_id,nicpu*njcpu)
iid=mod(sid,nicpu)
jid=int(sid/nicpu)
call MPI_Comm_split(comm, gid, sid, s_comm, ierr)
call MPI_Comm_split(comm, sid, gid, g_comm, ierr)
!write(*,'(a,5i3.2,a)') 'pid,sid,iid,jid,gid=', my_proc_id,sid,iid,jid,gid,'     '

!-- allocate wrf variables, 2d->x2, 3d->x
nv = 0
do m=1,30 
  if(len_trim(enkfvar(m))>=1) nv=nv+1
enddo
ni=int((ix+1)/nicpu)+1
nj=int((jx+1)/njcpu)+1
nk=kx+1
ntot=ni*nj*nk*nv
if(my_proc_id==0) then
   write(*,'(a,i4,a,i4,a,i4,a,i3,a,i3,a,i10)') ' On each cpu: ni*nj*nk*nv*nm=',ni,'*',nj,'*',nk,'*',nv,'*',nm,'=',ntot*nm
   if(ntot*nm>50000000) write(*,*) '***** MEMORY WARNING: too much memory to be used, suggest use less cpu per node! *****'
  write(*,*) '---------------------------------------------------'
  write(*,*) ' '
endif
allocate(x(ni,nj,nk,nv,nm))
allocate(xm(ni,nj,nk,nv))
x=0.
xm=0.

! Read in initial Ensemble in x
if ( my_proc_id == 0 ) write(*,*)'Read in initial Ensemble in x'
call read_ensemble(i_unit,ix,jx,kx,ni,nj,nk,nv,nm,gid,sid,iid,jid,x)

! Calculate the ensemble mean and output
call MPI_Allreduce(sum(x,5),xm,ntot,MPI_REAL,MPI_SUM,g_comm,ierr)
xm=xm/real(numbers_en)
if(gid==0) &
  call output(i_unit+numbers_en+1,ix,jx,kx,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,xm,times)

! Give a finish flag
if ( my_proc_id == 0 ) write(*,'(a)')' Successful completion of EnKF.'
finish_file_flg = times(1:4)//times(6:7)//times(9:10)//times(12:13)//times(15:16)//'.finish_flag'
do k = 1, 12
  if ( finish_file_flg(k:k) .eq. " " ) finish_file_flg(k:k) ='0'
enddo
if ( my_proc_id == 0 ) then
  open(10,file=finish_file_flg)
  write(10,*)times
  close(10)
endif

! Clean up
deallocate(x)
deallocate(xm)
call MPI_Comm_free(g_comm,ierr)
call MPI_Comm_free(s_comm,ierr)
call parallel_finish()

end program mean
