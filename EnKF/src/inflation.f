  program inflation
  use constants
  use namelist_define
  use mpi_module
  use mapinfo_define
  use netcdf
  use obs_define
  use map_utils
  use wrf_tools
  use radar
  integer            :: u_unit=50010, i_unit=90010, o_unit=70010
  integer            :: unit
  integer            :: iob, ie, ic, it, ix, jx, kx, len, i, j, k, n, ntot, ni,nj,nk
  integer            :: ii, jj, kk, m, nv, nm, gid, sid, iid,jid, g_comm,s_comm
  real               :: member_per_cpu
  character (len=10) :: wrf_file
  double precision   :: time1, time2
  real, allocatable, dimension(:,:,:,:)   :: xm, xom, inf_mean, inf_sd
  real, allocatable, dimension(:,:,:,:,:) :: x, xo
  real, allocatable, dimension(:) :: randnum
  integer,allocatable,dimension(:) :: ind
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
  write(wrf_file,'(a5,i5.5)')'fort.',i_unit+numbers_en+1
  call get_all_obs(wrf_file, ix, jx, kx, times, proj)
  if ( print_detail > 1 .and. my_proc_id == 0 ) write(*,*)'All observations are loaded'

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
allocate(inf_mean(ni,nj,nk,nv))
allocate(inf_sd(ni,nj,nk,nv))
allocate(x(ni,nj,nk,nv,nm))
allocate(xm(ni,nj,nk,nv))
allocate(xo(ni,nj,nk,nv,nm))
allocate(xom(ni,nj,nk,nv))
x=0.
xm=0.
xo=0.
xom=0.
inf_mean=1.
inf_sd=1.

! Read in initial Ensemble in x
call read_ensemble(i_unit,ix,jx,kx,ni,nj,nk,nv,nm,gid,sid,iid,jid,xo)
call read_ensemble(u_unit,ix,jx,kx,ni,nj,nk,nv,nm,gid,sid,iid,jid,x)

! Calculate the ensemble mean and output
call MPI_Allreduce(sum(x,5),xm,ntot,MPI_REAL,MPI_SUM,g_comm,ierr)
call MPI_Allreduce(sum(xo,5),xom,ntot,MPI_REAL,MPI_SUM,g_comm,ierr)
xm=xm/real(numbers_en)
xom=xom/real(numbers_en)
do n=1,nm
  ie=(n-1)*nmcpu+gid+1
  if(ie==numbers_en+1) then
    x(:,:,:,:,n)=xm
    xo(:,:,:,:,n)=xom
  endif
enddo
if(gid==0) then
  call output(i_unit+numbers_en+1,ix,jx,kx,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,xom,times)
  call output(u_unit+numbers_en+1,ix,jx,kx,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,xm,times)
endif

allocate(ind(obs%num))
if(my_proc_id==0) then
  do it=1,obs%num
    ind(it)=it
  enddo
  if(random_order) then
    allocate(randnum(obs%num))
    call random_seed()
    call random_number(randnum)
    call quicksort(obs%num,randnum,ind)
    deallocate(randnum)
  endif
endif
call MPI_Bcast(ind,obs%num,MPI_INTEGER,0,comm,ierr)

write(wrf_file,'(a5,i5.5)')'fort.',i_unit+1
call update_inflation(wrf_file,ix,jx,kx,ni,nj,nk,nv,nm,g_comm,s_comm,gid,sid,iid,jid,xom,xo,xm,x,inf_mean,inf_sd,ind,proj,times)
if(gid==0) then
  call output(20011,ix,jx,kx,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,inf_mean,times)
endif

deallocate(ind,xo,xom,inf_mean,inf_sd)

! Save the Analyses
if ( my_proc_id == 0 ) write(*,*)'Output members and mean now...'
do n =1, nm
  ie = (n-1)*nmcpu+gid+1
  if ( ie <= numbers_en+1 ) then
    call output(o_unit+ie,ix,jx,kx,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,x(:,:,:,:,n),times)
  endif
enddo

! Clean up
deallocate(x,xm)
call MPI_Comm_free(g_comm,ierr)
call MPI_Comm_free(s_comm,ierr)
call parallel_finish()
if ( my_proc_id == 0 ) write(*,'(a)')' Successful completion of inflation.'

end program inflation
