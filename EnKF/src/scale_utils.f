module scale_utils

contains

function scale_response(k,krange,s) result(r)
  integer :: s
  real,dimension(:,:) :: k
  integer,dimension(:) :: krange
  real,dimension(size(k,1),size(k,2)) :: r
  ns=size(krange)
  if(s<1 .or. s>ns) then
    print *,'band not within range!'
    stop
  end if

  r=0.0
  if(s==ns) then
    where(k>=krange(s)) r=1.0
  else
    where(k>=krange(s) .and. k<krange(s+1)) r=1.0
  end if
end function scale_response

function scale_response1(k,krange,s) result(r)
  integer :: s,center_k,left_k,right_k
  real,dimension(:,:) :: k
  integer,dimension(:) :: krange
  real,dimension(size(k,1),size(k,2)) :: r
  real,parameter :: pi=3.14159265358979
  ns=size(krange)
  if(s<1 .or. s>ns) then
    print *,'band not within range!'
    stop
  end if

  r=0.0
  center_k=krange(s)
  if(s.eq.1) then
    where(k<=center_k) r=1.0
  else
    left_k=krange(s-1)
    where(k<=center_k .and. k>=left_k) &
      r=cos((k-center_k)*(0.5*pi/(left_k-center_k)))**2
  end if
  if(s.eq.ns) then
    where(k>=center_k) r=1.0
  else
    right_k=krange(s+1)
    where(k>=center_k .and. k<=right_k) &
      r=cos((k-center_k)*(0.5*pi/(right_k-center_k)))**2
  end if
end function scale_response1

!subroutine separate_scales(u,krange,ums)
!  use,intrinsic :: iso_c_binding
!  use fft_pack
!  real,dimension(:,:,:) :: u
!  integer,dimension(:) :: krange
!  real,dimension(size(u,1),size(u,2),size(u,3),size(krange)) :: ums
!  integer :: nx,ny,nz,ns, n,i,z,s
!  integer,dimension(size(u,1)) :: kx1
!  integer,dimension(size(u,2)) :: ky1
!  real,dimension(size(u,1),size(u,2)) :: kx,ky
!  real,dimension(size(u,1),size(u,2)) :: keff,flt
!  complex(C_DOUBLE_COMPLEX),dimension(size(u,1),size(u,2)) :: uspec
!  nx=size(u,1)
!  ny=size(u,2)
!  nz=size(u,3)
!  ns=size(krange)
!  if(ns.eq.1) then
!    ums(:,:,:,1)=u
!  else
!    if(mod(nx,2).eq.0) then
!      kx1=(/(i,i=0,ceiling(real(nx-1)/2)),(i,i=-ceiling(real(nx-1)/2)+1,-1)/)
!    else
!      kx1=(/(i,i=0,ceiling(real(nx-1)/2)),(i,i=-ceiling(real(nx-1)/2),-1)/)
!    end if
!    if(mod(ny,2).eq.0) then
!      ky1=(/(i,i=0,ceiling(real(ny-1)/2)),(i,i=-ceiling(real(ny-1)/2)+1,-1)/)
!    else
!      ky1=(/(i,i=0,ceiling(real(ny-1)/2)),(i,i=-ceiling(real(ny-1)/2),-1)/)
!    end if
!    call grid2d(real(kx1),real(ky1),kx,ky)
!    n=max(nx,ny)
!    keff=sqrt((kx*(n/nx))**2+(ky*(n/ny))**2)
!    do s=1,ns
!      flt=scale_response1(keff,krange,s)
!      do z=1,nz
!        uspec=fft2(dcmplx(u(:,:,z),0.0))/(nx*ny)
!        uspec=uspec*flt
!        ums(:,:,z,s)=real(ifft2(uspec))
!      end do
!    end do
!  end if
!end subroutine separate_scales

subroutine scale_bandpass(u,krange,s,us)
  use,intrinsic :: iso_c_binding
  use fft_pack
  real,dimension(:,:) :: u
  integer,dimension(:) :: krange
  real,dimension(size(u,1),size(u,2)),intent(out) :: us
  integer :: nx,ny,ns, n,i,z,s
  integer,dimension(size(u,1)) :: kx1
  integer,dimension(size(u,2)) :: ky1
  real,dimension(size(u,1),size(u,2)) :: kx,ky
  real,dimension(size(u,1),size(u,2)) :: keff,flt
  complex(C_DOUBLE_COMPLEX),dimension(size(u,1),size(u,2)) :: uspec
  nx=size(u,1)
  ny=size(u,2)
  ns=size(krange)
  call init_fft(nx,ny)
  if(ns.eq.1) then
    us=u
  else
    if(mod(nx,2).eq.0) then
      kx1=(/(i,i=0,ceiling(real(nx-1)/2)),(i,i=-ceiling(real(nx-1)/2)+1,-1)/)
    else
      kx1=(/(i,i=0,ceiling(real(nx-1)/2)),(i,i=-ceiling(real(nx-1)/2),-1)/)
    end if
    if(mod(ny,2).eq.0) then
      ky1=(/(i,i=0,ceiling(real(ny-1)/2)),(i,i=-ceiling(real(ny-1)/2)+1,-1)/)
    else
      ky1=(/(i,i=0,ceiling(real(ny-1)/2)),(i,i=-ceiling(real(ny-1)/2),-1)/)
    end if
    call grid2d(real(kx1),real(ky1),kx,ky)
    n=max(nx,ny)
    keff=sqrt((kx*(n/nx))**2+(ky*(n/ny))**2)
    flt=scale_response(keff,krange,s)
    uspec=fft2(dcmplx(u,0.0))/(nx*ny)
    uspec=uspec*flt
    us=real(ifft2(uspec))
  end if
end subroutine scale_bandpass

subroutine grid2d(xind,yind,x,y)
  real,dimension(:),intent(in) :: xind,yind
  real,dimension(size(xind),size(yind)),intent(out) :: x,y
  integer :: i,j,nx,ny
  nx=size(xind)
  ny=size(yind)
  y=spread(yind,1,nx)
  x=spread(xind,2,ny)
end subroutine grid2d

subroutine grid3d(xind,yind,zind,x,y,z)
  real,dimension(:),intent(in) :: xind,yind,zind
  real,dimension(size(xind),size(yind),size(zind)),intent(out) :: x,y,z
  integer :: i,j,k,nx,ny,nz
  nx=size(xind)
  ny=size(yind)
  nz=size(zind)
  x=spread(spread(xind,2,ny),3,nz)
  y=spread(spread(yind,1,nx),3,nz)
  z=spread(spread(zind,1,nx),2,ny)
end subroutine grid3d

end module scale_utils
