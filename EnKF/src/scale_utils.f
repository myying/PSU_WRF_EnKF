module scale_utils

contains

function scale_response(k,krange,s) result(r)
  integer :: s
  real,dimension(:,:) :: k
  integer,dimension(:) :: krange
  real,dimension(size(k,1),size(k,2)) :: r
  ns=size(krange)+1
  if(s<1 .or. s>ns) then
    print *,'band not within range!'
    stop
  end if
  r=0.0
  if(s==1) then
    where(k<krange(s)) r=1.0
  else if(s==ns) then
    where(k>=krange(s-1)) r=1.0
  else
    where(k>=krange(s-1) .and. k<krange(s)) r=1.0
  end if
end function scale_response

function scale_response1(k,krange,s) result(r)
  integer :: s
  real,dimension(:,:) :: k
  integer,dimension(:) :: krange
  real,dimension(size(k,1),size(k,2)) :: r
  ns=size(krange)+1
  if(s<1 .or. s>ns) then
    print *,'band not within range!'
    stop
  end if
  r=0.0
  if(s==ns) then
    where(k>=krange(s-1)) r=1.0
  else
    where(k>=krange(s)) r=1.0
  end if
end function scale_response1

subroutine scale_bandpass(u,krange,s,us,us1)
  use,intrinsic :: iso_c_binding
  use fft_pack
  real,dimension(:,:) :: u
  integer,dimension(:) :: krange
  real,dimension(size(u,1),size(u,2)),intent(out) :: us,us1
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
    uspec=fft2(dcmplx(u,0.0))/(nx*ny)
    flt=scale_response(keff,krange,s)
    us=real(ifft2(uspec*flt))
    flt=scale_response1(keff,krange,s)
    us1=real(ifft2(uspec*flt))
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

end module scale_utils
