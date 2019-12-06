module multiscale_utils

contains

subroutine scale_response(k,krange,s,r)  !!!!response func for scale s
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
    where(k<krange(1)) r=1.0
  else if(s==ns) then
    where(k>=krange(ns-1)) r=1.0
  else
    where(k>=krange(s-1) .and. k<krange(s)) r=1.0
  end if
end subroutine scale_response

subroutine scale_response1(k,krange,s,r)  !!!scales smaller than s
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
  if(s<ns) then
    where(k>=krange(s)) r=1.0
  end if
end subroutine scale_response1

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
  ns=size(krange)+1
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
    call scale_response(keff,krange,s,flt)
    us=real(ifft2(uspec*flt))
    call scale_response1(keff,krange,s,flt)
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


function interp2d(x,i,j) result(xout)
  real,dimension(:,:),intent(in) :: x
  real :: i,j,xout,di,dj
  integer :: i1,i2,j1,j2
  i1=floor(i); i2=i1+1; di=i-real(i1)
  j1=floor(j); j2=j1+1; dj=j-real(j1)
  xout=(1-di)*(1-dj)*x(i1,j1)+di*(1-dj)*x(i2,j1)+(1-di)*dj*x(i1,j2)+di*dj*x(i2,j2)
end function interp2d

subroutine optical_flow_HS(xb,xa,u,v)
  real :: w1=100, w2=100
  real,dimension(:,:),intent(in) :: xb,xa
  real,dimension(:,:),intent(inout) :: u,v
  real,dimension(size(xb,1),size(xb,2)) :: Ix,Iy,It,ubar1,vbar1,ubar2,vbar2,uxy,vxy
  integer :: niter,nx,ny
  nx=size(xb,1); ny=size(xb,2)
  niter=100
  Ix=0.5*(deriv_x(xb)+deriv_x(xa))
  Iy=0.5*(deriv_y(xb)+deriv_y(xa))
  It=xa-xb
  u=0.; v=0.
  do i=1,niter
    u(1,:)=0.; u(nx,:)=0.; u(:,1)=0.; u(:,ny)=0. !!boundary condition
    v(1,:)=0.; v(nx,:)=0.; v(:,1)=0.; v(:,ny)=0.
    ubar2=laplacian(u)+u
    vbar2=laplacian(v)+v
    ubar1=deriv_x(deriv_x(u))+u
    vbar1=deriv_y(deriv_y(v))+v
    uxy=deriv_x(deriv_y(u))
    vxy=deriv_x(deriv_y(v))
    u = (w1*ubar2 + w2*(ubar1+vxy))/(w1+w2) - Ix*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)
    v = (w1*vbar2 + w2*(vbar1+uxy))/(w1+w2) - Iy*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)
  enddo
end subroutine optical_flow_HS

function deriv_x(u) result(dudx)
  real,dimension(:,:) :: u
  real,dimension(size(u,1),size(u,2)) :: dudx
  integer :: nx
  nx=size(u,1)
  dudx(2:nx-1,:)=0.5*(u(3:nx,:)-u(1:nx-2,:))
  dudx(1,:)=u(2,:)-u(1,:)
  dudx(nx,:)=u(nx,:)-u(nx-1,:)
end function deriv_x

function deriv_y(u) result(dudy)
  real,dimension(:,:) :: u
  real,dimension(size(u,1),size(u,2)) :: dudy
  integer :: ny
  ny=size(u,2)
  dudy(:,2:ny-1)=0.5*(u(:,3:ny)-u(:,1:ny-2))
  dudy(:,1)=u(:,2)-u(:,1)
  dudy(:,ny)=u(:,ny)-u(:,ny-1)
end function deriv_y

function laplacian(u) result(lapl_u)
  real,dimension(:,:) :: u
  real,dimension(size(u,1),size(u,2)) :: lapl_u
  lapl_u=deriv_x(deriv_x(u))+deriv_y(deriv_y(u))
end function laplacian

end module multiscale_utils
