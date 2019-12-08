module multiscale_utils

contains

subroutine scale_response(k,krange,s,r)  !!!!response func for scale s
  integer :: s
  real,dimension(:,:) :: k
  real,dimension(:) :: krange
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
  real,dimension(:) :: krange
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
  real,dimension(:) :: krange
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
  integer :: i1,i2,j1,j2,ni,nj
  ni=size(x,1); nj=size(x,2)
  i1=floor(i); i2=i1+1; di=i-real(i1)
  j1=floor(j); j2=j1+1; dj=j-real(j1)
  xout=(1-di)*(1-dj)*x(i1,j1)+di*(1-dj)*x(i2,j1)+(1-di)*dj*x(i1,j2)+di*dj*x(i2,j2)
end function interp2d


subroutine optical_flow_HS(xb,xa,u,v,n)
  real :: w=100
  integer :: niter=100,nlevel=5,buffer=4
  integer :: i,j,n,nn,lev,itr
  real,dimension(n,n),intent(in) :: xb,xa
  real,dimension(n,n),intent(inout) :: u,v
  real,dimension(n,n) :: x1,x2,Im1c,Im2c,Ix,Iy,It,du,dv,ubar,vbar
  !print*,size(xb,1),size(xb,2)
  u=0.; v=0.
  x1=xb
  do lev=nlevel+1,1,-1
    x2=x1
    do i=1+buffer,n-buffer
    do j=1+buffer,n-buffer
      x1(i,j) = interp2d(x2,real(i)-u(i,j),real(j)-v(i,j))
    enddo
    enddo
    Im1c = coarsen(x1,lev)
    Im2c = coarsen(xa,lev)
    Ix=0.5*(deriv_x(Im1c)+deriv_x(Im2c))
    Iy=0.5*(deriv_y(Im1c)+deriv_y(Im2c))
    It=Im2c-Im1c
    du=0.; dv=0.
    do itr=1,niter
      du(1,:)=0.; du(nn,:)=0.; du(:,1)=0.; du(:,nn)=0. !!boundary condition
      dv(1,:)=0.; dv(nn,:)=0.; dv(:,1)=0.; dv(:,nn)=0.
      ubar=laplacian(du)+du
      vbar=laplacian(dv)+dv
      du=ubar-Ix*(Ix*ubar+Iy*vbar+It)/(w+Ix**2+Iy**2)
      dv=vbar-Iy*(Ix*ubar+Iy*vbar+It)/(w+Ix**2+Iy**2)
    enddo
    u=u+du
    v=v+dv
    u(1,:)=0.; u(n,:)=0.; u(:,1)=0.; u(:,n)=0. !!boundary condition
    v(1,:)=0.; v(n,:)=0.; v(:,1)=0.; v(:,n)=0.
  enddo
end subroutine optical_flow_HS

function coarsen(Im,lev) result(Imc)
  integer :: i,n,lev
  real,dimension(:,:) :: Im
  real,dimension(size(Im,1),size(Im,2)) :: Imc
  Imc=Im
end function coarsen

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
