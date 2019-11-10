module fft_pack
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'

type(C_PTR) :: plan_f,plan_b

contains

subroutine init_fft(nx,ny)
  integer,intent(in) :: nx,ny
  complex(C_DOUBLE_COMPLEX),dimension(:,:),allocatable :: f,fr
  allocate(f(nx,ny),fr(nx,ny))
  plan_f=fftw_plan_dft_2d(nx,ny,f,fr,FFTW_FORWARD,FFTW_PATIENT)
  plan_b=fftw_plan_dft_2d(nx,ny,f,fr,FFTW_BACKWARD,FFTW_PATIENT)
  deallocate(f,fr)
end subroutine init_fft

function fft2(f) result(fr)
  complex(C_DOUBLE_COMPLEX), dimension(:,:) :: f
  complex(C_DOUBLE_COMPLEX), dimension(size(f,1),size(f,2)) :: fr
  call fftw_execute_dft(plan_f,f,fr)
end function fft2

function ifft2(f) result(fr)
  complex(C_DOUBLE_COMPLEX), dimension(:,:) :: f
  complex(C_DOUBLE_COMPLEX), dimension(size(f,1),size(f,2)) :: fr
  call fftw_execute_dft(plan_b,f,fr)
end function ifft2

!function fullspec(hspec) result(fspec)
!  complex(C_DOUBLE_COMPLEX),dimension(:,:) :: hspec
!  complex(C_DOUBLE_COMPLEX),dimension(size(hspec,1)+1,2*size(hspec,2)) :: fspec
!  kmax=size(hspec,2)-1
!  fspec=0.
!  hspec(1:kmax+1,1)=conjg(hspec(2*kmax+1:kmax+1:-1,1))
!  fspec(2:2*kmax+2,kmax+2:2*kmax+2)=hspec
!  fspec(2*kmax+2:2:-1,kmax+2:2:-1)=conjg(hspec)
!end function fullspec

!function halfspec(fspec) result(hspec)
!  complex(C_DOUBLE_COMPLEX),dimension(:,:) :: fspec
!  complex(C_DOUBLE_COMPLEX),dimension(size(fspec,1)-1,size(fspec,2)/2) :: hspec
!  kmax=size(fspec,1)/2-1
!  hspec=fspec(2:2*kmax+2,kmax+2:2*kmax+2)
!  hspec(1:kmax,1)=0.
!end function halfspec

function fftshift(f) result(fr)
  complex(C_DOUBLE_COMPLEX),dimension(:,:) :: f
  complex(C_DOUBLE_COMPLEX),dimension(size(f,1),size(f,2)) :: fr
  integer,dimension(size(f,1)) :: xind
  integer,dimension(size(f,2)) :: yind
  integer :: nx,ny,i,j
  nx=size(f,1)
  ny=size(f,2)
  xind=(/ (i, i=ceiling(real(nx)/2)+1,nx), (i, i=1,ceiling(real(nx)/2)) /)
  yind=(/ (j, j=ceiling(real(ny)/2)+1,ny), (j, j=1,ceiling(real(ny)/2)) /)
  fr=f(xind,yind)
end function fftshift

function ifftshift(f) result(fr)
  complex(C_DOUBLE_COMPLEX),dimension(:,:) :: f
  complex(C_DOUBLE_COMPLEX),dimension(size(f,1),size(f,2)) :: fr
  integer,dimension(size(f,1)) :: xind
  integer,dimension(size(f,2)) :: yind
  integer :: nx,ny,i,j
  nx=size(f,1)
  ny=size(f,2)
  xind=(/ (i, i=floor(real(nx)/2)+1,nx), (i, i=1,floor(real(nx)/2)) /)
  yind=(/ (j, j=floor(real(ny)/2)+1,ny), (j, j=1,floor(real(ny)/2)) /)
  fr=f(xind,yind)
end function ifftshift

end module fft_pack

