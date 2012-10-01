!    -*- f90 -*-

! Construct the Hamiltonian matrix H_ij according to Marston's
! "Fourier Grid Hamiltonian method"

! cmp. BALINT-KURTI et al., INTERNATIONAL REVIEWS IN PHYSICAL CHEMISTRY, VOL. 11, No. 2, 317-344
! (http://www.tandfonline.com/toc/trpc20/11/2)

! compile to a Python module with
! f2py --f90flags="-ftree-vectorize -ffree-line-length-none -fmax-errors=2 -Wall -Wextra -Warray-temporaries" -m Hamiltonian -c Hamiltonian.f90

subroutine H(N, L, V, Hij)
  ! compute the Hamiltonian NxN matrix in atomic units (hbar = 1; m = 1)
  implicit none
  integer, intent(in) :: N  ! number of samples for the potential, also: output dimension dim([H_ij]) = N
  double precision, intent(in) :: L  ! length of interval
  double precision, dimension(N), intent(in) :: V  ! sampled potential V(x_i)
  double precision, dimension(N,N), intent(out) :: Hij

  integer :: i, j
  double precision :: diagConst, Nreal
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)

  ! note that we chose hbar to be 1, not h. thus, we have a few h = 2pi factors

  Nreal = real(N, kind=kind(1.0d0))
  diagConst = (2*pi)**2/(4.0d0*L**2)*( (Nreal-1.0d0)*(Nreal-2.0d0)/6.0d0 + Nreal/2.0d0 )  ! note the N/2 instead of the 1 in the reference

  do i=1,N
     do j=1,i
        if(i == j) then
           Hij(i, j) = diagConst + V(i)
        else
           Hij(i, j) = (-1)**(i-j) * (2*pi/(2.0d0*L*sin(pi*(i-j)/N)))**2
           Hij(j,i) = Hij(i,j)  ! use Hermitian symmetry
        end if
     end do
  end do

end subroutine H
