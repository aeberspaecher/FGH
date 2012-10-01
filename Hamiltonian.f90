!    -*- f90 -*-

! Construct the Hamiltonian matrix H_ij according to Marston's
! "Fourier Grid Hamiltonian method"

! cmp. BALINT-KURTI et al., INTERNATIONAL REVIEWS IN PHYSICAL CHEMISTRY, VOL. 11, No. 2, 317-344
! (http://www.tandfonline.com/toc/trpc20/11/2)

! The implementations actually used are taken from D. Tannor's book
! "Introdoction to Quantum Mechanics: a time-dependent perspective"

! compile to a Python module with
! f2py --f90flags="-fopenmp -ftree-vectorize -ffree-line-length-none -fmax-errors=2 -Wall -Wextra -Warray-temporaries" -lgomp -m Hamiltonian -c Hamiltonian.f90

subroutine H(N, L, V, Hij)
  ! compute the Hamiltonian NxN matrix in atomic units (hbar = 1; m = 1)
  ! according to Tannor's book, sec. 11.6.4
  implicit none
  integer, intent(in) :: N  ! number of samples for the potential, also: output dimension dim([H_ij]) = N
  double precision, intent(in) :: L  ! length of interval
  double precision, dimension(N), intent(in) :: V  ! sampled potential V(x_i)
  double precision, dimension(N,N), intent(out) :: Hij

  integer :: i, j
  double precision :: K  ! wavenumbers are in [-K, +K]
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)

  K = pi/(L/N)

  do i=1,N
     !$OMP PARALLEL DO PRIVATE(j) SHARED(Hij, i, K, N)
     do j=1,i
        if(i == j) then
           Hij(i, j) = 0.5d0*(K**2/3.0d0*(1.0d0 + 2.0d0/N**2)) + V(i)
        else
           Hij(i, j) = K**2/N**2 * (-1)**(j-i) / (sin(pi*(j-i)/N))**2
           Hij(j, i) = Hij(i, j)  ! use Hermitian symmetry
        end if
     end do
     !$OMP END PARALLEL DO
  end do

end subroutine H

subroutine HCont(N, L, V, Hij)
  ! compute the Hamiltonian NxN matrix in atomic units (hbar = 1; m = 1)
  ! according to Tannor's book, sec. 11.6.4
  ! This routine uses a continuous k variable - thus, no assumptions on periodic
  ! boundary conditions are made
  implicit none
  integer, intent(in) :: N  ! number of samples for the potential, also: output dimension dim([H_ij]) = N
  double precision, intent(in) :: L  ! length of interval
  double precision, dimension(N), intent(in) :: V  ! sampled potential V(x_i)
  double precision, dimension(N,N), intent(out) :: Hij

  integer :: i, j
  double precision :: K  ! wavenumbers are in [-K, +K]
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)

  K = pi/(L/N)

  do i=1,N
     !$OMP PARALLEL DO PRIVATE(j) SHARED(Hij, i, K, N)
     do j=1,i
        if(i == j) then
           Hij(i, j) = 0.5d0*(K**2/3.0d0) + V(i)
        else
           Hij(i, j) = K**2/pi**2 * (-1)**(j-i) / (j-i)**2
           Hij(j, i) = Hij(i, j)  ! use Hermitian symmetry
        end if
     end do
     !$OMP END PARALLEL DO
  end do

end subroutine HCont

