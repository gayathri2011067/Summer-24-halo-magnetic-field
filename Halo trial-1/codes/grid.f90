module time_grid

  use parameters
  implicit none
  
  integer, parameter :: Nt= 5000  !points per diffusion time
  integer, parameter :: n1= 500  !Number of snapshots
  double precision, parameter :: total_t= 27. !unit diffusion time
  integer :: n2 = total_t*Nt/n1!Number of timesteps between snapshots
  double precision, parameter :: dt= 1./Nt !time step
  double precision :: t=0.
  double precision :: first=0.  !for Runge-Kutta routine


end module time_grid


!**************************************************************************************************************************************************


module physical_grid

  use parameters
  use time_grid
  implicit none

  integer, parameter :: nxphys = 276 !CHANGE: old value was 51 here and 225 as a seperate nxvacuum
  !integer, parameter :: nxvacuum = 225  !Resolution in z for the vacuum halo region on one side !CHANGED!
  integer, parameter :: nxghost = 3  !Number of ghost cells at each end in z
  integer, parameter :: nx = nxphys + 2*nxghost !+ 2*nxvacuum  !CHANGE: no more nxvacuum
  double precision, dimension(nx) :: x
  double precision :: dx, alp, beta
end module physical_grid
!TODO_LATER: Make the number of gz automatic wrt fd order.

!******************************************************************************************************************************************************


module make_a_grid

  implicit none
  integer :: i
  contains
  subroutine construct_grid
    use parameters
    use time_grid
    use physical_grid

    ! integer :: i

    double precision, parameter :: len = 20.*h  !CHANGED! !CHANGE: kept the new change
    double precision, dimension(nx) :: spac

    dx=len/(nxphys-1)  !x corresponds to z !CHANGE: old: dx=len/(nxphys+nxvacuum-1)
    alp=dt/dx**2  !alpha for finite difference scheme
    beta= dt/dx !beta for fin diff 

    do i=1,nx
      x(i)= -(10.*h +nxghost*dx) +(i-1)*dx +0.01 !NOTE:+0.01 to avoid 0 !!CHANGE: kept the new change
      x(i)= x(i) !dimensionless now
    enddo
  endsubroutine construct_grid

end module make_a_grid