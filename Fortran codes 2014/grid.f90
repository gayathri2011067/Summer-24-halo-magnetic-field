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

  integer, parameter :: nxphys= 51 !Resolution in z
  integer, parameter :: nxghost= 3  !Number of ghost cells at each end in z
  integer, parameter :: nx= nxphys +2*nxghost  !Resolution in z
  double precision, dimension(nx) :: x
  double precision :: dx
end module physical_grid
!TODO_LATER: Make the number of gz automatic wrt fd order.

!******************************************************************************************************************************************************


module make_a_grid

  implicit none

contains
  subroutine construct_grid
    use parameters
    use time_grid
    use physical_grid

    integer :: i
    double precision, parameter :: len= 2.*h  
    double precision, dimension(nx) :: spac

    dx=len/(nxphys-1)  !x corresponds to z
    do i=1,nx
      x(i)= -(h +nxghost*dx) +(i-1)*dx +0.01 !NOTE:+0.01 to avoid 0
      x(i)= x(i) !dimensionless now
    enddo
  endsubroutine construct_grid

end module make_a_grid