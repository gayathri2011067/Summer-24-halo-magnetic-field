module time_grid

  use parameters
  implicit none
  
  integer, parameter :: n1= 800
  integer, parameter :: n2= 1000  !long(1./dt)+1;	;Replot n1 times (after n2 timesteps)
  double precision, parameter :: tsnap= 0.025/td0_Gyr !0.1d0  !Time between successive snapshots
  double precision, parameter :: dt= tsnap/n2 !0.0339/30 !0.0339/3500 !1.d-3!1.d-4!5.d-4  !Timestep in units of td=h^2/C_etat
  double precision :: t=0.
  double precision :: first=0.  !for Runge-Kutta routine
end module time_grid


!**************************************************************************************************************************************************


module physical_grid

  use parameters
  use time_grid
  implicit none

  integer, parameter :: nxphys= 201  !Resolution in z
  integer, parameter :: nxghost= 3  !Number of ghost cells at each end in z
  integer, parameter :: nx= nxphys +2*nxghost  !Resolution in z
  double precision, dimension(nx) :: x
  double precision :: dx
end module physical_grid


!******************************************************************************************************************************************************


module make_a_grid

  implicit none

contains
  subroutine construct_grid
    use parameters
    use time_grid
    use physical_grid

    integer :: i
    double precision, parameter :: len= 2.*h0  !Simulation domain (in units of h0)
    double precision, dimension(nx) :: spac

    dx=len/(nxphys-1)  !x corresponds to z coordinate
    do i=1,nx
      x(i)= -(h0 +nxghost*dx) +(i-1)*dx
    enddo
  endsubroutine construct_grid

end module make_a_grid