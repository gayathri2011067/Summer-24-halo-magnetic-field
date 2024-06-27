module velocity_profile

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use omega_profile
  use spatial_derivatives
!
  implicit none
!
 
  double precision, dimension(nx) :: U_z_cap,U_z, d_U_z_cap

 

! 
contains
    subroutine construct_velocity_profile
        character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
        U_z_cap = x
        U_z = R_U*U_z_cap
        d_U_z_cap = 1.

    end subroutine construct_velocity_profile
!
  end module velocity_profile