module omega_profile

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use eta_profile
!
  implicit none
!
  double precision, dimension(nx) :: omega_fzr
  


! 
contains
    subroutine construct_omega_profile
      
        ! omega_fzr = om0_kmskpc*(1.+(r_kpc/r0_kpc)**2)**(-1/2)
        omega_fzr = 1

    end subroutine construct_omega_profile
!
  end module omega_profile