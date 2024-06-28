module initial_field

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use eta_profile
!
  implicit none
!
  double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi, B_eq
  


! 
contains
    subroutine field_initialization
      
        B_r = 0.0001*(1.0-x**2.)*exp(-x**2.)
        B_phi = 0.0
        B_eq = exp(-radius/R - x**2./2.)

    end subroutine field_initialization
!
  end module initial_field