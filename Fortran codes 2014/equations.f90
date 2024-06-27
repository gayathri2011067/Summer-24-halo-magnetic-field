module equations
    use parameters
    use eta_profile
    use velocity_profile
    use alpha_profile
    use omega_profile
    use initial_field
    use time_grid
    use physical_grid
    use make_a_grid

    implicit none
    double precision, dimension(nx) :: dBrdt, dBphidt 
    !   NOTE: Assumed to be predefined: B_r, B_phi,dBr, d2Br, dBphi, d2Bphi, d_alpha
    contains
    subroutine diff_equations(B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy)
      double precision, intent(in), dimension(nx) :: B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy
      dBrdt = R_alpha*B_phi_dummy*d_alpha_cap + R_alpha*alpha_cap*dBphi_dummy &
      + d2Br_dummy - R_U*U_z_cap*dBr_dummy - R_U*B_r_dummy*d_U_z_cap
      dBphidt = R_omega*B_r_dummy + R_alpha*B_r_dummy*d_alpha_cap + R_alpha*alpha_cap*dBr_dummy &
      + d2Bphi_dummy - R_U*U_z_cap*dBphi_dummy - R_U*B_phi_dummy*d_U_z_cap
    end subroutine diff_equations



 




end module equations