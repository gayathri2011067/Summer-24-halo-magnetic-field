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
    double precision, dimension(nx) :: dBrdt, dBphidt ,d_alpha_Bphi, d_alpha_Br, d_Uz_Bphi, d2_Uz_Br
    !   NOTE: Assumed to be predefined: B_r, B_phi,dBr, d2Br, dBphi, d2Bphi, d_alpha
    contains
    subroutine diff_equations_split(B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy)
      double precision, intent(in), dimension(nx) :: B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy
      dBrdt = R_alpha*B_phi_dummy*d_alpha_cap + R_alpha*alpha_cap2*dBphi_dummy &
      + d2Br_dummy - R_U*U_z_cap*dBr_dummy - R_U*B_r_dummy*d_U_z_cap
      dBphidt = R_omega*B_r_dummy + R_alpha*B_r_dummy*d_alpha_cap + R_alpha*alpha_cap2*dBr_dummy &
      + d2Bphi_dummy - R_U*U_z_cap*dBphi_dummy - R_U*B_phi_dummy*d_U_z_cap

    end subroutine diff_equations_split

    subroutine diff_equations_no_split(B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy, d_alpha_Bphi_dummy, d_alpha_Br_dummy, d_Uz_Bphi_dummy, d_Uz_Br_dummy)
      double precision, intent(in), dimension(nx) :: B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy, d_alpha_Bphi_dummy, d_alpha_Br_dummy, d_Uz_Bphi_dummy, d_Uz_Br_dummy
      dBrdt = R_alpha*d_alpha_Bphi_dummy &
      + d2Br_dummy - R_U*d_Uz_Br_dummy
      dBphidt = R_omega*B_r_dummy + R_alpha*d_alpha_Br_dummy &
      + d2Bphi_dummy - R_U*d_Uz_Bphi_dummy

    end subroutine diff_equations_no_split



 




end module equations