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
    double precision, dimension(nx) :: dBrdt, dBphidt !,d_alpha_Bphi, d_alpha_Br, d_Uz_Bphi, d2_Uz_Br
    !   NOTE: Assumed to be predefined: B_r, B_phi,dBr, d2Br, dBphi, d2Bphi, d_alpha
    contains
    ! subroutine diff_equations_split(B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy)
    !   double precision, intent(in), dimension(nx) :: B_r_dummy, B_phi_dummy, dBr_dummy, d2Br_dummy, dBphi_dummy, d2Bphi_dummy
    !   dBrdt = R_alpha*B_phi_dummy*d_alpha_cap + R_alpha*alpha_cap2*dBphi_dummy &
    !   + d2Br_dummy - R_U*U_z_cap*dBr_dummy - R_U*B_r_dummy*d_U_z_cap
    !   dBphidt = R_omega*B_r_dummy + R_alpha*B_r_dummy*d_alpha_cap + R_alpha*alpha_cap2*dBr_dummy &
    !   + d2Bphi_dummy - R_U*U_z_cap*dBphi_dummy - R_U*B_phi_dummy*d_U_z_cap

    ! end subroutine diff_equations_split
   
    subroutine diff_equations_no_split(B_r_dummy, B_phi_dummy)
      double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy
      double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
      double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
      double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
      double precision,  dimension(nx) :: Uz_Bphi, Uz_Br
      double precision,  dimension(nx) :: alpha_Bphi, alpha_Br

      character(len=30) :: ghost_zone_type = 'anti-symmetric'
      call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
      call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)

      alpha_cap2 = alpha_cap / (1.0 + (B_r_dummy**2 + B_phi_dummy**2) / B_eq**2)
      alpha_Bphi = alpha_cap2 * B_phi_dummy
      alpha_Br = alpha_cap2 * B_r_dummy
      Uz_Br = B_r_dummy * U_z_cap
      Uz_Bphi = B_phi_dummy * U_z_cap

      call spatial_derivative(B_r_dummy, 2, dBr, d2Br)
      call spatial_derivative(B_phi_dummy, 2, dBphi, d2Bphi)
      call spatial_derivative(alpha_Br, 2, d_alpha_Br, d2_alpha_Br)
      call spatial_derivative(alpha_Bphi, 2, d_alpha_Bphi, d2_alpha_Bphi)
      call spatial_derivative(Uz_Br, 2, d_Uz_Br, d2_Uz_Br)
      call spatial_derivative(Uz_Bphi, 2, d_Uz_Bphi, d2_Uz_Bphi)


      dBrdt = - R_alpha*d_alpha_Bphi &
      + d2Br - R_U*d_Uz_Br
      dBphidt = R_omega*B_r_dummy + R_alpha*d_alpha_Br &
      + d2Bphi - R_U*d_Uz_Bphi

    end subroutine diff_equations_no_split
  

 

     


 




end module equations