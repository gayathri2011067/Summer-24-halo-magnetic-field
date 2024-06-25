module equations
    use parameters
    use eta_profile
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
        dBrdt = -(alpha_fz)*dbphi_dummy -B_phi_dummy*(d_alpha)+ d2Br_dummy 

        dBphidt=R_omega*B_r_dummy + (alpha_fz)*dBr_dummy + B_r_dummy*(d_alpha) + d2Bphi_dummy
    end subroutine diff_equations



 




end module equations