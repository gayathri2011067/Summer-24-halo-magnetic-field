module timestepping
    use parameters
    use eta_profile
    use alpha_profile
    use omega_profile
    use initial_field
    use time_grid
    use physical_grid
    use make_a_grid
    use equations

    implicit none
    

    double precision, dimension(nx) :: k1r,k1phi,k2r,k2phi,k3r,k3phi,k4r,k4phi
    contains
        subroutine RK4
            character(len=30) :: ghost_zone_type = 'anti-symmetric'
            call impose_boundary_conditions(B_r, ghost_zone_type)
            call impose_boundary_conditions(B_phi, ghost_zone_type)

            alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)


            call spatial_derivative(B_r,6,dBr,d2Br)
            call spatial_derivative(B_phi,6,dBphi,d2Bphi)
            call spatial_derivative(alpha_cap2,6,d_alpha_cap,d2alpha_cap)

            call diff_equations(B_r,B_phi,dBr,dBphi,d2Br,d2Bphi)
            k1r = dt*dBrdt
            k1phi = dt*dBphidt

            call diff_equations(B_r+0.5*k1r,B_phi+0.5*k1phi,dBr,dBphi,d2Br,d2Bphi)
            k2r = dt*dBrdt
            k2phi = dt*dBphidt

            call diff_equations(B_r+0.5*k2r,B_phi+0.5*k2phi,dBr,dBphi,d2Br,d2Bphi)
            k3r = dt*dBrdt
            k3phi = dt*dBphidt

            call diff_equations(B_r+k3r,B_phi+k3phi,dBr,dBphi,d2Br,d2Bphi)
            k4r = dt*dBrdt
            k4phi = dt*dBphidt

            B_r = B_r + (k1r + 2.*k2r + 2.*k3r + k4r)/6.
            B_phi = B_phi + (k1phi + 2.*k2phi + 2.*k3phi + k4phi)/6.

            t = t + dt
        end subroutine RK4

    


end module timestepping