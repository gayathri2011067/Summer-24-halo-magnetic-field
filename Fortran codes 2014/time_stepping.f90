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
    double precision,dimension(nx) :: g_1, g_2, g_3, z_1, z_2
    contains
        ! subroutine RK4
        !     character(len=30) :: ghost_zone_type = 'anti-symmetric'
        !     call impose_boundary_conditions(B_r, ghost_zone_type)
        !     call impose_boundary_conditions(B_phi, ghost_zone_type)

        !     alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)


        !     call spatial_derivative(B_r,6,dBr,d2Br)
        !     call spatial_derivative(B_phi,6,dBphi,d2Bphi)
        !     call spatial_derivative(alpha_cap2,6,d_alpha_cap,d2alpha_cap)

        !     call diff_equations_split(B_r,B_phi,dBr,dBphi,d2Br,d2Bphi)
        !     k1r = dt*dBrdt
        !     k1phi = dt*dBphidt

        !     call diff_equations_split(B_r+0.5*k1r,B_phi+0.5*k1phi,dBr,dBphi,d2Br,d2Bphi)
        !     k2r = dt*dBrdt
        !     k2phi = dt*dBphidt

        !     call diff_equations_split(B_r+0.5*k2r,B_phi+0.5*k2phi,dBr,dBphi,d2Br,d2Bphi)
        !     k3r = dt*dBrdt
        !     k3phi = dt*dBphidt

        !     call diff_equations_split(B_r+k3r,B_phi+k3phi,dBr,dBphi,d2Br,d2Bphi)
        !     k4r = dt*dBrdt
        !     k4phi = dt*dBphidt

        !     B_r = B_r + (k1r + 2.*k2r + 2.*k3r + k4r)/6.
        !     B_phi = B_phi + (k1phi + 2.*k2phi + 2.*k3phi + k4phi)/6.

        !     t = t + dt
        ! end subroutine RK4

        subroutine forward_difference
            double precision, dimension(nx) :: temp_Br, temp_Bphi
            character(len=30) :: ghost_zone_type = 'anti-symmetric'

            call impose_boundary_conditions(B_r, ghost_zone_type)
            call impose_boundary_conditions(B_phi, ghost_zone_type)
            call construct_grid
            alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
            alpha_Bphi = alpha_cap2*B_phi
            alpha_Br = alpha_cap2*B_r
            Uz_Br = B_r*U_z_cap
            Uz_Bphi = B_phi*U_z_cap

            ! call spatial_derivative(B_r,2,dBr,d2Br)
            ! call spatial_derivative(B_phi,2,dBphi,d2Bphi)
            ! call spatial_derivative(alpha_Br,2,d_alpha_Br,d2_alpha_Br)
            ! call spatial_derivative(alpha_Bphi,2,d_alpha_Bphi,d2_alpha_Bphi)
            ! call spatial_derivative(Uz_Br,2,d_Uz_Br,d2_Uz_Br)
            ! call spatial_derivative(Uz_Bphi,2,d_Uz_Bphi,d2_Uz_Bphi)

            do i = nxghost+1,nxphys+nxghost
                temp_Br(i) = B_r(i) - (beta/2)*(alpha_Bphi(i+1)-alpha_Bphi(i-1)) &
                +alp*(B_r(i+1)-2*B_r(i)+B_r(i-1)) - (beta/2)*(Uz_Br(i+1) - Uz_Br(i-1))
                temp_Bphi(i) = B_phi(i) - (beta/2)*(alpha_Br(i+1)-alpha_Br(i-1)) &
                +alp*(B_phi(i+1)-2*B_phi(i)+B_phi(i-1)) - (beta/2)*(Uz_Bphi(i+1) - Uz_Bphi(i-1))&
                +G*B_r(i)

                
            end do
            B_r = temp_Br
            B_phi = temp_Bphi
            t = t + dt
        end subroutine forward_difference

        subroutine backward_difference
            character(len=30) :: ghost_zone_type = 'anti-symmetric'

            double precision, dimension(nx) :: temp_Br, temp_Bphi
            call impose_boundary_conditions(B_r, ghost_zone_type)
            call impose_boundary_conditions(B_phi, ghost_zone_type)
            call construct_grid
            alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
            alpha_Bphi = alpha_cap2*B_phi
            alpha_Br = alpha_cap2*B_r
            Uz_Br = B_r*U_z_cap
            Uz_Bphi = B_phi*U_z_cap

            ! call spatial_derivative(B_r,2,dBr,d2Br)
            ! call spatial_derivative(B_phi,2,dBphi,d2Bphi)
            ! call spatial_derivative(alpha_Br,2,d_alpha_Br,d2_alpha_Br)
            ! call spatial_derivative(alpha_Bphi,2,d_alpha_Bphi,d2_alpha_Bphi)
            ! call spatial_derivative(Uz_Br,2,d_Uz_Br,d2_Uz_Br)
            ! call spatial_derivative(Uz_Bphi,2,d_Uz_Bphi,d2_Uz_Bphi)

            do i = nxghost+1,nxphys+nxghost
                temp_Br(i) = B_r(i) + (beta/2)*(alpha_Bphi(i)-alpha_Bphi(i-2)) &
                -alp*(B_r(i+1)-2*B_r(i)+B_r(i-1)) + (beta/2)*(Uz_Br(i) - Uz_Br(i-2))
                temp_Bphi(i) = B_phi(i) - (beta/2)*(alpha_Br(i)-alpha_Br(i-2)) &
                -alp*(B_phi(i+1)-2*B_phi(i)+B_phi(i-1)) + (beta/2)*(Uz_Bphi(i) - Uz_Bphi(i-2))&
                -G*B_r(i)

            end do
            B_r = temp_Br
            B_phi = temp_Bphi
            t = t + dt

        end subroutine backward_difference

        subroutine central_difference
            double precision, dimension(nx) :: temp_Br, temp_Bphi
            character(len=30) :: ghost_zone_type = 'anti-symmetric'

            call impose_boundary_conditions(B_r, ghost_zone_type)
            call impose_boundary_conditions(B_phi, ghost_zone_type)
            call impose_boundary_conditions(old_Br, ghost_zone_type)
            call impose_boundary_conditions(old_Bphi, ghost_zone_type)
            call construct_grid
            alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
            alpha_Bphi = alpha_cap2*B_phi
            alpha_Br = alpha_cap2*B_r
            Uz_Br = B_r*U_z_cap
            Uz_Bphi = B_phi*U_z_cap


            do i = nxghost+1,nxphys+nxghost
                temp_Br(i) = (1+2*alp)*((beta)*(alpha_Bphi(i+1)-alpha_Bphi(i-1)) &
                +2*alp*(B_r(i+1)+B_r(i-1)) - (beta)*(Uz_Br(i+1) - Uz_Br(i-1)))+old_Br(i)

                temp_Bphi(i) = (1+2*alp+G)*((beta)*(alpha_Br(i+1)-alpha_Br(i-1)) &
                +2*alp*(B_phi(i+1)+B_phi(i-1)) - (beta)*(Uz_Bphi(i+1) - Uz_Bphi(i-1)))+old_Bphi(i) !NOT_SURE
                

            end do
            B_r = temp_Br
            B_phi = temp_Bphi
            t = t + dt

        end subroutine central_difference


        subroutine RK4_new
            character(len=30) :: ghost_zone_type = 'anti-symmetric'
            call impose_boundary_conditions(B_r, ghost_zone_type)
            call impose_boundary_conditions(B_phi, ghost_zone_type)

            alpha_cap2 = alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
            alpha_Bphi = alpha_cap2*B_phi
            alpha_Br = alpha_cap2*B_r
            Uz_Br = B_r*U_z_cap
            Uz_Bphi = B_phi*U_z_cap

            call spatial_derivative(B_r,2,dBr,d2Br)
            call spatial_derivative(B_phi,2,dBphi,d2Bphi)
            call spatial_derivative(alpha_Br,2,d_alpha_Br,d2_alpha_Br)
            call spatial_derivative(alpha_Bphi,2,d_alpha_Bphi,d2_alpha_Bphi)
            call spatial_derivative(Uz_Br,2,d_Uz_Br,d2_Uz_Br)
            call spatial_derivative(Uz_Bphi,2,d_Uz_Bphi,d2_Uz_Bphi)




            call diff_equations_no_split(B_r,B_phi,dBr,dBphi,d2Br,d2Bphi,d_alpha_Bphi, &
            d_alpha_Br,d_Uz_Bphi,d_Uz_Br)
            k1r = dt*dBrdt
            k1phi = dt*dBphidt

            call diff_equations_no_split(B_r+0.5*k1r,B_phi+0.5*k1phi,dBr,dBphi,d2Br,d2Bphi, &
            d_alpha_Bphi,d_alpha_Br,d_Uz_Bphi,d_Uz_Br)
            k2r = dt*dBrdt
            k2phi = dt*dBphidt

            call diff_equations_no_split(B_r+0.5*k2r,B_phi+0.5*k2phi,dBr,dBphi,d2Br, &
            d2Bphi,d_alpha_Bphi,d_alpha_Br,d_Uz_Bphi,d_Uz_Br)
            k3r = dt*dBrdt
            k3phi = dt*dBphidt

            call diff_equations_no_split(B_r+k3r,B_phi+k3phi,dBr,dBphi,d2Br,d2Bphi, &
            d_alpha_Bphi,d_alpha_Br,d_Uz_Bphi,d_Uz_Br)
            k4r = dt*dBrdt
            k4phi = dt*dBphidt

            B_r = B_r + (k1r + 2.*k2r + 2.*k3r + k4r)/6.
            B_phi = B_phi + (k1phi + 2.*k2phi + 2.*k3phi + k4phi)/6.

            t = t + dt



        end subroutine RK4_new

        subroutine RK3_implicit

            implicit none
            character(len=30) :: ghost_zone_type = 'anti-symmetric'
            g_1 = 8.0 / 15.0        !
            g_2 = 5.0 / 12.0        !
            g_3 = 3.0 / 4.0         !----> REFER: Brandenburg, A. (2001).ArXiv. 
            z_1 = -17.0 / 60.0      !             Computational aspects of astrophysical MHD and turbulence.  
            z_2 = -5.0 / 12.0       !             https://doi.org/10.1201/9780203493137.ch9
            
         
            call impose_boundary_conditions(B_r, ghost_zone_type)
            call impose_boundary_conditions(B_phi, ghost_zone_type)

           
            alpha_cap2 = alpha_cap / (1.0 + (B_r**2 + B_phi**2) / B_eq**2)
            alpha_Bphi = alpha_cap2 * B_phi
            alpha_Br = alpha_cap2 * B_r
            Uz_Br = B_r * U_z_cap
            Uz_Bphi = B_phi * U_z_cap

  
            call spatial_derivative(B_r, 2, dBr, d2Br)
            call spatial_derivative(B_phi, 2, dBphi, d2Bphi)
            call spatial_derivative(alpha_Br, 2, d_alpha_Br, d2_alpha_Br)
            call spatial_derivative(alpha_Bphi, 2, d_alpha_Bphi, d2_alpha_Bphi)
            call spatial_derivative(Uz_Br, 2, d_Uz_Br, d2_Uz_Br)
            call spatial_derivative(Uz_Bphi, 2, d_Uz_Bphi, d2_Uz_Bphi)


  

            
            call diff_equations_no_split(B_r, B_phi, dBr, dBphi, d2Br, d2Bphi, d_alpha_Bphi, d_alpha_Br, d_Uz_Bphi, d_Uz_Br)
            k1r = dt * dBrdt
            k1phi = dt * dBphidt
            temp_Br = B_r + g_1 * k1r
            temp_Bphi = B_phi + g_1 * k1phi

            B_r = temp_Br + z_1 * k1r
            B_phi = temp_Bphi + z_1 * k1phi
            t = t + g_1 * dt

            
            call diff_equations_no_split(B_r, B_phi, dBr, dBphi, d2Br, d2Bphi, d_alpha_Bphi, d_alpha_Br, d_Uz_Bphi, d_Uz_Br)
            k2r = dt * dBrdt
            k2phi = dt * dBphidt
            temp_Br = B_r + g_2 * k2r
            temp_Bphi = B_phi + g_2 * k2phi

            B_r = temp_Br + z_2 * k2r
            B_phi = temp_Bphi + z_2 * k2phi
            t = t + g_2 * dt

            ! 3rd step
            call diff_equations_no_split(B_r, B_phi, dBr, dBphi, d2Br, d2Bphi, d_alpha_Bphi, d_alpha_Br, d_Uz_Bphi, d_Uz_Br)
            k3r = dt * dBrdt
            k3phi = dt * dBphidt
            temp_Br = B_r + g_3 * k3r
            temp_Bphi = B_phi + g_3 * k3phi

            B_r = temp_Br
            B_phi = temp_Bphi
            t = t + g_3 * dt



        end subroutine RK3_implicit


         









    


end module timestepping