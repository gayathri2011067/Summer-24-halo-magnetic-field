!************************************************************************************************************************************************************
module courant
         use para
  !use physical_params
!
  implicit none
!
  integer, parameter :: Nzghost= 1, Nz=n1+2*Nzghost !Number of ghost cells at each end in z ! Nzphys=101
  integer :: k
  double precision :: length,dz,dt,c,cc,D !Simulation domain and grid spacing (in units of L_UNIT)
  double precision, dimension(Nz) ::z
!
contains
  subroutine courant_cond 
    length=    2.*h
    dz=     length/(n1-1)/h  !x corresponds to z coordinate
    D=R_alp*R_w       
    dt=tf/(Nt-1)!/t_d
    c=dt/dz
    cc=(dt/(dz**2))!*cc
  endsubroutine courant_cond 
end module courant
!**************************************************************************************************************************************************************
