- read the review paper (half)
- Read through 1993 paper (plotted profiles)

- didnt do python this time

***fortran***

- parameters (from sirs code)
- profiles (from paper)
  - currently simplified to avoid quenching
  - still matches tho

- generalized derviatives (finite diff)
  - multiple ghost zones option
  - multiple order of fd options
  - need to integrate fd order with ghost zone
- generalized time_stepping routine
  - implemented rk4
  - boundary condition maintained (not for alpha since no quenching)
- multiple files for readability
- bash file created and run simplified