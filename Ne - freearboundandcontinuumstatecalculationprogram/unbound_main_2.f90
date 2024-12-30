program calculate_energy
    implicit none

    integer, parameter :: n = 2000  ! Number of grid points
    integer :: i
    real*8 :: r(n), u(n), up(n)
    real*8 :: l, Z, q, energy
    real*8 :: energy_guess

    ! Set up parameters for the 3s orbital of Argon
    l = 1.d0          ! Orbital angular momentum quantum number for s-orbital
    Z = 18.d0         ! Atomic number of Argon
    q = 1.d0          ! Potential parameter (set as needed)
    
    ! Initialize radial grid
    do i = 1, n
        r(i) = i * 0.1d0  ! Example radial grid points, adjust as needed
    end do

    ! Guess for the energy (for example, -1.0 atomic units, adjust as needed)
    energy_guess = -1.0d0

    ! Call the subroutine to calculate the energy
    call unbound(n, u, up, r, l, Z, q, energy_guess)


    ! Output the result
    print*, 'Calculated energy for the unbound 3s orbital of Argon: ', energy_guess
    do i = 1,n
      print*, r(i), u(i)
    end do

end program calculate_energy

