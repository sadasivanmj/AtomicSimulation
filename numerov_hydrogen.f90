PROGRAM NumerovShrodinger
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 1000
  INTEGER :: i
  REAL, PARAMETER :: rmin = 1.0e-5
  REAL, PARAMETER :: rmax = 20.0
  REAL, PARAMETER :: dr = (rmax - rmin) / (n - 1)
  REAL, PARAMETER :: hbar = 1.0545718e-34
  REAL, PARAMETER :: m_e = 9.10938356e-31
  REAL, PARAMETER :: q = 1.602176634e-19
  REAL, PARAMETER :: eps0 = 8.8541878173e-12
  REAL, PARAMETER :: Z = 1.0
  REAL, PARAMETER :: E = -13.6 * q  ! Converting eV to Joules
  INTEGER, PARAMETER :: l = 0
  REAL, DIMENSION(n) :: r, u, V, k

  ! Initializing the radial grid
  r = [(rmin + i * dr, i = 0, n-1)]

  ! Potential energy V(r)
  V = -Z * q**2 / (4 * 3.141592653589793 * eps0 * r)

  ! Initialize wavefunction arrays
  u = 0.0
  u(1) = 0.0
  u(2) = dr**(l + 1)  ! Small value for the starting point

  ! Calculate k for all points
  k = 2 * m_e / hbar**2 * (E - V - l * (l + 1) * hbar**2 / (2 * m_e * r**2))

  ! Numerov integration loop
  DO i = 2, n-1
    u(i+1) = (2 * (1 - 5/12 * dr**2 * k(i)) * u(i) - (1 + dr**2 / 12 * k(i-1)) * u(i-1)) / (1 + dr**2 / 12 * k(i+1))
  END DO

  ! Normalize the wavefunction
  u = u / SQRT(SUM(u**2 * dr))

  ! Print the results
  PRINT *, "Radial Wavefunction:"
  DO i = 1, n
    PRINT *, r(i), u(i)
  END DO

END PROGRAM NumerovShrodinger

