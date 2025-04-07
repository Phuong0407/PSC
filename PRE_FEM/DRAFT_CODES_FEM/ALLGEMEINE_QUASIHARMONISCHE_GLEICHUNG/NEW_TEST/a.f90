MODULE FlowThroughCube
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    INTEGER, PARAMETER :: Nx = 50, Ny = 50, Nz = 50
    INTEGER, PARAMETER :: maxSteps = 10000
    REAL(dp), PARAMETER :: Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0
    REAL(dp), PARAMETER :: dt = 1.0D-4, delta = 1.0D0
    REAL(dp), PARAMETER :: nu = 0.01D0
    REAL(dp), DIMENSION(0:Nx+1, 0:Ny+1, 0:Nz+1) :: &
         u1, u2, u3, rho, u1_new, u2_new, u3_new, rho_new
  
  CONTAINS
  
    ! Subroutine to initialize fields for flow through a cube
    SUBROUTINE InitializeFields()
      INTEGER :: i, j, k
      DO k = 0, Nz+1
        DO j = 0, Ny+1
          DO i = 0, Nx+1
            ! Initial velocity field
            u1(i,j,k) = 0.0D0
            u2(i,j,k) = 0.0D0
            u3(i,j,k) = 0.0D0
            ! Initial density field
            rho(i,j,k) = 1.0D0
          END DO
        END DO
      END DO
  
      ! Set inlet velocity profile at x=0
      DO k = 1, Nz
        DO j = 1, Ny
          u1(0,j,k) = 1.0D0  ! Uniform inlet velocity
          u2(0,j,k) = 0.0D0
          u3(0,j,k) = 0.0D0
          rho(0,j,k) = 1.0D0
        END DO
      END DO
    END SUBROUTINE InitializeFields
  
    ! Subroutine to apply boundary conditions for the cube
    SUBROUTINE ApplyBoundaryConditions()
      INTEGER :: i, j, k
  
      ! Inlet (x = 0)
      DO k = 1, Nz
        DO j = 1, Ny
          u1(0,j,k) = 1.0D0
          u2(0,j,k) = 0.0D0
          u3(0,j,k) = 0.0D0
          rho(0,j,k) = 1.0D0
        END DO
      END DO
  
      ! Outlet (x = Lx)
      DO k = 1, Nz
        DO j = 1, Ny
          u1(Nx+1,j,k) = u1(Nx,j,k)  ! Extrapolated velocity
          u2(Nx+1,j,k) = u2(Nx,j,k)
          u3(Nx+1,j,k) = u3(Nx,j,k)
          rho(Nx+1,j,k) = rho(Nx,j,k)
        END DO
      END DO
  
      ! No-slip walls (y = 0, Ly and z = 0, Lz)
      DO i = 0, Nx+1
        DO k = 0, Nz+1
          u1(i,0,k) = 0.0D0
          u2(i,0,k) = 0.0D0
          u3(i,0,k) = 0.0D0
          u1(i,Ny+1,k) = 0.0D0
          u2(i,Ny+1,k) = 0.0D0
          u3(i,Ny+1,k) = 0.0D0
        END DO
  
        DO j = 0, Ny+1
          u1(i,j,0) = 0.0D0
          u2(i,j,0) = 0.0D0
          u3(i,j,0) = 0.0D0
          u1(i,j,Nz+1) = 0.0D0
          u2(i,j,Nz+1) = 0.0D0
          u3(i,j,Nz+1) = 0.0D0
        END DO
      END DO
    END SUBROUTINE ApplyBoundaryConditions
  
    ! Main simulation loop
    SUBROUTINE SimulateFlow()
      INTEGER :: nStep
      CALL InitializeFields()
  
      DO nStep = 1, maxSteps
        CALL ApplyBoundaryConditions()
        CALL UpdateVelocities()
        CALL UpdateDensity()
  
        ! Swap new and old fields
        u1 = u1_new
        u2 = u2_new
        u3 = u3_new
        rho = rho_new
  
        ! Output results at regular intervals
        IF (MOD(nStep, 100) == 0) THEN
          PRINT *, 'Step:', nStep
        END IF
      END DO
    END SUBROUTINE SimulateFlow
  
  END MODULE FlowThroughCube
  