program streamline_function_inlet_solution
    use streamline_solution_triangular
    use csr_sparse_matrix
    implicit none
    ! real :: v_x
    integer :: unit_number
    integer :: i
    real(8) :: VELOCITY(2)

    call INIT_SOLUTION(3.0D0, 3.0D0, 3.0D0, 4.0D0, 2.0D0, 2, 2, 2, 2, 1.0D0)
    
    ! Assign a unit number to the file
    unit_number = 10

    ! Open the file for writing
    OPEN(unit=unit_number, file="mod_global_stiff_mat.txt", status='replace', action='write', iostat=i)
    IF (i /= 0) THEN
        PRINT *, "Error opening file!"
        STOP
    END IF

    ! Write the specific matrix line by line
    DO i = 1, SIZE(STIFF_MAT_GLO, 1)
        WRITE(unit_number, '(21F8.2)') STIFF_MAT_GLO(i, :)
    END DO

    ! Close the file
    CLOSE(unit_number)

    call INIT_SPARSE_SOLVER(STIFF_MAT_GLO)

    NODAL_VALS = SOLVE_SPARSE_SYSTEM(FORCE_MAT_GLO)

    ! do i = 1, size(GRID%NODE_CONN, 1)
    !     print *, GRID%NODE_CONN(i, :)
    ! end do

    ! do i = 1, size(NODAL_VALS)
    !     print *, i,  "   ", NODAL_VALS(i)
    ! end do

    ! do i = 1, size(NODAL_VALS)
    !     print *, i,  "   ", FORCE_MAT_GLO(i)
    ! end do

    ! v_x = CALC_DEV_Y(1.6, 4.0)
    ! print* , v_x
    ! v_x = CALC_DEV_Y(1.7, 4.0)
    ! print* , v_x
    ! v_x = CALC_DEV_Y(1.8, 4.0)
    ! print* , v_x
    ! v_x = CALC_DEV_Y(1.9, 4.0)
    ! print* , v_x
    ! v_x = CALC_DEV_Y(1.6, 3.5)
    ! print* , v_x

    ! do i = 1, size(NODAL_VALS)
    !     print *, NODAL_VALS(i)
    ! end do

    call CALC_NODE_GRAD_RECOVERY()
    ! print*, SIZE(NODE_GRAD_RECOVERY,1)
    ! do i = 1, SIZE(NODE_GRAD_RECOVERY,1)
    !     ! print *, "RECOVERY DERIVATIVE AT NODE", i, "   Psi_xx", NODE_GRAD_RECOVERY(i, 4), "   Psi_yy", NODE_GRAD_RECOVERY(i, 5)
    !     print *, "RECOVERY DERIVATIVE AT NODE", i, "Psi_x = ", NODE_GRAD_RECOVERY(i, 1), "Psi_y = ", NODE_GRAD_RECOVERY(i, 2)
    ! end do

    unit_number = 10  ! You can choose any valid unit number
    open(unit=unit_number, file="gradient_out.txt", status="unknown", action="write")

    do i = 1, size(NODE_GRAD_RECOVERY, 1)
        ! Properly formatted node coordinates with compact spacing
        WRITE(unit_number, '(A,F24.16,A,F24.16,A)') "RECOVERY DERIVATIVE AT NODE (", GRID%COORD(i, 1), ", ", GRID%COORD(i, 2), ")"
        ! Use exactly 4 spaces before each Psi entry
        WRITE(unit_number, '(4X,A,F24.16)') "Psi_z =", NODE_GRAD_RECOVERY(i, 1)
        WRITE(unit_number, '(4X,A,F24.16)') "Psi_r =", NODE_GRAD_RECOVERY(i, 2)
        WRITE(unit_number, '(4X,A,F24.16)') "Psi_zz =", NODE_GRAD_RECOVERY(i, 3)
        WRITE(unit_number, '(4X,A,F24.16)') "Psi_zr =", NODE_GRAD_RECOVERY(i, 4)
        WRITE(unit_number, '(4X,A,F24.16)') "Psi_rr =", NODE_GRAD_RECOVERY(i, 5)
    end do    

    ! Close the file
    close(unit_number)

    unit_number = 10  ! You can choose any valid unit number
    open(unit=unit_number, file="psi_out.txt", status="unknown", action="write")

    do i = 1, size(NODE_GRAD_RECOVERY, 1)
        ! Properly formatted node coordinates with compact spacing
        WRITE(unit_number, '(A,F24.16,A,F24.16,A)') "STREAMFUNCTION AT(", GRID%COORD(i, 1), ", ", GRID%COORD(i, 2), ")"
        ! Use exactly 4 spaces before Psi
        WRITE(unit_number, '(4X,A,F24.16)') "Psi =", NODAL_VALS(i)
    end do
    

    ! Close the file
    close(unit_number)
    call CALCULATE_VELOCITY_AT_POINT(0.1D0,0.1D0, VELOCITY)
    print*, "VELOCITY AT POINT (0.1, 0.1) = ", VELOCITY(1), VELOCITY(2)

    call CALCULATE_VELOCITY_AT_POINT(0.5D0,0.0D0, VELOCITY)
    print*, "VELOCITY AT POINT (0.5, 0.0) = ", VELOCITY(1), VELOCITY(2)

    call CALCULATE_VELOCITY_AT_POINT(0.0D0,0.0D0, VELOCITY)
    print*, "VELOCITY AT POINT (0.0, 0.0) = ", VELOCITY(1), VELOCITY(2)

    call CALCULATE_VELOCITY_AT_POINT(0.0D0,2.0D0, VELOCITY)
    print*, "VELOCITY AT POINT (0.0, 2.0) = ", VELOCITY(1), VELOCITY(2)

    call CALCULATE_VELOCITY_AT_POINT(3.6D0,0.4D0, VELOCITY)
    print*, "VELOCITY AT POINT (3.6, 0.4) = ", VELOCITY(1), VELOCITY(2)
end program streamline_function_inlet_solution