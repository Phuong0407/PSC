program streamline_function_inlet_solution
    use streamline_solution_triangular
    use csr_sparse_matrix
    implicit none
    ! real :: v_x
    integer :: unit_number
    integer :: i
    call INIT_SOLUTION(3.0D0, 3.0D0, 3.0D0, 4.0D0, 2.0D0, 200, 5, 200, 5, 1.0D0)
    call INIT_SPARSE_SOLVER(STIFF_MAT_GLO)

    NODAL_VALS = SOLVE_SPARSE_SYSTEM(FORCE_MAT_GLO)

    ! do i = 1, size(GRID%NODE_CONN, 1)
    !     print *, GRID%NODE_CONN(i, :)
    ! end do

    do i = 1, size(NODAL_VALS)
        print *, i, "   ", NODAL_VALS(i)
    end do


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
    open(unit=unit_number, file="output.txt", status="unknown", action="write")

    do i = 1, size(NODE_GRAD_RECOVERY, 1)
        WRITE(unit_number, '(A,F24.16,A,F24.16)') "RECOVERY DERIVATIVE AT NODE (", GRID%COORD(i, 1), ",", GRID%COORD(i, 2), ")"
        WRITE(unit_number, '(A,F24.16)') "   Psi_x = ", NODE_GRAD_RECOVERY(i, 1)
        WRITE(unit_number, '(A,F24.16)') "   Psi_y = ", NODE_GRAD_RECOVERY(i, 2)
        WRITE(unit_number, '(A,F24.16)') "   Psi_xy = ", NODE_GRAD_RECOVERY(i, 3)
        WRITE(unit_number, '(A,F24.16)') "   Psi_xx = ", NODE_GRAD_RECOVERY(i, 4)
        WRITE(unit_number, '(A,F24.16)') "   Psi_yy = ", NODE_GRAD_RECOVERY(i, 5)
    end do    
    ! Close the file
    close(unit_number)
end program streamline_function_inlet_solution