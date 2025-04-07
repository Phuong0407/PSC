program streamline_function_inlet_solution
    use streamline_solution_triangular
    use csr_sparse_matrix
    implicit none
    real :: result
    integer :: i
    ! integer :: j
    ! integer :: unit
    ! character(len=100) :: filename
    
    call INIT_SOLUTION(3.0, 3.0, 3.0, 4.0, 2.0, 2, 2, 2, 2, 1.0)
    call INIT_SPARSE_SOLVER(STIFF_MAT_GLO)
    VALS = SOLVE_SPARSE_SYSTEM(FORCE_MAT_GLO)
    ! do i = 1, GRID%n_Tri_Elem
    !     print *, "ELEMENT", i
    !     print *, GRID%COORD(GRID%CONN(i, 1), 1), "  ", GRID%COORD(GRID%CONN(i, 1), 2)
    !     print *, GRID%COORD(GRID%CONN(i, 2), 1), "  ", GRID%COORD(GRID%CONN(i, 2), 2)
    !     print *, GRID%COORD(GRID%CONN(i, 3), 1), "  ", GRID%COORD(GRID%CONN(i, 3), 2)
    ! end do

    ! print *, INTERPOLATE_STREAMLINE_FUNC(3.0,3.0)
    ! print *, INTERPOLATE_STREAMLINE_FUNC(2.8,1.2)
    ! print *, "Solution vector:"
    ! do i = 1, size(VALS)
    !     print *, "PSI(", GRID%COORD(i, 1), ",", GRID%COORD(i, 2), ") = ", VALS(i)
    ! end do

    ! do i =1, size(GRID%COORD,1)
    !     print *, GRID%COORD(i, 1), ",", GRID%COORD(i, 2), 
    ! end do

    ! print *, size(GRID%CONN, 1)
    ! ! print *, GRID%n_Tri_Elem
    ! do i =1, size(GRID%CONN, 1)
    !     print *, "[", GRID%CONN(i, 1), ",", GRID%CONN(i, 2), ",", GRID%CONN(i, 3), "]"
    ! end do

    ! do i = 1, size(FORCE_MAT_GLO)
    !     print *, VALS(i)
    ! end do


    ! do i = 1, size(GRID%CONN, 1)
    !     print *, (GRID%CONN(i, j), j = 1, size(GRID%CONN, 2))
    ! end do
    result = INTERPOLATE_STREAMLINE_FUNC(0.2,0.3)
    print *, "Result: ", result




    ! print *, INTERPOLATE_STREAMLINE_FUNC(0.2, 0.1)

    ! print *, (INTERPOLATE_STREAMLINE_FUNC(2.0, 2.5) - INTERPOLATE_STREAMLINE_FUNC(2.0, 2.49999))/(2.5-2.49999)
    ! print *, (INTERPOLATE_STREAMLINE_FUNC(2.0, 2.5) - INTERPOLATE_STREAMLINE_FUNC(2.0, 2.4999999))/(2.5-2.4999999)
    ! print *, (INTERPOLATE_STREAMLINE_FUNC(2.0, 2.5) - INTERPOLATE_STREAMLINE_FUNC(2.0, 2.499))/(2.5-2.499)



    ! filename = "streamline_values.txt"
    ! unit = 10  ! Logical unit number

    ! ! Open the file
    ! open(unit=unit, file=filename, status="replace", action="write", &
    !      position="rewind", form="formatted")

    ! ! Write data to the file
    ! do i = 1, GRID%nNode
    !     ! write(unit, '(A,F0.3,A,F0.3,A,F0.3)') "PSI(", GRID%COORD(i, 1), ",", GRID%COORD(i, 2), ")=", VALS(i)
    !     write(unit, '(A,F10.5,A,F10.5,A,F15.5)') "PSI(", GRID%COORD(i, 1), ",", GRID%COORD(i, 2), ")=", VALS(i)
    ! end do

    ! ! Close the file
    ! close(unit)

    
    ! filename = "grid_configuration_MATLAB.txt"
    ! unit = 10  ! Logical unit number

    ! ! Open the file
    ! open(unit=unit, file=filename, status="replace", action="write", &
    !      position="rewind", form="formatted")

    ! ! Write data to the file
    ! do i =1, size(VALS)
    !     print *, GRID%COORD(i, 1), ",", GRID%COORD(i, 2), ";"
    ! end do

    ! do i =1, size(GRID%CONN)
    !     print *, GRID%CONN(i, 1), ",", GRID%CONN(i, 2), ",", GRID%CONN(i, 3), ";"
    ! end do

    ! ! Close the file
    ! close(unit)

end program streamline_function_inlet_solution