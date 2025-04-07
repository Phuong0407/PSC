program Laplace
    use quadrilateral
    implicit none

    ! Declare variables for node coordinates
    real :: x1, y1, x2, y2, x3, y3, x4, y4
    real, dimension(4, 4) :: stiffness_matrix

    ! Initialize coordinates of the quadrilateral element
    x1 = 1.0; y1 = 0.0
    x2 = 0.0; y2 = 0.0
    x3 = 2.0; y3 = 0.5
    x4 = 2.0; y4 = 1.0

    ! Call the cal_stiff_loc function to compute the stiffness matrix
    stiffness_matrix = cal_stiff_loc(x1, y1, x2, y2, x3, y3, x4, y4)

    ! Print the resulting stiffness matrix
    print *, "Stiffness Matrix:"
    print *, stiffness_matrix
end program Laplace