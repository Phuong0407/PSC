! This is the program use to solve the equation
!   k(u_xx + u_yy + u_zz) + Q = 0
! On Dirichlet Boundary
!   u = u_D
! On Newman Boundary
!

program axisym_quad_stiff_local_test
    use axisym_quad_stiff_loc
    implicit none

    integer :: j
    real, allocatable :: axisym_K_loc(:,:)
    real, allocatable :: coord_mat(:,:)

    allocate(axisym_K_loc(4, 4))
    allocate(coord_mat(4, 2))

    axisym_K_loc = cal_axisym_quad_stiff_loc(0.0, 1.0, 2.0, 1.0, 2.0, 2.0, 0.0, 2.0)

    print *, "Axisymmetric Local Stiffness Matrix:"
    do j = 1, 4
        write(*, '(4F10.4, A)') axisym_K_loc(j, :), "  "
    end do
end program axisym_quad_stiff_local_test
