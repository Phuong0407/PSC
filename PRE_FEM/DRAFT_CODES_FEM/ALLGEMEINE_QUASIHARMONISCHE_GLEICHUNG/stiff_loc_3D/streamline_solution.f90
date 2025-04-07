module streamline_solution
    use axisym_stiff_glo
    implicit none
    private :: init_grid, gen_grid_coord, gen_grid_con!, gen_bound_mat!, gen_bound_ind
    public :: grid_struct, grid_dat, init_solution_data, stiff_glo

    type :: grid_struct
        integer :: N                                        ! Total number of vertical nodes
        integer :: M                                        ! Total number of horizontal nodes
        real, allocatable :: x(:,:), y(:,:)                 ! Coordinate matrices
        integer, allocatable :: node_con_mat(:,:)           ! Connection Data for Matrix Assemble
        ! integer, allocatable :: lower_bound_node_ind(:,:)   ! Index of Lower Boundary Node, in this project, it is all about Dirichlet Boundary
        ! integer, allocatable :: upper_bound_node_ind(:,:)   ! Index of Upper Boundary Node, in this project, it is all about Dirichlet Boundary
        ! integer, allocatable :: left_bound_node_ind(:,:)    ! Index of Left Boundary Node, in this project, it is all about Dirichlet Boundary
        ! integer, allocatable :: right_bound_node_ind(:,:)   ! Index of Right Boundary Node, in this project, it is all about Dirichlet Boundary
        ! Add Neumann boundary condition and mixed condition in the future
    contains
    end type grid_struct

    real :: H_1, H_2, H_3           ! Horizontal dimension of each region
    real :: dx_1, dx_2, dx_3        ! Finite-Diffrent horizontal dimension of each region
    real :: L_1, L_2                ! Vertical dimension of each region
    real :: dy_1, dy_2              ! Finite-Diffrent vertical dimension of each region
    integer :: N                    ! Vertical number of element
    integer :: M_1, M_2, M_3        ! Horizontal number of element for each region
    real :: U_0                     ! Velocity of the Entrance of Engine
    
    type(grid_struct) :: grid_dat
    real, allocatable :: force_mat(:,:)
    real, allocatable :: stiff_glo(:,:)
    real, allocatable :: nodal_val(:,:)

contains
    subroutine init_grid(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3)
        ! class(grid_struct), intent(inout) :: this
        real, intent(in) :: iH_1, iH_2, iH_3, iL_1, iL_2
        integer, intent(in) :: iN, iM_1, iM_2, iM_3

        H_1 = iH_1
        H_2 = iH_2
        H_3 = iH_3
        L_1 = iL_1
        L_2 = iL_2
        dx_1 = H_1 / M_1
        dx_2 = H_2 / M_2
        dx_3 = H_3 / M_3
        dy_1 = L_1 / N
        dy_2 = L_2 / N
        N = iN
        M_1 = iM_1
        M_2 = iM_2
        M_3 = iM_3

        grid_dat%N = N + 1
        grid_dat%M = M_1 + M_2 + M_3 + 1

        allocate(grid_dat%x(N + 1, M_1 + M_2 + M_3 + 1))
        allocate(grid_dat%y(N + 1, M_1 + M_2 + M_3 + 1))
    end subroutine init_grid

    subroutine gen_grid_coord()
        ! class(grid_struct), intent(inout) :: this
        integer :: i, j
        real :: x_2

        do i = 0, N
            ! Region 1: First Rectangular Region
            do j = 0, M_1
                grid_dat%x(i+1, j+1) = H_1 / M_1 * j
                grid_dat%y(i+1, j+1) = L_1 / N * i
            end do

            ! Region 2: Oblique Region
            do j = M_1 + 1, M_1 + M_2
                x_2 = H_1 + (H_2 / M_2) * (j - M_1)
                grid_dat%x(i+1, j+1) = x_2
                grid_dat%y(i+1, j+1) = ((L_2 - L_1) / (H_2 * N)) * i * x_2 + &
                                   (L_1 / N) * i - ((L_2 - L_1) / (H_2 * N)) * i * H_1
            end do

            ! Region 3: Second Rectangular Region
            do j = M_1 + M_2 + 1, M_1 + M_2 + M_3
                grid_dat%x(i+1, j+1) = H_1 + H_2 + (H_3 / M_3) * (j - M_1 - M_2)
                grid_dat%y(i+1, j+1) = L_2 / N * i
            end do
        end do
    end subroutine gen_grid_coord

    subroutine gen_grid_con()
        integer :: i, j, ind_elem

        allocate(grid_dat%node_con_mat((grid_dat%N - 1) * (grid_dat%M - 1), 4))

        do i = 1, grid_dat%N - 1
            do j = 1, grid_dat%M - 1
                ind_elem = (grid_dat%M - 1) * (i - 1) + j
                print *, ind_elem
                grid_dat%node_con_mat(ind_elem, 1) = grid_dat%M * i + j - grid_dat%M
                grid_dat%node_con_mat(ind_elem, 2) = grid_dat%M * i + j - grid_dat%M + 1
                grid_dat%node_con_mat(ind_elem, 3) = grid_dat%M * i + j + 1
                grid_dat%node_con_mat(ind_elem, 4) = grid_dat%M * i + j
            end do
        end do
    end subroutine gen_grid_con

    ! subroutine gen_bound_mat()
    !     integer :: i
    !     integer :: l_bound_ind, r_bound_ind
    !     real :: r_1, r_2
    !     real :: U_1

    !     U_1 = U_0 * (L_1/(2.0*L_1 - L_2)) * (L_1/L_2)
    !     allocate(force_mat(grid_dat%M * grid_dat%N, 1))
    !     force_mat = 0.0

    !     do i = 1, grid_dat%M
    !         force_mat(1, i) = 0                                                 ! Lower Boundary Value
    !         force_mat(grid_dat%N, i) = 0.5 * U_0 * L_1**2.0                     ! Upper Boundary Value
    !     end do
    !     do i = 2, grid_dat%N - 1
    !         l_bound_ind = i * grid_dat%M + 1
    !         r_bound_ind = i * grid_dat%M + grid_dat%M
    !         r_1 = grid_dat%y(l_bound_ind, 1)
    !         r_2 = grid_dat%y(l_bound_ind, grid_dat%M)
    !         force_mat(l_bound_ind, 1) = 0.5 * U_0 * r_1**2.0                    ! Left Boundary Value
    !         force_mat(r_bound_ind, 1) = 0.5 * U_1 * r_2**2.0                    ! Right Boundary Value
    !     end do
    ! end subroutine gen_bound_mat
    ! subroutine gen_bound_ind()
    !     integer :: i, j
    !     allocate(grid_dat%lower_bound_node_ind(grid_dat%M, 1))
    !     allocate(grid_dat%upper_bound_node_ind(grid_dat%M, 1))
    !     allocate(grid_dat%left_bound_node_ind(grid_dat%M, 1))
    !     allocate(grid_dat%right_bound_node_ind(grid_dat%M, 1))
        
    !     do i = 1, grid_dat%M
    !         lower_bound_node_ind(i,
    !     end do
    ! end subroutine gen_bound_ind

    subroutine init_solution_data(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3)
        integer :: i, j
        real, intent(in) :: iH_1, iH_2, iH_3, iL_1, iL_2
        integer, intent(in) :: iN, iM_1, iM_2, iM_3
        real, allocatable :: grid_coord(:,:)

        call init_grid(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3)
        call gen_grid_coord()
        call gen_grid_con()

        allocate(nodal_val(grid_dat%M * grid_dat%N, 1))
        allocate(stiff_glo(grid_dat%M * grid_dat%N, grid_dat%M * grid_dat%N))
        allocate(grid_coord(grid_dat%N * grid_dat%M, 2))
        
        do i = 1, grid_dat%N
            do j = 1, grid_dat%M
                grid_coord(i * grid_dat%M + j, 1) = grid_dat%x(i, j)
                grid_coord(i * grid_dat%M + j, 2) = grid_dat%y(i, j)
            end do
        end do

        call init_gsm(grid_dat%N, grid_dat%M, grid_dat%N - 1, grid_dat%M - 1, grid_coord, grid_dat%node_con_mat)
    end subroutine init_solution_data
end module streamline_solution