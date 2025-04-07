module axisym_stiff_glo
    use axisym_quad_elem
    implicit none
    private
    integer :: n_node, nv_node, nh_node, n_elem, nh_elem, nv_elem
    real, parameter :: pi = 3.141592653589793

    public :: init_gsm, assemble_gsm
    real, allocatable :: grid_coord(:,:)
    real, allocatable :: stiff_glo_mat(:,:)
    ! real, allocatable :: heat_gen_glo_mat(:,:)                                ! TO BE DONE IN THE FUTURE
    integer, allocatable :: node_con_mat(:,:)

contains
    subroutine init_gsm(nv_node_in, nh_node_in, nv_elem_in, nh_elem_in, grid_coord_in, node_con_in)
        implicit none
        integer, intent(in) :: nh_node_in, nv_node_in
        integer, intent(in) :: nh_elem_in, nv_elem_in
        real, intent(in) :: grid_coord_in(:,:)
        integer, intent(in) :: node_con_in(:,:)

        allocate(grid_coord(n_node, 2))
        allocate(node_con_mat(n_node, 4))

        if (nh_node_in /= nh_elem_in + 1 .or. nv_node_in /= nv_elem_in + 1) then
            stop
        endif

        nh_node = nh_node_in
        nv_node = nv_node_in
        nh_elem = nh_elem_in
        nv_elem = nv_elem_in
        n_node = nh_node * nv_node
        n_elem = nh_elem_in * nv_elem_in

        allocate(grid_coord(n_node, 2))
        allocate(node_con_mat(n_elem, 4))
        allocate(stiff_glo_mat(n_node, n_node))
        ! allocate(heat_gen_glo_mat(n_node, 1))                                 ! TO BE DONE IN THE FUTURE

        stiff_glo_mat = 0.0
        ! heat_gen_glo_mat = 0.0                                                ! TO BE DONE IN THE FUTURE
        grid_coord = grid_coord_in
        node_con_mat = node_con_in
    end subroutine init_gsm

    subroutine assemble_gsm()
        implicit none
        integer :: node_1, node_2, node_3, node_4
        integer :: i, j, elem
        real :: x1, y1, x2, y2, x3, y3, x4, y4

        do elem = 1, n_elem
            ! TO BE DONE IN THE FUTURE
            ! if (size(node_con_mat, 2) == 4) then
            !     ! USE 4-node elements
            ! else if (size(node_con_mat, 2) == 3) then
            !     ! USE 3-node elements
            ! end if

            node_1 = node_con_mat(elem, 1)
            node_2 = node_con_mat(elem, 2)
            node_3 = node_con_mat(elem, 3)
            node_4 = node_con_mat(elem, 4)
            x1 = grid_coord(node_1, 1)
            y1 = grid_coord(node_1, 2)
            x2 = grid_coord(node_2, 1)
            y2 = grid_coord(node_2, 2)
            x3 = grid_coord(node_3, 1)
            y3 = grid_coord(node_3, 2)
            x4 = grid_coord(node_4, 1)
            y4 = grid_coord(node_4, 2)

            call cal_axisym_quad_stiff_loc(x1, y1, x2, y2, x3, y3, x4, y4)

            do i = 1, 4
                do j = 1, 4
                    stiff_glo_mat(node_con_mat(elem, i), node_con_mat(elem, j)) = stiff_glo_mat(node_con_mat(elem, i), node_con_mat(elem, j)) + stiff_loc(i, j)
                end do
                ! TO BE DONE IN THE FUTURE
                ! heat_gen_glo_mat(node_con_mat(elem, i), 1) = heat_gen_glo_mat(node_con_mat(elem, i), 1) + heat_gen_loc(i, 1)
            end do
        end do
        stiff_glo_mat = 2 * pi * stiff_glo_mat
    end subroutine assemble_gsm
end module axisym_stiff_glo