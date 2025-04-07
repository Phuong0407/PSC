module global_stiff
    use four_node_stiff_loc
    ! use tri_stiff_loc, TO BE DONE IN THE FUTURE
    implicit none
    private
    integer :: n_node, nv_node, nh_node, n_elem, nh_elem, nv_elem
    real, allocatable :: grid_coord(:,:)
    real, allocatable :: K_g(:,:)
    integer, allocatable :: node_con_mat(:,:)
    
    public :: init_gsm, assemble_gsm
    
contains
    subroutine init_gsm(n_node_in, nv_node_in, nh_node_in, n_elem_in, nh_elem_in, nv_elem_in, grid_coord_in, node_con_in)
        implicit none
        integer, intent(in) :: n_node_in, nh_node_in, nv_node_in
        integer, intent(in) :: n_elem_in, nh_elem_in, nv_elem_in
        real, intent(in) :: grid_coord_in(:,:)
        integer, intent(in) :: node_con_in(:,:)
        logical :: is_success

        allocate(grid_coord(n_node, 2))
        allocate(node_con_mat(n_node, 4))

        if (nh_node_in /= nh_elem_in + 1 .or. nv_node_in /= nv_elem_in + 1) then
            is_success = .false.
        else if (n_node_in /= nh_node_in * nv_node_in .or. n_elem_in /= nh_elem_in * nv_elem_in) then
            is_success = .false.
        else if (size(grid_coord_in, 1) /= n_node_in .or. size(node_con_in, 1) /= n_elem_in) then
            is_success = .false.
        else
            is_success = .true.
        endif

        n_node = n_node_in
        nv_node = nv_node_in
        nh_node = nh_node_in
        n_elem = n_elem_in
        nh_elem = nh_elem_in
        nv_elem = nv_elem_in

        allocate(grid_coord(n_node, 2))
        allocate(node_con_mat(n_elem, 4))
        allocate(K_g(n_node, n_node))
        
        grid_coord = grid_coord_in
        node_con_mat = node_con_in
        K_g = 0.0
    end subroutine init_gsm

    subroutine assemble_gsm()
        implicit none
        integer :: node_1, node_2, node_3, node_4
        integer :: i, j, elem
        real :: K_loc(4, 4)
        real :: x1, y1, x2, y2, x3, y3, x4, y4
        logical :: is_success

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

            K_loc = cal_four_node_stiff_loc(x1, y1, x2, y2, x3, y3, x4, y4)
            do i = 1, 4
                do j = 1, 4
                    K_g(node_con_mat(elem, i), node_con_mat(elem, j)) = K_g(node_con_mat(elem, i), node_con_mat(elem, j)) + K_loc(i, j)
                end do
            end do
        end do
    end subroutine assemble_gsm
end module global_stiff