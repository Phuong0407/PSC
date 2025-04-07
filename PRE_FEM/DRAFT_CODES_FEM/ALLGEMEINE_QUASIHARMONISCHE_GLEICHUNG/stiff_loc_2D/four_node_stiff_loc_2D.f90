module four_node_stiff_loc
    use rect_stiff_loc_2D
    use quad_stiff_loc_2D
    implicit none
    private :: is_rect_grid
    public :: cal_four_node_stiff_loc

contains
    function is_rect_grid(x1, y1, x2, y2, x3, y3, x4, y4) result(is_rect)
        real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
        logical :: is_rect

        if (x1 == x4 .and. x2 == x3 .and. y1 == y2 .and. y3 == y4 ) then
            is_rect = .true.
        else
            is_rect = .false.
        endif
    end function is_rect_grid

    function cal_four_node_stiff_loc(x1, y1, x2, y2, x3, y3, x4, y4) result(K_loc)
        real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
        real :: K_loc(4, 4)

        if (is_rect_grid(x1, y1, x2, y2, x3, y3, x4, y4)) then
            K_loc = cal_rec_stiff_loc(x2 - x1, y3 - y2)
        else
            K_loc = cal_quad_stiff_loc(x1, y1, x2, y2, x3, y3, x4, y4)
        endif
    end function cal_four_node_stiff_loc
end module four_node_stiff_loc