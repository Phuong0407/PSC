module rect_stiff_loc
    implicit none
    private
    public :: cal_rec_stiff_loc

contains
    function cal_rec_stiff_loc(x, y) result(K_loc_rect)
        implicit none
        real, intent(in) :: x, y
        real, dimension(4, 4) :: K_loc_rect

        K_loc_rect(1, 1) = (x / y + y / x) / 3.0
        K_loc_rect(1, 2) = (x / y - 2.0 * y / x) / 6.0
        K_loc_rect(1, 3) = -(x / y + y / x) / 6.0
        K_loc_rect(1, 4) = (y / x - 2.0 * x / y) / 6.0
        
        K_loc_rect(2, 1) = (x / y - 2.0 * y / x) / 6.0
        K_loc_rect(2, 2) = (x / y + y / x) / 3.0
        K_loc_rect(2, 3) = (y / x - 2.0 * x / y) / 6.0
        K_loc_rect(2, 4) = -(x / y + y / x) / 6.0
        
        K_loc_rect(3, 1) = -(x / y + y / x) / 6.0
        K_loc_rect(3, 2) = (y / x - 2.0 * x / y) / 6.0
        K_loc_rect(3, 3) = (x / y + y / x) / 3.0
        K_loc_rect(3, 4) = (x / y - 2.0 * y / x) / 6.0
        
        K_loc_rect(4, 1) = (y / x - 2.0 * x / y) / 6.0
        K_loc_rect(4, 2) = -(x / y + y / x) / 6.0
        K_loc_rect(4, 3) = (x / y - 2.0 * y / x) / 6.0
        K_loc_rect(4, 4) = (x / y + y / x) / 3.0
    end function cal_rec_stiff_loc
end module rect_stiff_loc