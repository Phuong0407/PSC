module axisym_quad_elem
    implicit none
    private :: init_grad_mat, init_coord_mat, cal_DET22, cal_Mat22_INV, cal_sf_weight!, ! cal_sf_mat ! TO BE DONE IN THE FUTURE
    public :: cal_axisym_quad_stiff_loc, stiff_loc !, heat_gen_loc, heat_gen_intern, bound_grad ! TO BE DONE IN THE FUTURE

    ! real :: heat_gen_intern                                                   ! TO BE DONE IN THE FUTURE
    real, dimension(:, :), allocatable :: stiff_loc
    ! real, dimension(:, :), allocatable :: heat_gen_loc                        ! TO BE DONE IN THE FUTURE
    ! real, dimension(:, :), allocatable :: bound_grad # TO BE DONE IN THE FUTURE

contains
    function init_grad_mat(xi, eta) result(grad_mat)
        real, intent(in) :: xi, eta
        real, dimension(2, 4) :: grad_mat

        grad_mat(1, 1) = (eta - 1) / 4.0
        grad_mat(1, 2) = (1 - eta) / 4.0
        grad_mat(1, 3) = (1 + eta) / 4.0
        grad_mat(1, 4) = (-eta - 1) / 4.0

        grad_mat(2, 1) = (xi - 1) / 4.0
        grad_mat(2, 2) = (-xi - 1) / 4.0
        grad_mat(2, 3) = (1 + xi) / 4.0
        grad_mat(2, 4) = (1 - xi) / 4.0
    end function init_grad_mat

    function init_coord_mat(x1, y1, x2, y2, x3, y3, x4, y4) result(coord_mat)
        real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
        real, dimension(4, 2) :: coord_mat

        coord_mat(1, 1) = x1
        coord_mat(1, 2) = y1

        coord_mat(2, 1) = x2
        coord_mat(2, 2) = y2

        coord_mat(3, 1) = x3
        coord_mat(3, 2) = y3

        coord_mat(4, 1) = x4
        coord_mat(4, 2) = y4
    end function init_coord_mat

    function cal_sf_weight(y1, y2, y3, y4, xi, eta) result(SFW)
        real, intent(in) :: y1, y2, y3, y4, xi, eta
        real :: SFW

        SFW = 0.0
        SFW = SFW + (1 - xi) * (1 - eta) / 4.0 * y1
        SFW = SFW + (1 + xi) * (1 - eta) / 4.0 * y2
        SFW = SFW + (1 + xi) * (1 + eta) / 4.0 * y3
        SFW = SFW + (1 - xi) * (1 + eta) / 4.0 * y4
    end function cal_sf_weight

    ! TO BE DONE IN THE FUTURE
    ! function cal_sf_mat(xi, eta) result(SF_MAT)
    !     real, intent(in) :: xi, eta
    !     real, dimension(4, 1) :: SF_MAT

    !     SF_MAT(1, 1) = (1 - xi) * (1 - eta) / 4.0
    !     SF_MAT(2, 1) = (1 + xi) * (1 - eta) / 4.0
    !     SF_MAT(3, 1) = (1 + xi) * (1 + eta) / 4.0
    !     SF_MAT(4, 1) = (1 - xi) * (1 + eta) / 4.0
    ! end function cal_sf_mat

    function cal_DET22(Mat) result(DET22)
        real, dimension(2, 2), intent(in) :: Mat
        real :: DET22

        DET22 = Mat(1, 1) * Mat(2, 2) - Mat(1, 2) * Mat(2, 1)
    end function cal_DET22

    function cal_Mat22_INV(JMat) result(JMat_INV)
        real, dimension(2, 2), intent(in) :: JMat
        real, dimension(2, 2) :: JMat_INV
        real :: DET22, JMat11

        DET22 = cal_DET22(JMat)
        JMat11 = JMat(1, 1)
        JMat_INV(1, 1) = JMat(2, 2) / DET22
        JMat_INV(1, 2) = -JMat(1, 2) / DET22
        JMat_INV(2, 1) = -JMat(2, 1) / DET22
        JMat_INV(2, 2) = JMat11 / DET22
    end function cal_Mat22_INV

    subroutine cal_axisym_quad_stiff_loc(x1, y1, x2, y2, x3, y3, x4, y4)
        implicit none
        real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
    
        integer :: i
        real :: DET_JMat, SFW
        ! real, dimension(4, 1) :: SF_MAT                                       ! TO BE DONE IN THE FUTURE
        real, dimension(4, 2) :: coord_mat, gauss_points
        real, dimension(2, 4) :: grad_mat, B
        real, dimension(2, 2) :: JMat, INV_JMat
    
        gauss_points(1, 1) = -1.0/sqrt(3.0)
        gauss_points(1, 2) = -1.0/sqrt(3.0)

        gauss_points(2, 1) = -1.0/sqrt(3.0)
        gauss_points(2, 2) = 1.0/sqrt(3.0)

        gauss_points(3, 1) = 1.0/sqrt(3.0)
        gauss_points(3, 2) = -1.0/sqrt(3.0)

        gauss_points(4, 1) = 1.0/sqrt(3.0)
        gauss_points(4, 2) = 1.0/sqrt(3.0)
        
        coord_mat = init_coord_mat(x1, y1, x2, y2, x3, y3, x4, y4)

        allocate(stiff_loc(4, 4))
        ! allocate(heat_gen_loc(4, 1))                                          ! TO BE DONE IN THE FUTURE
        
        stiff_loc = 0.0
        ! heat_gen_loc = 0.0                                                    ! TO BE DONE IN THE FUTURE

        do i = 1, 4
            grad_mat = init_grad_mat(gauss_points(i, 1), gauss_points(i, 2))
            SFW = cal_sf_weight(y1, y2, y3, y4, gauss_points(i, 1), gauss_points(i, 2))
            JMat = matmul(grad_mat, coord_mat)
            DET_JMat = cal_DET22(JMat)
            INV_JMat = cal_Mat22_INV(JMat)
            B = matmul(INV_JMat, grad_mat)

            stiff_loc = stiff_loc + SFW * DET_JMat * matmul(transpose(B), B)

            ! SF_MAT = cal_sf_mat(gauss_points(i, 1), gauss_points(i, 2))       ! TO BE DONE IN THE FUTURE
            ! heat_gen_loc = heat_gen_loc + heat_gen_intern * SFW * SF_MAT      ! TO BE DONE IN THE FUTURE
        end do
    end subroutine cal_axisym_quad_stiff_loc
end module axisym_quad_elem