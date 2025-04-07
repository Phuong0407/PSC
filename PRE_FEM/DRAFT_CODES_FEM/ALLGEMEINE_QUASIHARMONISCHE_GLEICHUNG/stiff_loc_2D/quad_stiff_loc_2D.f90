module quad_stiff_loc
    implicit none
    private :: init_grad_mat, init_coord_mat, cal_DET22, cal_Mat22_INV
    public :: cal_quad_stiff_loc

contains

    function init_grad_mat(xi, eta) result(result)
        real, intent(in) :: xi, eta
        real, dimension(2, 4) :: result

        result(1, 1) = (eta - 1) / 4.0
        result(1, 2) = (1 - eta) / 4.0
        result(1, 3) = (1 + eta) / 4.0
        result(1, 4) = (-eta - 1) / 4.0

        result(2, 1) = (xi - 1) / 4.0
        result(2, 2) = (-xi - 1) / 4.0
        result(2, 3) = (1 + xi) / 4.0
        result(2, 4) = (1 - xi) / 4.0
    end function init_grad_mat

    function init_coord_mat(x1, y1, x2, y2, x3, y3, x4, y4) result(result)
        real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
        real, dimension(4, 2) :: result

        result(1, 1) = x1
        result(1, 2) = y1

        result(2, 1) = x2
        result(2, 2) = y2

        result(3, 1) = x3
        result(3, 2) = y3

        result(4, 1) = x4
        result(4, 2) = y4
    end function init_coord_mat

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

    function cal_quad_stiff_loc(x1, y1, x2, y2, x3, y3, x4, y4) result(K_loc_quad)
        implicit none
        real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
        real, dimension(4, 4) :: K_loc_quad
    
        real, dimension(4, 2) :: coord_mat
        real, dimension(4, 2) :: gauss_points
        real, dimension(2, 4) :: grad_mat, B
        real, dimension(2, 2) :: JMat, INV_JMat
        real :: DET_JMat
        integer :: i, j
    
        gauss_points(1, 1) = -1.0/sqrt(3.0)
        gauss_points(1, 2) = -1.0/sqrt(3.0)

        gauss_points(2, 1) = -1.0/sqrt(3.0)
        gauss_points(2, 2) = 1.0/sqrt(3.0)

        gauss_points(3, 1) = 1.0/sqrt(3.0)
        gauss_points(3, 2) = -1.0/sqrt(3.0)

        gauss_points(4, 1) = 1.0/sqrt(3.0)
        gauss_points(4, 2) = 1.0/sqrt(3.0)

    
        coord_mat = init_coord_mat(x1, y1, x2, y2, x3, y3, x4, y4)
        
        print *, "Coordinate Matrix Matrix:"
        do j = 1, 4
            write(*, '(4F10.4)') coord_mat(j, :)
        end do

        K_loc_quad = 0.0
    
        do i = 1, 4
            print *, "Gaussian Point: ", gauss_points(i, :)

            grad_mat = init_grad_mat(gauss_points(i, 1), gauss_points(i, 2))

            print *, "Gradient Matrix:"
            do j = 1, 2
                write(*, '(4F10.4)') grad_mat(j, :)
            end do

            JMat = matmul(grad_mat, coord_mat)
            DET_JMat = cal_DET22(JMat)
            INV_JMat = cal_Mat22_INV(JMat)
            B = matmul(INV_JMat, grad_mat)
            
            K_loc_quad = K_loc_quad + DET_JMat * matmul(transpose(B), B)
        end do
    end function cal_quad_stiff_loc

end module quad_stiff_loc