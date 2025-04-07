!**************************************************************************************************
!**************************************************************************************************
! THIS IS THE PROGRAM TO SOLVE THE LAPLACE EQUATION FOR STREAMLINE FUNCTION
!**************************************************************************************************
!**************************************************************************************************
module streamline_solution
    implicit none
    private :: INIT_GRID, GEN_GRID_COORD, GEN_GRID_CONN, CALC_AXISYM_LSM_LOC, CALC_STIFF_MAT_GLO, CALC_DIRICHLET_BOUND, &
                MODIFY_STIFF_MAT_AND_FORCE_MAT, CALC_BELONGING_ELEM, IS_POINT_IN_QUAD
    real, parameter :: PI = 3.141592653589793
!**************************************************************************************************
!**************************************************************************************************
! INPUT DATA
    type :: grid_struct
        real :: H_1, H_2, H_3                               ! Horizontal dimension of each region
        real :: L_1, L_2                                    ! Vertical dimension of each region
        integer :: eH_1, eH_2, eH_3                         ! Horizontal NUMBER OF ELEMENT for each region
        integer :: indH_1, indH_2, indH_3                   ! Horizontal INDEX OF NODE for each region
        integer :: eV                                       ! Total vertical number of element
        integer :: eH                                       ! Total number of horizontal element
        integer :: nElem                                    ! Total number of element
        integer :: nV                                       ! Total number of vertical nodes
        integer :: nH                                       ! Total number of horizontal nodes
        integer :: nNode                                    ! Total number of nodes
        real, allocatable :: COORD(:,:)                     ! Coordinate matrices, NODE_IND[X, Y]
        integer, allocatable :: CONN(:,:)                   ! Connection Data for Assemble: ELEMENT[NODE1_IND, NODE_2_IND, NODE_3_IND, NODE_4_IND]
        integer, allocatable :: IND_DIRICHLET_BOUND(:)      ! Connection Data for Assemble: ELEMENT[NODE1_IND, NODE_2_IND, NODE_3_IND, NODE_4_IND]
        real :: dx_1, dx_2, dx_3                            ! Finite-Different horizontal dimension of each region
        real :: dy_1, dy_2                                  ! Finite-Different vertical dimension of each region
    contains
    end type grid_struct
    real :: U_0, U_1                                        ! Velocity at the Entrance and the Outlet of the Inlet
!**************************************************************************************************
!**************************************************************************************************
    public :: INIT_SOLUTION, INTERPOLATE_STREAMLINE_FUNC
    type(grid_struct) :: GRID
    real, allocatable :: STIFF_MAT_GLO(:,:)
    real, allocatable :: DIRICHLET_BOUND(:)
    real, allocatable :: FORCE_MAT_GLO(:)
    real, allocatable :: VALS(:)

contains
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! GRID DATA GENERATION
    subroutine INIT_GRID(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3)
        implicit none
        real, intent(in) :: iH_1, iH_2, iH_3, iL_1, iL_2
        integer, intent(in) :: iN, iM_1, iM_2, iM_3

        GRID%H_1 = iH_1
        GRID%H_2 = iH_2
        GRID%H_3 = iH_3
        GRID%L_1 = iL_1
        GRID%L_2 = iL_2
        GRID%eH_1 = iM_1
        GRID%eH_2 = iM_2
        GRID%eH_3 = iM_3

        GRID%eV = iN
        GRID%eH = GRID%eH_1 + GRID%eH_2 + GRID%eH_3

        GRID%dx_1 = GRID%H_1 / GRID%eH_1
        GRID%dx_2 = GRID%H_1 / GRID%eH_2
        GRID%dx_3 = GRID%H_1 / GRID%eH_3
        GRID%dy_1 = GRID%L_1 / GRID%eV
        GRID%dy_2 = GRID%L_2 / GRID%eV

        GRID%indH_1 = GRID%eH_1 + 1
        GRID%indH_2 = GRID%eH_1 + GRID%eH_2
        GRID%indH_3 = GRID%eH_1 + GRID%eH_2 + GRID%eH_3 + 1

        GRID%nV = GRID%eV + 1
        GRID%nH = GRID%eH + 1
        GRID%nElem = GRID%eH * GRID%eV
        GRID%nNode = GRID%nH * GRID%nV
    end subroutine INIT_GRID

    subroutine GEN_GRID_COORD()
        implicit none
        integer :: i, j, ind_node
        real :: x, y, a, b

        allocate(GRID%COORD(GRID%nNode, 2))
        do i = 1, GRID%nV
            do j = 1, GRID%indH_1
                ind_node = (i - 1) * GRID%nH + j
                x = GRID%dx_1 * (j - 1)
                y = GRID%dy_1 * (i - 1)
                GRID%COORD(ind_node, 1) = x
                GRID%COORD(ind_node, 2) = y
            end do

            a = ((GRID%dy_2 - GRID%dy_1) * (i - 1) + GRID%L_1 - GRID%L_2)/GRID%H_2
            b = GRID%dy_1 * (i - 1) - a * GRID%H_1
            do j = GRID%indH_1 + 1, GRID%indH_2
                ind_node = (i - 1) * GRID%nH + j
                x = GRID%H_1 + GRID%dx_2 * (j - GRID%indH_1)
                y = a * x + b
                GRID%COORD(ind_node, 1) = x
                GRID%COORD(ind_node, 2) = y
            end do

            a = GRID%H_1 + GRID%H_2
            b = GRID%L_1 - GRID%L_2
            do j = GRID%indH_2 + 1, GRID%indH_3
                ind_node = (i - 1) * GRID%nH + j
                x = a + GRID%dx_3 * (j - 1 - GRID%indH_2)
                y = b + GRID%dy_2 * (i - 1)
                GRID%COORD(ind_node, 1) = x
                GRID%COORD(ind_node, 2) = y
            end do
        end do
    end subroutine GEN_GRID_COORD

    subroutine GEN_GRID_CONN()
        implicit none
        integer :: i, j, ind_elem, IND_FIRST_NODE, IND_SECOND_NODE, IND_THIRD_NODE, IND_FOURTH_NODE

        allocate(GRID%CONN(GRID%nElem, 4))
        do i = 1, GRID%eV
            do j = 1, GRID%eH
                ind_elem = (i - 1) * GRID%eH + j
                IND_FIRST_NODE = GRID%nH * i + j - GRID%nH
                IND_SECOND_NODE = IND_FIRST_NODE + 1
                IND_THIRD_NODE = IND_SECOND_NODE + GRID%nH
                IND_FOURTH_NODE = IND_FIRST_NODE + GRID%nH
                GRID%CONN(ind_elem, 1) = IND_FIRST_NODE
                GRID%CONN(ind_elem, 2) = IND_SECOND_NODE
                GRID%CONN(ind_elem, 3) = IND_THIRD_NODE
                GRID%CONN(ind_elem, 4) = IND_FOURTH_NODE
            end do
        end do
    end subroutine GEN_GRID_CONN
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! VELOCITY INITIALIZATION
    subroutine INIT_VELOCITY(U)
        implicit none
        real, intent(in) :: U
        U_0 = U
        U_1 = U_0 * (GRID%L_1**2.0)/(GRID%L_1**2.0 - (GRID%L_1 - GRID%L_2)**2.0)
    end subroutine  INIT_VELOCITY

!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! LOCAL STIFFNESS MATRIX CALCULATION
    function CALC_AXISYM_LSM_LOC(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4) result(STIFF_MAT_LOC)
        implicit none
        integer :: i, j
        real, intent(in) :: x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4
        real :: JACOBIAN_DET, SFW
        real, dimension(4, 2) :: ELEM_COORD_MAT
        real, dimension(2, 4) :: GRAD_MAT, SF_DEV
        real, dimension(2, 2) :: JACOBIAN_MAT, JACOBIAN_INV
        real, dimension(4, 4) :: STIFF_MAT_LOC
        real, dimension(2) :: XI, ETA
        
        XI(1) = -1.0/sqrt(3.0)
        XI(2) = 1.0/sqrt(3.0)
        ETA(1) = -1.0/sqrt(3.0)
        ETA(2) = 1.0/sqrt(3.0)
        
        ELEM_COORD_MAT(1, 1) = x_1
        ELEM_COORD_MAT(1, 2) = y_1
        ELEM_COORD_MAT(2, 1) = x_2
        ELEM_COORD_MAT(2, 2) = y_2
        ELEM_COORD_MAT(3, 1) = x_3
        ELEM_COORD_MAT(3, 2) = y_3
        ELEM_COORD_MAT(4, 1) = x_4
        ELEM_COORD_MAT(4, 2) = y_4

        STIFF_MAT_LOC = 0.0

        do i = 1, 2
            do j = 1, 2
                GRAD_MAT(1, 1) = (ETA(j) - 1) / 4.0
                GRAD_MAT(1, 2) = (1 - ETA(j)) / 4.0
                GRAD_MAT(1, 3) = (1 + ETA(j)) / 4.0
                GRAD_MAT(1, 4) = (-ETA(j) - 1) / 4.0
                GRAD_MAT(2, 1) = (XI(i) - 1) / 4.0
                GRAD_MAT(2, 2) = (-XI(i) - 1) / 4.0
                GRAD_MAT(2, 3) = (1 + XI(i)) / 4.0
                GRAD_MAT(2, 4) = (1 - XI(i)) / 4.0

                SFW = - (1 - XI(i)) * GRAD_MAT(1, 1) * y_1 &
                        + (1 + XI(i)) * GRAD_MAT(1, 2) * y_2 &
                        + (1 + XI(i)) * GRAD_MAT(1, 3) * y_3 &
                        - (1 - XI(i)) * GRAD_MAT(1, 4) * y_4

                JACOBIAN_MAT = matmul(GRAD_MAT, ELEM_COORD_MAT)
                JACOBIAN_DET = JACOBIAN_MAT(1, 1) * JACOBIAN_MAT(2, 2) &
                                - JACOBIAN_MAT(1, 2) * JACOBIAN_MAT(2, 1)
                
                JACOBIAN_INV(1, 1) = JACOBIAN_MAT(2, 2) / JACOBIAN_DET
                JACOBIAN_INV(1, 2) = - JACOBIAN_MAT(1, 2) / JACOBIAN_DET
                JACOBIAN_INV(2, 1) = - JACOBIAN_MAT(2, 1) / JACOBIAN_DET
                JACOBIAN_INV(2, 2) = JACOBIAN_MAT(1, 1) / JACOBIAN_DET

                SF_DEV = matmul(JACOBIAN_INV, GRAD_MAT)

                STIFF_MAT_LOC = STIFF_MAT_LOC + 2 * PI * SFW * JACOBIAN_DET * matmul(transpose(SF_DEV), SF_DEV)
            end do
        end do
    end function CALC_AXISYM_LSM_LOC
!**************************************************************************************************
!**************************************************************************************************
! GLOBAL SYSTEM OF EQUATIONS
    subroutine CALC_STIFF_MAT_GLO()
        implicit none
        integer :: i, j, LOC_I, LOC_J, GLO_I, GLO_J
        real, allocatable :: STIFF_MAT_LOC(:,:), FORCE_MAT_LOC(:)
        integer :: ind_elem, ind_node_1, ind_node_2, ind_node_3, ind_node_4
        real :: x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4

        allocate(STIFF_MAT_LOC(4,4))
        allocate(FORCE_MAT_LOC(4))
        allocate(STIFF_MAT_GLO(GRID%nNode, GRID%nNode))

        STIFF_MAT_GLO = 0.0
        STIFF_MAT_LOC = 0.0
        FORCE_MAT_LOC = 0.0

        do i = 1, GRID%eV
            do j = 1, GRID%eH
                ind_elem = (i - 1) * GRID%eH + j
                ind_node_1 = GRID%CONN(ind_elem, 1)
                ind_node_2 = GRID%CONN(ind_elem, 2)
                ind_node_3 = GRID%CONN(ind_elem, 3)
                ind_node_4 = GRID%CONN(ind_elem, 4)

                x_1 = GRID%COORD(ind_node_1, 1)
                y_1 = GRID%COORD(ind_node_1, 2)
                x_2 = GRID%COORD(ind_node_2, 1)
                y_2 = GRID%COORD(ind_node_2, 2)
                x_3 = GRID%COORD(ind_node_3, 1)
                y_3 = GRID%COORD(ind_node_3, 2)
                x_4 = GRID%COORD(ind_node_4, 1)
                y_4 = GRID%COORD(ind_node_4, 2)
                
                STIFF_MAT_LOC = CALC_AXISYM_LSM_LOC(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4)

                do LOC_I = 1, 4
                    GLO_I = GRID%CONN(ind_elem, LOC_I)
                    do LOC_J = 1, 4
                        GLO_J = GRID%CONN(ind_elem, LOC_J)
                        STIFF_MAT_GLO(GLO_I, GLO_J) = STIFF_MAT_GLO(GLO_I, GLO_J) + STIFF_MAT_LOC(LOC_I, LOC_J)
                    end do
                end do
            end do
        end do
    end subroutine CALC_STIFF_MAT_GLO
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! CALCULATE DIRICHLET BOUNDARY CONDITION
    subroutine CALC_DIRICHLET_BOUND()
        implicit none
        integer :: i, k, IND_LOWER_BOUND, IND_UPPER_BOUND, IND_LEFT_BOUND, IND_RIGHT_BOUND
        real :: PSI_UPPER, PSI_LEFT, PSI_RIGHT, R_1, R_2, R_3

        allocate(DIRICHLET_BOUND(GRID%nNode))
        allocate(GRID%IND_DIRICHLET_BOUND(2 * GRID%nH + 2 * GRID%nV - 4))

        DIRICHLET_BOUND = 0.0
        PSI_UPPER = 0.5 * U_0 * (GRID%L_1)**2.0
        k = 1
        do i = 1, GRID%nH
            IND_LOWER_BOUND = i
            GRID%IND_DIRICHLET_BOUND(k) = IND_LOWER_BOUND
            k = k + 1
        end do

        R_3 = GRID%L_1 - GRID%L_2
        do i = 2, GRID%eV
            IND_LEFT_BOUND = (i - 1) * GRID%nH + 1
            R_1 = GRID%dy_1 * (i - 1)
            PSI_LEFT = 0.5 * U_0 * R_1**2.0
            DIRICHLET_BOUND(IND_LEFT_BOUND) = PSI_LEFT
            GRID%IND_DIRICHLET_BOUND(k) = IND_LEFT_BOUND
            k = k + 1
            
            IND_RIGHT_BOUND = i * GRID%nH
            R_2 = GRID%dy_2 * (i - 1) + R_3
            PSI_RIGHT = 0.5 * U_1 * (R_2**2.0 - R_3**2.0)
            DIRICHLET_BOUND(IND_RIGHT_BOUND) = PSI_RIGHT
            GRID%IND_DIRICHLET_BOUND(k) = IND_RIGHT_BOUND
            k = k + 1
        end do

        do i = 1, GRID%nH
            IND_UPPER_BOUND = (GRID%nV - 1) * GRID%nH + i
            DIRICHLET_BOUND(IND_UPPER_BOUND) = PSI_UPPER
            GRID%IND_DIRICHLET_BOUND(k) = IND_UPPER_BOUND
            k = k + 1
        end do
        if (k /= size(GRID%IND_DIRICHLET_BOUND) + 1) then
            print *, "FAIL TO INITIALIZE IND_DIRICHLET_BOUND"
            stop
        end if
    end subroutine CALC_DIRICHLET_BOUND
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! MODIFY GLOBAL STIFFNESS MATRIX AND FORCE MATRIX
    subroutine MODIFY_STIFF_MAT_AND_FORCE_MAT()
        implicit none
        integer :: IND, i, j

        allocate(FORCE_MAT_GLO(GRID%nNode))
        FORCE_MAT_GLO = 0.0

        do i = 1, size(GRID%IND_DIRICHLET_BOUND)
            IND = GRID%IND_DIRICHLET_BOUND(i) ! INDEX OF COLUMN TO BE REMOVED
            do j = 1, GRID%nNode
                FORCE_MAT_GLO(j) = FORCE_MAT_GLO(j) - DIRICHLET_BOUND(IND) * STIFF_MAT_GLO(j, IND)
            end do
            STIFF_MAT_GLO(:, IND) = 0.0
            STIFF_MAT_GLO(IND, :) = 0.0
            STIFF_MAT_GLO(IND, IND) = 1.0
            FORCE_MAT_GLO(IND) = DIRICHLET_BOUND(IND)
        end do
    end subroutine MODIFY_STIFF_MAT_AND_FORCE_MAT
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
    subroutine INIT_SOLUTION(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3, iU_0)
        implicit none
        real, intent(in) :: iH_1, iH_2, iH_3, iL_1, iL_2, iU_0
        integer, intent(in) :: iN, iM_1, iM_2, iM_3
        call INIT_GRID(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3)
        call GEN_GRID_COORD()
        call GEN_GRID_CONN()
        call INIT_VELOCITY(iU_0)
        call CALC_STIFF_MAT_GLO()
        call CALC_DIRICHLET_BOUND()
        call MODIFY_STIFF_MAT_AND_FORCE_MAT()
    end subroutine INIT_SOLUTION
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! STREAMLINE FUNCTION INTERPOLATION
    function CALC_BELONGING_ELEM(x, y) result(IND_ELEM)
        implicit none
        real, intent(in) :: x, y
        integer :: ind_elem
        integer :: ind_node_1, ind_node_2, ind_node_3, ind_node_4
        integer :: i, j
        logical :: is_inside
        real :: x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4
    
        ind_elem = -1
    
        do i = 1, GRID%eV
            do j = 1, GRID%nH
                ind_elem = (i - 1) * GRID%nH + j
    
                ind_node_1 = GRID%CONN(ind_elem, 1)
                ind_node_2 = GRID%CONN(ind_elem, 2)
                ind_node_3 = GRID%CONN(ind_elem, 3)
                ind_node_4 = GRID%CONN(ind_elem, 4)
    
                x_1 = GRID%COORD(ind_node_1, 1)
                y_1 = GRID%COORD(ind_node_1, 2)
                x_2 = GRID%COORD(ind_node_2, 1)
                y_2 = GRID%COORD(ind_node_2, 2)
                x_3 = GRID%COORD(ind_node_3, 1)
                y_3 = GRID%COORD(ind_node_3, 2)
                x_4 = GRID%COORD(ind_node_4, 1)
                y_4 = GRID%COORD(ind_node_4, 2)
    
                if (x_1 <= x .and. x <= x_2) then
                    if (y_1 == y_2 .and. y_3 == y_4 .and. y_1 <= y .and. y <= y_4) then
                        return
                    else
                        is_inside = IS_POINT_IN_QUAD(x, y, x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4)
                        if (is_inside) then
                            return
                        end if
                    end if
                end if
            end do
        end do
    end function CALC_BELONGING_ELEM
    
    logical function IS_POINT_IN_QUAD(x, y, x1, y1, x2, y2, x3, y3, x4, y4)
        implicit none
        real, intent(in) :: x, y
        real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
        real :: cross(4)

        cross(1) = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
        cross(2) = (x3 - x2) * (y - y2) - (y3 - y2) * (x - x2)
        cross(3) = (x4 - x3) * (y - y3) - (y4 - y3) * (x - x3)
        cross(4) = (x1 - x4) * (y - y4) - (y1 - y4) * (x - x4)

        IS_POINT_IN_QUAD = all(cross >= 0.0) .or. all(cross <= 0.0)
    end function IS_POINT_IN_QUAD

    function INTERPOLATE_STREAMLINE_FUNC(x, y) result(STREAMLINE_FUNC)
        real, intent(in) :: x, y
        real :: STREAMLINE_FUNC
        real :: x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4
        integer :: ind_elem
        integer :: ind_node_1, ind_node_2, ind_node_3, ind_node_4
        real :: PSI_1, PSI_2, PSI_3, PSI_4
        real :: xi, eta

        ind_elem = CALC_BELONGING_ELEM(x, y)
        ind_node_1 = GRID%CONN(ind_elem, 1)
        ind_node_2 = GRID%CONN(ind_elem, 2)
        ind_node_3 = GRID%CONN(ind_elem, 3)
        ind_node_4 = GRID%CONN(ind_elem, 4)

        x_1 = GRID%COORD(ind_node_1, 1)
        y_1 = GRID%COORD(ind_node_1, 2)
        x_2 = GRID%COORD(ind_node_2, 1)
        y_2 = GRID%COORD(ind_node_2, 2)
        x_3 = GRID%COORD(ind_node_3, 1)
        y_3 = GRID%COORD(ind_node_3, 2)
        x_4 = GRID%COORD(ind_node_4, 1)
        y_4 = GRID%COORD(ind_node_4, 2)

        PSI_1 = VALS(ind_node_1)
        PSI_2 = VALS(ind_node_2)
        PSI_3 = VALS(ind_node_3)
        PSI_4 = VALS(ind_node_4)

        if (.not. (x == x_1 .and. y == y_1) .and. & 
            .not. (x == x_2 .and. y == y_2) .and. &
            .not. (x == x_3 .and. y == y_3) .and. &
            .not. (x == x_4 .and. y == y_4)) then

            xi = ((x - x_1) * (y_4 - y_1) - (y - y_1) * (x_4 - x_1)) / &
                ((x_2 - x_1) * (y_4 - y_1) - (y_2 - y_1) * (x_4 - x_1))
            eta = ((x - x_1) * (y_2 - y_1) - (y - y_1) * (x_2 - x_1)) / &
                    ((x_4 - x_1) * (y_2 - y_1) - (y_4 - y_1) * (x_2 - x_1))
            STREAMLINE_FUNC = PSI_1 * (1 - xi) * (1 - eta) + &
                                PSI_2 * xi * (1 - eta) + &
                                PSI_3 * xi * eta + &
                                PSI_4 * (1 - xi) * eta
        else
            if (x == x_1 .and. y == y_1) then
                STREAMLINE_FUNC = PSI_1
            else if (x == x_2 .and. y == y_2) then
                STREAMLINE_FUNC = PSI_2
            else if (x == x_3 .and. y == y_3) then
                STREAMLINE_FUNC = PSI_3
            else
                STREAMLINE_FUNC = PSI_4
            end if
        end if
    end function INTERPOLATE_STREAMLINE_FUNC
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
end module streamline_solution