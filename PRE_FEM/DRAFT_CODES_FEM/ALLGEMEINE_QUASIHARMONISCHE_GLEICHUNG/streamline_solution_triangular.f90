!**************************************************************************************************
!**************************************************************************************************
! THIS IS THE PROGRAM TO SOLVE THE LAPLACE EQUATION FOR STREAMLINE FUNCTION
!**************************************************************************************************
!**************************************************************************************************
module streamline_solution_triangular
    implicit none
    private :: INIT_GRID, GEN_GRID_COORD, GEN_GRID_CONN, &
                CALC_AXISYM_LSM_LOC, CALC_LSM_LOC, CALC_STIFF_MAT_GLO, &
                CALC_DIRICHLET_BOUND, CALC_DIRICHLET_BOUND_AXISYM, &
                MODIFY_STIFF_MAT_AND_FORCE_MAT, &
                CALC_BELONGING_ELEM, IS_POINT_IN_TRIANGLE, CALC_TRIANGLE_AREA, &
                INIT_VELOCITY, INIT_VELOCITY_AXISYM
    real, parameter :: PI = 3.141592653589793
!**************************************************************************************************
!**************************************************************************************************
! INPUT DATA
    type :: grid_struct
        real :: H_1, H_2, H_3                      ! Horizontal dimension of each region
        real :: L_1, L_2                           ! Vertical dimension of each region
        integer :: eH_1, eH_2, eH_3                         ! Horizontal NUMBER OF ELEMENT for each region
        integer :: indH_1, indH_2, indH_3                   ! Horizontal INDEX OF NODE for each region
        integer :: eV                                       ! Total vertical number of element
        integer :: eH                                       ! Total number of horizontal element
        integer :: nElem                                    ! Total number of element
        integer :: n_Tri_Elem                               ! Total number of TRIANGLULAR element
        integer :: nV                                       ! Total number of vertical nodes
        integer :: nH                                       ! Total number of horizontal nodes
        integer :: nNode                                    ! Total number of nodes
        real, allocatable :: COORD(:,:)            ! Coordinate matrices, NODE_IND[X, Y]
        integer, allocatable :: CONN(:,:)                   ! Connection Data for Assemble: ELEMENT[NODE1_IND, NODE_2_IND, NODE_3_IND, NODE_4_IND]
        integer, allocatable :: IND_DIRICHLET_BOUND(:)      ! Connection Data for Assemble: ELEMENT[NODE1_IND, NODE_2_IND, NODE_3_IND, NODE_4_IND]
        real :: dx_1, dx_2, dx_3                   ! Finite-Different horizontal dimension of each region
        real :: dy_1, dy_2                         ! Finite-Different vertical dimension of each region
    contains
    end type grid_struct
    real :: U_0, U_1                               ! Velocity at the Entrance and the Outlet of the Inlet
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
        GRID%nNode = GRID%nH * GRID%nV

        GRID%nElem = GRID%eH * GRID%eV
        ! UPDATE HERE AS WE ADAPT THE QUADRILATERAL MESH TO THE TRIANGULAR MESH,
        ! THEN THE NUMBER OF NODE IS UNCHANGED, BUT THE NUMBER OF ELEMENT INCREASE DOUBLY
        GRID%n_Tri_Elem = 2 * GRID%nElem
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
        integer :: i, j, ind_elem
        integer :: IND_FIRST_NODE, IND_SECOND_NODE, IND_THIRD_NODE, IND_FOURTH_NODE
        integer :: TRI_ELEM_IND

        allocate(GRID%CONN(GRID%n_Tri_Elem, 3))

        TRI_ELEM_IND = 1

        do i = 1, GRID%eV
            do j = 1, GRID%eH
                ind_elem = (i - 1) * GRID%eH + j

                IND_FIRST_NODE = GRID%nH * i + j - GRID%nH
                IND_SECOND_NODE = IND_FIRST_NODE + 1
                IND_THIRD_NODE = IND_SECOND_NODE + GRID%nH
                IND_FOURTH_NODE = IND_FIRST_NODE + GRID%nH

                GRID%CONN(TRI_ELEM_IND, 1) = IND_FIRST_NODE
                GRID%CONN(TRI_ELEM_IND, 2) = IND_SECOND_NODE
                GRID%CONN(TRI_ELEM_IND, 3) = IND_THIRD_NODE
                TRI_ELEM_IND = TRI_ELEM_IND + 1

                GRID%CONN(TRI_ELEM_IND, 1) = IND_FIRST_NODE
                GRID%CONN(TRI_ELEM_IND, 2) = IND_THIRD_NODE
                GRID%CONN(TRI_ELEM_IND, 3) = IND_FOURTH_NODE
                TRI_ELEM_IND = TRI_ELEM_IND + 1
            end do
        end do
        if (TRI_ELEM_IND /= GRID%n_Tri_Elem + 1) then
            print *, "FAIL TO INITIALIZE CONNECTION MATRIX"
            stop
        end if
    end subroutine GEN_GRID_CONN
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! VELOCITY INITIALIZATION
    subroutine INIT_VELOCITY_AXISYM(U)
        implicit none
        real, intent(in) :: U
        U_0 = U
        U_1 = U_0 * (GRID%L_1**2.0)/(GRID%L_1**2.0 - (GRID%L_1 - GRID%L_2)**2.0)
    end subroutine INIT_VELOCITY_AXISYM

    subroutine INIT_VELOCITY(U)
        implicit none
        real, intent(in) :: U
        U_0 = U
        U_1 = U_0 * GRID%L_1/GRID%L_2
    end subroutine INIT_VELOCITY

!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! LOCAL STIFFNESS MATRIX CALCULATION
    function CALC_AXISYM_LSM_LOC(z_1, r_1, z_2, r_2, z_3, r_3) result(STIFF_MAT_LOC)
        implicit none
        real, intent(in) :: z_1, r_1, z_2, r_2, z_3, r_3
        real, allocatable :: STIFF_MAT_LOC(:,:)
        real :: c(3), d(3)
        real :: A_elem, r_centroid
        integer :: i, j

        ! b(1) = r_2 * z_3 - r_3 * z_2
        ! b(2) = r_3 * z_1 - r_1 * z_3
        ! b(3) = r_1 * z_2 - r_2 * z_1
        c(1) = z_2 - z_3
        c(2) = z_3 - z_1
        c(3) = z_1 - z_2
        d(1) = r_3 - r_1
        d(2) = r_1 - r_3
        d(3) = r_2 - r_1

        A_elem = CALC_TRIANGLE_AREA(z_1, r_1, z_2, r_2, z_3, r_3)
        r_centroid = (r_1 + r_2 + r_3) / 3.0
        allocate(STIFF_MAT_LOC(3, 3))

        do i = 1, 3
            do j = 1, 3
                STIFF_MAT_LOC(i, j) = PI * r_centroid / (2.0 * A_elem) * (c(i) * c(j) + d(i) * d(j))
            end do
        end do
    end function CALC_AXISYM_LSM_LOC

    function CALC_LSM_LOC(x_1, y_1, x_2, y_2, x_3, y_3) result(STIFF_MAT_LOC)
        implicit none
        real, intent(in) :: x_1, y_1, x_2, y_2, x_3, y_3
        real, allocatable :: STIFF_MAT_LOC(:,:)
        real, dimension(3) :: b, c
        real :: AREA
        integer :: i, j

        AREA = CALC_TRIANGLE_AREA(x_1, y_1, x_2, y_2, x_3, y_3)
        b(1) = y_2 - y_3
        b(2) = y_3 - y_1
        b(3) = y_1 - y_2
        c(1) = x_3 - x_2
        c(2) = x_1 - x_3
        c(3) = x_2 - x_1
        
        allocate(STIFF_MAT_LOC(3,3))
        do i = 1, 3
            do j = 1, 3
                STIFF_MAT_LOC(i,j) = (b(i) * b(j) + c(i) * c(j))/(4.0 * AREA)
            end do
        end do
    end function CALC_LSM_LOC
!**************************************************************************************************
!**************************************************************************************************
! GLOBAL SYSTEM OF EQUATIONS
    subroutine CALC_STIFF_MAT_GLO()
        implicit none
        integer :: i, LOC_I, LOC_J, GLO_I, GLO_J
        ! integer :: j, k
        real, allocatable :: STIFF_MAT_LOC(:,:)
        integer :: ind_node_1, ind_node_2, ind_node_3, ind_node_4
        real :: x_1, y_1, x_2, y_2, x_3, y_3

        allocate(STIFF_MAT_LOC(3,3))
        allocate(STIFF_MAT_GLO(GRID%nNode, GRID%nNode))

        STIFF_MAT_GLO = 0.0
        STIFF_MAT_LOC = 0.0

        do i = 1, GRID%n_Tri_Elem
            ind_node_1 = GRID%CONN(i, 1)
            ind_node_2 = GRID%CONN(i, 2)
            ind_node_3 = GRID%CONN(i, 3)

            x_1 = GRID%COORD(ind_node_1, 1)
            y_1 = GRID%COORD(ind_node_1, 2)
            x_2 = GRID%COORD(ind_node_2, 1)
            y_2 = GRID%COORD(ind_node_2, 2)
            x_3 = GRID%COORD(ind_node_3, 1)
            y_3 = GRID%COORD(ind_node_3, 2)
            
            ! STIFF_MAT_LOC = CALC_AXISYM_LSM_LOC(z_1, r_1, z_2, r_2, z_3, r_3)
            STIFF_MAT_LOC = CALC_LSM_LOC(x_1, y_1, x_2, y_2, x_3, y_3)
            
            ! print *, "Local Matrix Element", i
            ! do k = 1, size(STIFF_MAT_LOC, 1)  ! Loop over rows
            !     do j = 1, size(STIFF_MAT_LOC, 2)  ! Loop over columns
            !         write(*, '(F10.4)', advance='no') STIFF_MAT_LOC(k, j)  ! Print each element
            !     end do
            !     write(*, *)  ! Move to the next line after printing a row
            ! end do
            ! print *, "End Local Matrix Element", i

            do LOC_I = 1, 3
                GLO_I = GRID%CONN(i, LOC_I)
                do LOC_J = 1, 3
                    GLO_J = GRID%CONN(i, LOC_J)
                    STIFF_MAT_GLO(GLO_I, GLO_J) = STIFF_MAT_GLO(GLO_I, GLO_J) + STIFF_MAT_LOC(LOC_I, LOC_J)
                end do
            end do
        end do

        ! do i = 1, size(STIFF_MAT_GLO, 1)  ! Loop over rows
        !     do j = 1, size(STIFF_MAT_GLO, 2)  ! Loop over columns
        !         write(*, '(F10.4)', advance='no') STIFF_MAT_GLO(i, j)  ! Print each element
        !     end do
        !     write(*, *)  ! Move to the next line after printing a row
        ! end do
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
            PSI_LEFT = U_0 * R_1
            DIRICHLET_BOUND(IND_LEFT_BOUND) = PSI_LEFT
            GRID%IND_DIRICHLET_BOUND(k) = IND_LEFT_BOUND
            k = k + 1
            
            IND_RIGHT_BOUND = i * GRID%nH
            R_2 = GRID%dy_2 * (i - 1) + R_3
            PSI_RIGHT = U_1 * (R_2 - R_3)
            DIRICHLET_BOUND(IND_RIGHT_BOUND) = PSI_RIGHT
            GRID%IND_DIRICHLET_BOUND(k) = IND_RIGHT_BOUND
            k = k + 1
        end do

        PSI_UPPER = U_0 * GRID%L_1
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

    subroutine CALC_DIRICHLET_BOUND_AXISYM()
        implicit none
        integer :: i, k, IND_LOWER_BOUND, IND_UPPER_BOUND, IND_LEFT_BOUND, IND_RIGHT_BOUND
        real :: PSI_UPPER, PSI_LEFT, PSI_RIGHT, R_1, R_2, R_3

        allocate(DIRICHLET_BOUND(GRID%nNode))
        allocate(GRID%IND_DIRICHLET_BOUND(2 * GRID%nH + 2 * GRID%nV - 4))

        DIRICHLET_BOUND = 0.0
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

        PSI_UPPER = 0.5 * U_0 * (GRID%L_1)**2.0
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
    end subroutine CALC_DIRICHLET_BOUND_AXISYM
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
        ! call INIT_VELOCITY_AXISYM(iU_0)
        call CALC_STIFF_MAT_GLO()
        call CALC_DIRICHLET_BOUND()
        ! call CALC_DIRICHLET_BOUND_AXISYM()
        call MODIFY_STIFF_MAT_AND_FORCE_MAT()
    end subroutine INIT_SOLUTION
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
! STREAMLINE FUNCTION INTERPOLATION
    function CALC_TRIANGLE_AREA(x1, y1, x2, y2, x3, y3) result(triangle_area)
        implicit none
        real, intent(in) :: x1, y1, x2, y2, x3, y3
        real :: triangle_area

        triangle_area = 0.5 * ABS((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)))
    end function CALC_TRIANGLE_AREA

    function CALC_BELONGING_ELEM(x, y) result(ind_elem)
        implicit none
        real, intent(in) :: x, y
        integer :: i
        integer :: ind_elem, ind_node_1, ind_node_2, ind_node_3
        logical :: is_inside
        real :: x_1, y_1, x_2, y_2, x_3, y_3

        ind_elem = -1

        do i = 1, GRID%n_Tri_Elem
            ind_node_1 = GRID%CONN(i, 1)
            ind_node_2 = GRID%CONN(i, 2)
            ind_node_3 = GRID%CONN(i, 3)

            x_1 = GRID%COORD(ind_node_1, 1)
            y_1 = GRID%COORD(ind_node_1, 2)
            x_2 = GRID%COORD(ind_node_2, 1)
            y_2 = GRID%COORD(ind_node_2, 2)
            x_3 = GRID%COORD(ind_node_3, 1)
            y_3 = GRID%COORD(ind_node_3, 2)

            is_inside = IS_POINT_IN_TRIANGLE(x, y, x_1, y_1, x_2, y_2, x_3, y_3)
            if (is_inside) then
                ind_elem = i
                return
            end if
        end do
    end function CALC_BELONGING_ELEM

    logical function IS_POINT_IN_TRIANGLE(x, y, x1, y1, x2, y2, x3, y3)
        implicit none
        real, intent(in) :: x, y
        real, intent(in) :: x1, y1, x2, y2, x3, y3
        real :: area, area1, area2, area3

        area = CALC_TRIANGLE_AREA(x1, y1, x2, y2, x3, y3)
        area1 = CALC_TRIANGLE_AREA(x, y, x2, y2, x3, y3)
        area2 = CALC_TRIANGLE_AREA(x1, y1, x, y, x3, y3)
        area3 = CALC_TRIANGLE_AREA(x1, y1, x2, y2, x, y)

        print *, "Debug: area =", area, ", area1 =", area1, ", area2 =", area2, ", area3 =", area3

        IS_POINT_IN_TRIANGLE = (ABS(area - (area1 + area2 + area3)) < 1.0E-6)
    end function IS_POINT_IN_TRIANGLE

    function INTERPOLATE_STREAMLINE_FUNC(x, y) result(STREAMLINE_FUNC)
        implicit none
        real, intent(in) :: x, y
        real :: STREAMLINE_FUNC
        real :: x_1, y_1, x_2, y_2, x_3, y_3
        integer :: ind_elem, ind_node_1, ind_node_2, ind_node_3
        real :: PSI_1, PSI_2, PSI_3
        real :: N_1, N_2, N_3, AREA

        ind_elem = CALC_BELONGING_ELEM(x, y)
        print *, ind_elem

        if (ind_elem == -1) then
            print *, "Error: Point lies outside the mesh."
            STREAMLINE_FUNC = -10000.0
            return
        end if

        ind_node_1 = GRID%CONN(ind_elem, 1)
        ind_node_2 = GRID%CONN(ind_elem, 2)
        ind_node_3 = GRID%CONN(ind_elem, 3)

        x_1 = GRID%COORD(ind_node_1, 1)
        y_1 = GRID%COORD(ind_node_1, 2)
        x_2 = GRID%COORD(ind_node_2, 1)
        y_2 = GRID%COORD(ind_node_2, 2)
        x_3 = GRID%COORD(ind_node_3, 1)
        y_3 = GRID%COORD(ind_node_3, 2)

        PSI_1 = VALS(ind_node_1)
        PSI_2 = VALS(ind_node_2)
        PSI_3 = VALS(ind_node_3)

        AREA = CALC_TRIANGLE_AREA(x_1, y_1, x_2, y_2, x_3, y_3)

        N_1 = (1 / (2*AREA)) * ((x_2*y_3 - x_3*y_2) + (y_2 - y_3)*x + (x_3 - x_2)*y);
        N_2 = (1 / (2*AREA)) * ((x_3*y_1 - x_1*y_3) + (y_3 - y_1)*x + (x_1 - x_3)*y);
        N_3 = (1 / (2*AREA)) * ((x_1*y_2 - x_2*y_1) + (y_1 - y_2)*x + (x_2 - x_1)*y);

        STREAMLINE_FUNC = PSI_1 * N_1 + PSI_2 * N_2 + PSI_3 * N_3
    end function INTERPOLATE_STREAMLINE_FUNC
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
end module streamline_solution_triangular