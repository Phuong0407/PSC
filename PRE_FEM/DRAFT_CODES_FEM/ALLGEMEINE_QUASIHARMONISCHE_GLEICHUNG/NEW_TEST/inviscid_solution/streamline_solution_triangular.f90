!***************************************************************************
!***************************************************************************
! THIS IS THE PROGRAM TO SOLVE THE LAPLACE EQUATION FOR STREAMLINE FUNCTION
!***************************************************************************
!***************************************************************************
module streamline_solution_triangular
    use DENSE_MATRIX_SOLVER
    implicit none
    private :: INIT_GRID, GEN_GRID_COORD, GEN_ELEM_CONN, GEN_NODE_CONN, &
                ADD_NODE_CONNECTION, &
                CALC_AXISYM_LSM_LOC, CALC_LSM_LOC, CALC_STIFF_MAT_GLO, &
                CALC_DIRICHLET_BOUND, CALC_DIRICHLET_BOUND_AXISYM, &
                MODIFY_STIFF_MAT_AND_FORCE_MAT, &
                CALC_BELONGING_ELEM, IS_POINT_IN_TRIANGLE, CALC_TRIANGLE_AREA, &
                INIT_VELOCITY, INIT_VELOCITY_AXISYM, &
                CALC_2_CORNER_NODE_GRAD_RECOVERY, &
                CALC_3_CORNER_NODE_GRAD_RECOVERY, &
                CALC_BOUND_NODE_GRAD_RECOVERY, &
                CALC_INTERNAL_NODE_GRAD_RECOVERY
    real(8), parameter :: PI = 3.141592653589793
!***************************************************************************
!***************************************************************************
! INPUT DATA
    type :: grid_struct
        real(8) :: H_1, H_2, H_3                        ! Horizontal dimension of each region
        real(8) :: L_1, L_2                             ! Vertical dimension of each region
        integer :: eH_1, eH_2, eH_3                     ! Horizontal NUMBER OF ELEMENT for each region
        integer :: indH_1, indH_2, indH_3               ! Horizontal INDEX OF NODE for each region
        integer :: eV                                   ! Total vertical number of element
        integer :: eH                                   ! Total number of horizontal element
        integer :: nElem                                ! Total number of element
        integer :: n_Tri_Elem                           ! Total number of TRIANGLULAR element
        integer :: nV                                   ! Total number of vertical nodes
        integer :: nH                                   ! Total number of horizontal nodes
        integer :: nNode                                ! Total number of nodes
        real(8), allocatable :: COORD(:,:)              ! Coordinate matrices, NODE_IND[X, Y]
        integer, allocatable :: ELEM_CONN(:,:)          ! Connection Data for Assemble: ELEMENT[NODE1_IND, NODE_2_IND, NODE_3_IND, NODE_4_IND]
        integer, allocatable :: IND_DIRICHLET_BOUND(:)  ! Connection Data for Assemble: ELEMENT[NODE1_IND, NODE_2_IND, NODE_3_IND, NODE_4_IND]
        integer, allocatable :: EDGE_CONN(:,:)          ! Connection Data for Element Boundary Flux: BOUND[NODE1_IND, NODE_2_IND]
        integer, allocatable :: NODE_CONN(:,:)          ! Connection Data for Element Boundary Flux: BOUND[NODE1_IND, NODE_2_IND]
        real(8) :: dx1, dx2, dx3                        ! Finite-Different horizontal dimension of each region
        real(8) :: dy1, dy2                             ! Finite-Different vertical dimension of each region
    contains
    end type grid_struct
    real(8) :: U_0, U_1                               ! Velocity at the Entrance and the Outlet of the Inlet
!***************************************************************************
!***************************************************************************
    public :: INIT_SOLUTION, INTERPOLATE_STREAMLINE_FUNC, &
            CALC_DEV_X, CALC_DEV_Y, CALC_NODE_GRAD_RECOVERY
    type(grid_struct) :: GRID
    real(8), allocatable :: STIFF_MAT_GLO(:,:)
    real(8), allocatable :: DIRICHLET_BOUND(:)
    real(8), allocatable :: FORCE_MAT_GLO(:)
    real(8), allocatable :: NODAL_VALS(:)
    ! real(8), allocatable :: EDGE_DEV(:)
    real(8), allocatable :: NODE_GRAD_RECOVERY(:,:)
contains
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
! GRID DATA GENERATION
    subroutine INIT_GRID(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3)
        implicit none
        real(8), intent(in) :: iH_1, iH_2, iH_3, iL_1, iL_2
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

        GRID%dx1 = GRID%H_1 / GRID%eH_1
        GRID%dx2 = GRID%H_1 / GRID%eH_2
        GRID%dx3 = GRID%H_1 / GRID%eH_3
        GRID%dy1 = GRID%L_1 / GRID%eV
        GRID%dy2 = GRID%L_2 / GRID%eV

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
        real(8) :: x, y, a, b

        allocate(GRID%COORD(GRID%nNode, 2))
        do i = 1, GRID%nV
            do j = 1, GRID%indH_1
                ind_node = (i - 1) * GRID%nH + j
                x = GRID%dx1 * (j - 1)
                y = GRID%dy1 * (i - 1)
                GRID%COORD(ind_node, 1) = x
                GRID%COORD(ind_node, 2) = y
            end do

            a = ((GRID%dy2 - GRID%dy1) * (i - 1) + GRID%L_1 - GRID%L_2)/GRID%H_2
            b = GRID%dy1 * (i - 1) - a * GRID%H_1
            do j = GRID%indH_1 + 1, GRID%indH_2
                ind_node = (i - 1) * GRID%nH + j
                x = GRID%H_1 + GRID%dx2 * (j - GRID%indH_1)
                y = a * x + b
                GRID%COORD(ind_node, 1) = x
                GRID%COORD(ind_node, 2) = y
            end do

            a = GRID%H_1 + GRID%H_2
            b = GRID%L_1 - GRID%L_2
            do j = GRID%indH_2 + 1, GRID%indH_3
                ind_node = (i - 1) * GRID%nH + j
                x = a + GRID%dx3 * (j - 1 - GRID%indH_2)
                y = b + GRID%dy2 * (i - 1)
                GRID%COORD(ind_node, 1) = x
                GRID%COORD(ind_node, 2) = y
            end do
        end do
    end subroutine GEN_GRID_COORD

    subroutine GEN_ELEM_CONN()
        implicit none
        integer :: i, j, ind_elem
        integer :: IND_FIRST_NODE, IND_SECOND_NODE, IND_THIRD_NODE, IND_FOURTH_NODE
        integer :: TRI_ELEM_IND

        allocate(GRID%ELEM_CONN(GRID%n_Tri_Elem, 3))

        TRI_ELEM_IND = 1

        do i = 1, GRID%eV
            do j = 1, GRID%eH
                ind_elem = (i - 1) * GRID%eH + j

                IND_FIRST_NODE = GRID%nH * i + j - GRID%nH
                IND_SECOND_NODE = IND_FIRST_NODE + 1
                IND_THIRD_NODE = IND_SECOND_NODE + GRID%nH
                IND_FOURTH_NODE = IND_FIRST_NODE + GRID%nH

                GRID%ELEM_CONN(TRI_ELEM_IND, 1) = IND_FIRST_NODE
                GRID%ELEM_CONN(TRI_ELEM_IND, 2) = IND_SECOND_NODE
                GRID%ELEM_CONN(TRI_ELEM_IND, 3) = IND_THIRD_NODE
                TRI_ELEM_IND = TRI_ELEM_IND + 1

                GRID%ELEM_CONN(TRI_ELEM_IND, 1) = IND_FIRST_NODE
                GRID%ELEM_CONN(TRI_ELEM_IND, 2) = IND_THIRD_NODE
                GRID%ELEM_CONN(TRI_ELEM_IND, 3) = IND_FOURTH_NODE
                TRI_ELEM_IND = TRI_ELEM_IND + 1
            end do
        end do
        if (TRI_ELEM_IND /= GRID%n_Tri_Elem + 1) then
            print *, "FAIL TO INITIALIZE CONNECTION MATRIX"
            stop
        end if
    end subroutine GEN_ELEM_CONN

    subroutine GEN_NODE_CONN()
        implicit none
        integer :: i, NUM_ELEM, NUM_NODE, NUM_NODE_ELEM
        integer :: IND_FIRST_NODE, IND_SECOND_NODE, IND_THIRD_NODE
        integer, allocatable :: CURRENT_NUM_CONN(:)
        integer :: MAX_CONN_PErNODE
    
        NUM_NODE_ELEM = 3
        NUM_NODE = GRID%nNode
        NUM_ELEM = GRID%n_Tri_Elem
        MAX_CONN_PErNODE = 8
    
        allocate(CURRENT_NUM_CONN(NUM_NODE))
        allocate(GRID%NODE_CONN(NUM_NODE, MAX_CONN_PErNODE))
        CURRENT_NUM_CONN = 1
        GRID%NODE_CONN = 0
    
        do i = 1, NUM_NODE
            GRID%NODE_CONN(i, 1) = i
            CURRENT_NUM_CONN(i) = 2
        end do
    
        do i = 1, NUM_ELEM
            IND_FIRST_NODE = GRID%ELEM_CONN(i, 1)
            IND_SECOND_NODE = GRID%ELEM_CONN(i, 2)
            IND_THIRD_NODE = GRID%ELEM_CONN(i, 3)
    
            call ADD_NODE_CONNECTION(IND_FIRST_NODE, IND_SECOND_NODE, CURRENT_NUM_CONN)
            call ADD_NODE_CONNECTION(IND_FIRST_NODE, IND_THIRD_NODE, CURRENT_NUM_CONN)
    
            call ADD_NODE_CONNECTION(IND_SECOND_NODE, IND_FIRST_NODE, CURRENT_NUM_CONN)
            call ADD_NODE_CONNECTION(IND_SECOND_NODE, IND_THIRD_NODE, CURRENT_NUM_CONN)
    
            call ADD_NODE_CONNECTION(IND_THIRD_NODE, IND_FIRST_NODE, CURRENT_NUM_CONN)
            call ADD_NODE_CONNECTION(IND_THIRD_NODE, IND_SECOND_NODE, CURRENT_NUM_CONN)
        end do

        do i = 1, NUM_NODE
            GRID%NODE_CONN(i, 8) = CURRENT_NUM_CONN(i) - 1
        end do
    
        deallocate(CURRENT_NUM_CONN)
    end subroutine GEN_NODE_CONN
    
    subroutine ADD_NODE_CONNECTION(NODE_FROM, NODE_TO, CURRENT_NUM_CONN)
        implicit none
        integer, intent(in) :: NODE_FROM, NODE_TO
        integer, intent(inout) :: CURRENT_NUM_CONN(:)
        integer :: k, EXISTS
    
        EXISTS = 0
        do k = 2, CURRENT_NUM_CONN(NODE_FROM) - 1
            if (GRID%NODE_CONN(NODE_FROM, k) .eq. NODE_TO) then
                EXISTS = 1
                exit
            end if
        end do
    
        if (EXISTS .eq. 0) then
            GRID%NODE_CONN(NODE_FROM, CURRENT_NUM_CONN(NODE_FROM)) = NODE_TO
            CURRENT_NUM_CONN(NODE_FROM) = CURRENT_NUM_CONN(NODE_FROM) + 1
        end if
    end subroutine ADD_NODE_CONNECTION
    
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
! VELOCITY INITIALIZATION
    subroutine INIT_VELOCITY_AXISYM(U)
        implicit none
        real(8), intent(in) :: U
        U_0 = U
        U_1 = U_0 * (GRID%L_1**2.0)/(GRID%L_1**2.0 - (GRID%L_1 - GRID%L_2)**2.0)
    end subroutine INIT_VELOCITY_AXISYM

    subroutine INIT_VELOCITY(U)
        implicit none
        real(8), intent(in) :: U
        U_0 = U
        U_1 = U_0 * GRID%L_1/GRID%L_2
    end subroutine INIT_VELOCITY

!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
! LOCAL STIFFNESS MATRIX CALCULATION
    function CALC_AXISYM_LSM_LOC(z1, r1, z2, r2, z3, r3) result(STIFF_MAT_LOC)
        implicit none
        real(8), intent(in) :: z1, r1, z2, r2, z3, r3
        real(8), allocatable :: STIFF_MAT_LOC(:,:)
        real(8) :: c(3), d(3)
        real(8) :: A_elem, rcentroid
        integer :: i, j

        c(1) = z2 - z3
        c(2) = z3 - z1
        c(3) = z1 - z2
        d(1) = r3 - r1
        d(2) = r1 - r3
        d(3) = r2 - r1

        A_elem = CALC_TRIANGLE_AREA(z1, r1, z2, r2, z3, r3)
        rcentroid = (r1 + r2 + r3) / 3.0
        allocate(STIFF_MAT_LOC(3, 3))

        do i = 1, 3
            do j = 1, 3
                STIFF_MAT_LOC(i, j) = PI * rcentroid / (2.0 * A_elem) * (c(i) * c(j) + d(i) * d(j))
            end do
        end do
    end function CALC_AXISYM_LSM_LOC

    function CALC_LSM_LOC(x1, y1, x2, y2, x3, y3) result(STIFF_MAT_LOC)
        implicit none
        real(8), intent(in) :: x1, y1, x2, y2, x3, y3
        real(8), allocatable :: STIFF_MAT_LOC(:,:)
        real(8), dimension(3) :: b, c
        real(8) :: AREA
        integer :: i, j

        AREA = CALC_TRIANGLE_AREA(x1, y1, x2, y2, x3, y3)
        b(1) = y2 - y3
        b(2) = y3 - y1
        b(3) = y1 - y2
        c(1) = x3 - x2
        c(2) = x1 - x3
        c(3) = x2 - x1
        
        allocate(STIFF_MAT_LOC(3,3))
        do i = 1, 3
            do j = 1, 3
                STIFF_MAT_LOC(i,j) = (b(i) * b(j) + c(i) * c(j))/(4.0 * AREA)
            end do
        end do
    end function CALC_LSM_LOC
!***************************************************************************
!***************************************************************************
! GLOBAL SYSTEM OF EQUATIONS
    subroutine CALC_STIFF_MAT_GLO()
        implicit none
        integer :: i, LOC_I, LOC_J, GLO_I, GLO_J
        ! integer :: j, k
        real(8), allocatable :: STIFF_MAT_LOC(:,:)
        integer :: ind_node_1, ind_node_2, ind_node_3, ind_node_4
        ! real(8) :: x1, y1, x2, y2, x3, y3
        real(8) :: z1, r1, z2, r2, z3, r3

        allocate(STIFF_MAT_LOC(3,3))
        allocate(STIFF_MAT_GLO(GRID%nNode, GRID%nNode))

        STIFF_MAT_GLO = 0.0
        STIFF_MAT_LOC = 0.0

        do i = 1, GRID%n_Tri_Elem
            ind_node_1 = GRID%ELEM_CONN(i, 1)
            ind_node_2 = GRID%ELEM_CONN(i, 2)
            ind_node_3 = GRID%ELEM_CONN(i, 3)

            ! x1 = GRID%COORD(ind_node_1, 1)
            ! y1 = GRID%COORD(ind_node_1, 2)
            ! x2 = GRID%COORD(ind_node_2, 1)
            ! y2 = GRID%COORD(ind_node_2, 2)
            ! x3 = GRID%COORD(ind_node_3, 1)
            ! y3 = GRID%COORD(ind_node_3, 2)
            ! STIFF_MAT_LOC = CALC_LSM_LOC(x1, y1, x2, y2, x3, y3)
            
            z1 = GRID%COORD(ind_node_1, 1)
            r1 = GRID%COORD(ind_node_1, 2)
            z2 = GRID%COORD(ind_node_2, 1)
            r2 = GRID%COORD(ind_node_2, 2)
            z3 = GRID%COORD(ind_node_3, 1)
            r3 = GRID%COORD(ind_node_3, 2)

            STIFF_MAT_LOC = CALC_AXISYM_LSM_LOC(z1, r1, z2, r2, z3, r3)

            do LOC_I = 1, 3
                GLO_I = GRID%ELEM_CONN(i, LOC_I)
                do LOC_J = 1, 3
                    GLO_J = GRID%ELEM_CONN(i, LOC_J)
                    STIFF_MAT_GLO(GLO_I, GLO_J) = STIFF_MAT_GLO(GLO_I, GLO_J) + STIFF_MAT_LOC(LOC_I, LOC_J)
                end do
            end do
        end do
    end subroutine CALC_STIFF_MAT_GLO
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
! CALCULATE DIRICHLET BOUNDARY CONDITION
    subroutine CALC_DIRICHLET_BOUND()
        implicit none
        integer :: i, k, IND_LOWErBOUND, IND_UPPErBOUND, IND_LEFT_BOUND, IND_RIGHT_BOUND
        real(8) :: PSIUPPER, PSILEFT, PSIRIGHT, r1, r2, r3

        allocate(DIRICHLET_BOUND(GRID%nNode))
        allocate(GRID%IND_DIRICHLET_BOUND(2 * GRID%nH + 2 * GRID%nV - 4))

        DIRICHLET_BOUND = 0.0
        k = 1
        do i = 1, GRID%nH
            IND_LOWErBOUND = i
            GRID%IND_DIRICHLET_BOUND(k) = IND_LOWErBOUND
            k = k + 1
        end do

        r3 = GRID%L_1 - GRID%L_2
        do i = 2, GRID%eV
            IND_LEFT_BOUND = (i - 1) * GRID%nH + 1
            r1 = GRID%dy1 * (i - 1)
            PSILEFT = U_0 * r1
            DIRICHLET_BOUND(IND_LEFT_BOUND) = PSILEFT
            GRID%IND_DIRICHLET_BOUND(k) = IND_LEFT_BOUND
            k = k + 1
            
            IND_RIGHT_BOUND = i * GRID%nH
            r2 = GRID%dy2 * (i - 1) + r3
            PSIRIGHT = U_1 * (r2 - r3)
            DIRICHLET_BOUND(IND_RIGHT_BOUND) = PSIRIGHT
            GRID%IND_DIRICHLET_BOUND(k) = IND_RIGHT_BOUND
            k = k + 1
        end do

        PSIUPPER = U_0 * GRID%L_1
        do i = 1, GRID%nH
            IND_UPPErBOUND = (GRID%nV - 1) * GRID%nH + i
            DIRICHLET_BOUND(IND_UPPErBOUND) = PSIUPPER
            GRID%IND_DIRICHLET_BOUND(k) = IND_UPPErBOUND
            k = k + 1
        end do
        if (k /= size(GRID%IND_DIRICHLET_BOUND) + 1) then
            print *, "FAIL TO INITIALIZE IND_DIRICHLET_BOUND"
            stop
        end if
    end subroutine CALC_DIRICHLET_BOUND

    subroutine CALC_DIRICHLET_BOUND_AXISYM()
        implicit none
        integer :: i, k, IND_LOWErBOUND, IND_UPPErBOUND, IND_LEFT_BOUND, IND_RIGHT_BOUND
        real(8) :: PSIUPPER, PSILEFT, PSIRIGHT, r1, r2, r3

        allocate(DIRICHLET_BOUND(GRID%nNode))
        allocate(GRID%IND_DIRICHLET_BOUND(2 * GRID%nH + 2 * GRID%nV - 4))

        DIRICHLET_BOUND = 0.0
        k = 1
        do i = 1, GRID%nH
            IND_LOWErBOUND = i
            GRID%IND_DIRICHLET_BOUND(k) = IND_LOWErBOUND
            k = k + 1
        end do

        r3 = GRID%L_1 - GRID%L_2
        do i = 2, GRID%eV
            IND_LEFT_BOUND = (i - 1) * GRID%nH + 1
            r1 = GRID%dy1 * (i - 1)
            PSILEFT = 0.5 * U_0 * r1**2.0
            DIRICHLET_BOUND(IND_LEFT_BOUND) = PSILEFT
            GRID%IND_DIRICHLET_BOUND(k) = IND_LEFT_BOUND
            k = k + 1
            
            IND_RIGHT_BOUND = i * GRID%nH
            r2 = GRID%dy2 * (i - 1) + r3
            PSIRIGHT = 0.5 * U_1 * (r2**2.0 - r3**2.0)
            DIRICHLET_BOUND(IND_RIGHT_BOUND) = PSIRIGHT
            GRID%IND_DIRICHLET_BOUND(k) = IND_RIGHT_BOUND
            k = k + 1
        end do

        PSIUPPER = 0.5 * U_0 * (GRID%L_1)**2.0
        do i = 1, GRID%nH
            IND_UPPErBOUND = (GRID%nV - 1) * GRID%nH + i
            DIRICHLET_BOUND(IND_UPPErBOUND) = PSIUPPER
            GRID%IND_DIRICHLET_BOUND(k) = IND_UPPErBOUND
            k = k + 1
        end do
        if (k /= size(GRID%IND_DIRICHLET_BOUND) + 1) then
            print *, "FAIL TO INITIALIZE IND_DIRICHLET_BOUND"
            stop
        end if
    end subroutine CALC_DIRICHLET_BOUND_AXISYM
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
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
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
    subroutine INIT_SOLUTION(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3, iU_0)
        implicit none
        real(8), intent(in) :: iH_1, iH_2, iH_3, iL_1, iL_2, iU_0
        integer, intent(in) :: iN, iM_1, iM_2, iM_3
        ! integer :: i, j
        call INIT_GRID(iH_1, iH_2, iH_3, iL_1, iL_2, iN, iM_1, iM_2, iM_3)
        call GEN_GRID_COORD()
        call GEN_ELEM_CONN()
        call GEN_NODE_CONN()
        ! call GEN_EDGE_CONN()
        ! call INIT_VELOCITY(iU_0)
        call INIT_VELOCITY_AXISYM(iU_0)
        call CALC_STIFF_MAT_GLO()
        ! call CALC_DIRICHLET_BOUND()
        call CALC_DIRICHLET_BOUND_AXISYM()
        call MODIFY_STIFF_MAT_AND_FORCE_MAT()
    end subroutine INIT_SOLUTION
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
! STREAMLINE FUNCTION INTERPOLATION
    function CALC_TRIANGLE_AREA(x1, y1, x2, y2, x3, y3) result(triangle_area)
        implicit none
        real(8), intent(in) :: x1, y1, x2, y2, x3, y3
        real(8) :: triangle_area

        triangle_area = 0.5 * ABS((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)))
    end function CALC_TRIANGLE_AREA

    function CALC_BELONGING_ELEM(x, y) result(ind_elem)
        implicit none
        real(8), intent(in) :: x, y
        integer :: i
        integer :: ind_elem, ind_node_1, ind_node_2, ind_node_3
        logical :: is_inside
        real(8) :: x1, y1, x2, y2, x3, y3

        ind_elem = -1

        do i = 1, GRID%n_Tri_Elem
            ind_node_1 = GRID%ELEM_CONN(i, 1)
            ind_node_2 = GRID%ELEM_CONN(i, 2)
            ind_node_3 = GRID%ELEM_CONN(i, 3)

            x1 = GRID%COORD(ind_node_1, 1)
            y1 = GRID%COORD(ind_node_1, 2)
            x2 = GRID%COORD(ind_node_2, 1)
            y2 = GRID%COORD(ind_node_2, 2)
            x3 = GRID%COORD(ind_node_3, 1)
            y3 = GRID%COORD(ind_node_3, 2)

            is_inside = IS_POINT_IN_TRIANGLE(x, y, x1, y1, x2, y2, x3, y3)
            if (is_inside) then
                ind_elem = i
                return
            end if
        end do
    end function CALC_BELONGING_ELEM

    logical function IS_POINT_IN_TRIANGLE(x, y, x1, y1, x2, y2, x3, y3)
        implicit none
        real(8), intent(in) :: x, y
        real(8), intent(in) :: x1, y1, x2, y2, x3, y3
        real(8) :: area, area1, area2, area3

        area = CALC_TRIANGLE_AREA(x1, y1, x2, y2, x3, y3)
        area1 = CALC_TRIANGLE_AREA(x, y, x2, y2, x3, y3)
        area2 = CALC_TRIANGLE_AREA(x1, y1, x, y, x3, y3)
        area3 = CALC_TRIANGLE_AREA(x1, y1, x2, y2, x, y)

        IS_POINT_IN_TRIANGLE = (ABS(area - (area1 + area2 + area3)) < 1.0E-6)
    end function IS_POINT_IN_TRIANGLE

    function INTERPOLATE_STREAMLINE_FUNC(x, y) result(STREAMLINE_FUNC)
        implicit none
        real(8), intent(in) :: x, y
        real(8) :: STREAMLINE_FUNC
        real(8) :: x1, y1, x2, y2, x3, y3
        integer :: ind_elem, ind_node_1, ind_node_2, ind_node_3
        real(8) :: PSI1, PSI2, PSI3
        real(8) :: N_1, N_2, N_3, AREA

        ind_elem = CALC_BELONGING_ELEM(x, y)

        if (ind_elem .eq. -1) then
            print *, "Error: Point lies outside the mesh."
            STREAMLINE_FUNC = -10000.0
            return
        end if

        ind_node_1 = GRID%ELEM_CONN(ind_elem, 1)
        ind_node_2 = GRID%ELEM_CONN(ind_elem, 2)
        ind_node_3 = GRID%ELEM_CONN(ind_elem, 3)

        x1 = GRID%COORD(ind_node_1, 1)
        y1 = GRID%COORD(ind_node_1, 2)
        x2 = GRID%COORD(ind_node_2, 1)
        y2 = GRID%COORD(ind_node_2, 2)
        x3 = GRID%COORD(ind_node_3, 1)
        y3 = GRID%COORD(ind_node_3, 2)

        PSI1 = NODAL_VALS(ind_node_1)
        PSI2 = NODAL_VALS(ind_node_2)
        PSI3 = NODAL_VALS(ind_node_3)

        AREA = CALC_TRIANGLE_AREA(x1, y1, x2, y2, x3, y3)

        N_1 = (1 / (2*AREA)) * ((x2*y3 - x3*y2) + (y2 - y3)*x + (x3 - x2)*y);
        N_2 = (1 / (2*AREA)) * ((x3*y1 - x1*y3) + (y3 - y1)*x + (x1 - x3)*y);
        N_3 = (1 / (2*AREA)) * ((x1*y2 - x2*y1) + (y1 - y2)*x + (x2 - x1)*y);

        STREAMLINE_FUNC = PSI1 * N_1 + PSI2 * N_2 + PSI3 * N_3
    end function INTERPOLATE_STREAMLINE_FUNC
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
! CALCULATE DERIVATIVE OF PSI
    function CALC_DEV_X(x, y) result(xDEV)
        implicit none
        real(8), intent(in) :: x, y
        real(8) :: xDEV
        real(8) :: x1, y1, x2, y2, x3, y3
        integer :: ind_elem, ind_node_1, ind_node_2, ind_node_3
        real(8) :: PSI1, PSI2, PSI3
        real(8) :: AREA

        ind_elem = CALC_BELONGING_ELEM(x, y)

        if (ind_elem .eq. -1) then
            print *, "ERROR: POINT LIES OUTSIDE THE MESH."
            return
        end if

        ind_node_1 = GRID%ELEM_CONN(ind_elem, 1)
        ind_node_2 = GRID%ELEM_CONN(ind_elem, 2)
        ind_node_3 = GRID%ELEM_CONN(ind_elem, 3)

        x1 = GRID%COORD(ind_node_1, 1)
        y1 = GRID%COORD(ind_node_1, 2)
        x2 = GRID%COORD(ind_node_2, 1)
        y2 = GRID%COORD(ind_node_2, 2)
        x3 = GRID%COORD(ind_node_3, 1)
        y3 = GRID%COORD(ind_node_3, 2)

        PSI1 = NODAL_VALS(ind_node_1)
        PSI2 = NODAL_VALS(ind_node_2)
        PSI3 = NODAL_VALS(ind_node_3)
        AREA = CALC_TRIANGLE_AREA(x1, y1, x2, y2, x3, y3)

        xDEV = ((y2 - y3)*PSI1 + (y3 - y1)*PSI2 + (y1 - y2)*PSI3)/(2.0 * AREA)
    end function CALC_DEV_X

    function CALC_DEV_Y(x, y) result(yDEV)
        implicit none
        real(8), intent(in) :: x, y
        real(8) :: yDEV
        real(8) :: x1, y1, x2, y2, x3, y3
        integer :: ind_elem, ind_node_1, ind_node_2, ind_node_3
        real(8) :: PSI1, PSI2, PSI3
        real(8) :: AREA

        ind_elem = CALC_BELONGING_ELEM(x, y)

        if (ind_elem .eq. -1) then
            print *, "ERROR: POINT LIES OUTSIDE THE MESH."
            return
        end if

        ind_node_1 = GRID%ELEM_CONN(ind_elem, 1)
        ind_node_2 = GRID%ELEM_CONN(ind_elem, 2)
        ind_node_3 = GRID%ELEM_CONN(ind_elem, 3)

        x1 = GRID%COORD(ind_node_1, 1)
        y1 = GRID%COORD(ind_node_1, 2)
        x2 = GRID%COORD(ind_node_2, 1)
        y2 = GRID%COORD(ind_node_2, 2)
        x3 = GRID%COORD(ind_node_3, 1)
        y3 = GRID%COORD(ind_node_3, 2)

        PSI1 = NODAL_VALS(ind_node_1)
        PSI2 = NODAL_VALS(ind_node_2)
        PSI3 = NODAL_VALS(ind_node_3)

        AREA = CALC_TRIANGLE_AREA(x1, y1, x2, y2, x3, y3)

        yDEV = ((x3 - x2)*PSI1 + (x1 - x3)*PSI2 + (x2 - x1)*PSI3)/(2.0 * AREA)
    end function CALC_DEV_Y
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
! PATCH RECOVERY OF STREAMLINE FUCNTION GRADIENT
    function SEARCH_FOR_PATCH_ELEM(IND_NODE) result(IND_ELEM)
        implicit none
        integer, intent(in) :: IND_NODE
        integer, allocatable :: IND_ELEM(:) 

        integer :: N_NODE_IN_ELEM
        integer :: IND_NODE_IN_ELEM
        integer :: i, j
        integer :: temp_patch(1000)
        integer :: count

        N_NODE_IN_ELEM = size(GRID%ELEM_CONN, 2) 
        count = 0

        do i = 1, GRID%n_Tri_Elem
            do j = 1, N_NODE_IN_ELEM
                IND_NODE_IN_ELEM = GRID%ELEM_CONN(i, j)
                if (IND_NODE .eq. IND_NODE_IN_ELEM) then
                    count = count + 1
                    temp_patch(count) = i
                    exit
                end if
            end do
        end do

        allocate(IND_ELEM(count))
        IND_ELEM(:) = temp_patch(1:count)
    end function SEARCH_FOR_PATCH_ELEM

    subroutine CALC_NODE_GRAD_RECOVERY()
        implicit none
        integer :: NUM_NODE, i, NUM_NODE_CONN

        NUM_NODE = GRID%nNode
        allocate(NODE_GRAD_RECOVERY(NUM_NODE, 5))

        do i = 1, NUM_NODE
            NUM_NODE_CONN = GRID%NODE_CONN(i, 8)
            ! print *, "NODE ELEMENT", i, "    NUMBER OF CONNECTED NODES", NUM_NODE_CONN
            if (NUM_NODE_CONN .eq. 7) then
                call CALC_INTERNAL_NODE_GRAD_RECOVERY(i)
            else if (NUM_NODE_CONN .eq. 5) then
                call CALC_BOUND_NODE_GRAD_RECOVERY(i)
            else if (NUM_NODE_CONN .eq. 4) then
                call CALC_3_CORNER_NODE_GRAD_RECOVERY(i)
            else if (NUM_NODE_CONN .eq. 3) then
                call CALC_2_CORNER_NODE_GRAD_RECOVERY(i)
            end if
        end do
    end subroutine CALC_NODE_GRAD_RECOVERY

    subroutine CALC_INTERNAL_NODE_GRAD_RECOVERY(IND_INTERNAL_NODE)
        implicit none
        integer, intent(in) :: IND_INTERNAL_NODE
        integer :: IND_CONN_NODE
        real(8) :: LONGEST_LENGTH
        real(8) ::  x_i(7), y_i(7), xi(7), eta(7)
        real(8) :: FITTING_MATRIX(7, 6), NORMED_FITTING_MATRIX(6,6), RHS_FITTING(7), NORMED_RHS_FITTING(6), FITTING_COEFF(6)
        real(8) :: X_DEV, Y_DEV, XX_DEV, XY_DEV, YY_DEV
        integer :: i, DUM

        do i = 1, 7
            IND_CONN_NODE = GRID%NODE_CONN(IND_INTERNAL_NODE, i)
            x_i(i) = GRID%COORD(IND_CONN_NODE, 1)
            y_i(i) = GRID%COORD(IND_CONN_NODE, 2)
            RHS_FITTING(i) = NODAL_VALS(IND_CONN_NODE)
        end do
        
        LONGEST_LENGTH = CALC_LONGEST_LENGTH_OF_NODE_CONNECTED_TO_A_NODE(x_i, y_i, DUM)
        do i = 1, 7
            xi(i) = (x_i(i) - x_i(1))/LONGEST_LENGTH
            eta(i) = (y_i(i) - y_i(1))/LONGEST_LENGTH
        end do

        do i = 1, 7
            FITTING_MATRIX(i, 1) = 1.0
            FITTING_MATRIX(i, 2) = xi(i)
            FITTING_MATRIX(i, 3) = eta(i)
            FITTING_MATRIX(i, 4) = (xi(i))**2.0
            FITTING_MATRIX(i, 5) = xi(i)*eta(i)
            FITTING_MATRIX(i, 6) = (eta(i))**2.0
        end do

        NORMED_FITTING_MATRIX = MATMUL(transpose(FITTING_MATRIX), FITTING_MATRIX)
        NORMED_RHS_FITTING = MATMUL(transpose(FITTING_MATRIX), RHS_FITTING)
        FITTING_COEFF = SOLVE_DENSE_MATRIX_SYSTEM(NORMED_FITTING_MATRIX, NORMED_RHS_FITTING)

        X_DEV = FITTING_COEFF(2)/LONGEST_LENGTH
        Y_DEV = FITTING_COEFF(3)/LONGEST_LENGTH
        XX_DEV = 2.0 * FITTING_COEFF(4)/(LONGEST_LENGTH**2.0)
        XY_DEV = FITTING_COEFF(5)/(LONGEST_LENGTH**2.0)
        YY_DEV = 2.0 * FITTING_COEFF(6)/(LONGEST_LENGTH**2.0)

        NODE_GRAD_RECOVERY(IND_INTERNAL_NODE, 1) = X_DEV
        NODE_GRAD_RECOVERY(IND_INTERNAL_NODE, 2) = Y_DEV
        NODE_GRAD_RECOVERY(IND_INTERNAL_NODE, 3) = XX_DEV
        NODE_GRAD_RECOVERY(IND_INTERNAL_NODE, 4) = XY_DEV
        NODE_GRAD_RECOVERY(IND_INTERNAL_NODE, 5) = YY_DEV
    end subroutine CALC_INTERNAL_NODE_GRAD_RECOVERY

    subroutine CALC_BOUND_NODE_GRAD_RECOVERY(BOUND_NODE_IND)
        implicit none
        integer, intent(in) :: BOUND_NODE_IND
        integer :: IND_INTERNAL_NODE(2)                                         ! 2 IS A SPECIFIC NUMBER AS WE USE STRUCTURED GRID
        integer :: INTERNAL_NODE_COUNT, NUM_FITTING_IND, NUM_FITTING_POINT
        integer, allocatable:: IND_FITTING_NODE(:)
        integer :: i, j, k, DUM_IND
        real(8), allocatable :: x_i(:), y_i(:), RHS_FITTING(:), xi(:), eta(:)
        real(8), allocatable :: FITTING_MATRIX(:,:)
        real(8) :: NORMED_FITTING_MATRIX(6,6), NORMED_RHS_FITTING(6), FITTING_COEFF(6)
        real(8) :: X_DEV, Y_DEV, XX_DEV, XY_DEV, YY_DEV
        real(8) :: LONGEST_LENGTH
        logical :: IS_REPEATED

        IS_REPEATED = .false.
        INTERNAL_NODE_COUNT = 1
        do i = 2, 5
            if (GRID%NODE_CONN(GRID%NODE_CONN(BOUND_NODE_IND, i), 8) .eq. 7) then
                IND_INTERNAL_NODE(INTERNAL_NODE_COUNT) = GRID%NODE_CONN(BOUND_NODE_IND, i)
                INTERNAL_NODE_COUNT = INTERNAL_NODE_COUNT + 1
            end if
        end do
        INTERNAL_NODE_COUNT = INTERNAL_NODE_COUNT - 1

        if (INTERNAL_NODE_COUNT == 1) then
            NUM_FITTING_POINT = 8
        else if (INTERNAL_NODE_COUNT == 2) then
            NUM_FITTING_POINT = 10
        end if

        allocate(x_i(NUM_FITTING_POINT))
        allocate(y_i(NUM_FITTING_POINT))
        allocate(RHS_FITTING(NUM_FITTING_POINT))
        allocate(xi(NUM_FITTING_POINT))
        allocate(eta(NUM_FITTING_POINT))
        allocate(FITTING_MATRIX(NUM_FITTING_POINT,6))
        allocate(IND_FITTING_NODE(NUM_FITTING_POINT))

        IND_FITTING_NODE = 0

        do i = 1, 5
            ! print *, GRID%NODE_CONN(BOUND_NODE_IND, i)
            IND_FITTING_NODE(i) = GRID%NODE_CONN(BOUND_NODE_IND, i)
        end do

        NUM_FITTING_IND = 6
        do i = 1, INTERNAL_NODE_COUNT ! CONNECTED INTERNAL NODE INDICES
            ! print *, IND_INTERNAL_NODE(i), "IND_INTERNAL_NODE(i)"
            do j = 2, 7
                DUM_IND = GRID%NODE_CONN(IND_INTERNAL_NODE(i), j)
                do k = 1, NUM_FITTING_IND
                    IS_REPEATED = .false.
                    ! print *, "k = ", k, " IND_FITTING_NODE(k) = ", IND_FITTING_NODE(k), "DUM_IND = ", DUM_IND, "NUM_FITTING_IND = ", NUM_FITTING_IND
                    if (IND_FITTING_NODE(k) .eq. DUM_IND) then
                        IS_REPEATED = .true.
                        exit
                    end if
                end do
                if (.not. IS_REPEATED) then
                    IND_FITTING_NODE(NUM_FITTING_IND) = DUM_IND
                    NUM_FITTING_IND = NUM_FITTING_IND + 1
                    ! print *, NUM_FITTING_IND, "ASFDASFD", j
                end if
            end do
        end do

        ! print *, IND_FITTING_NODE(:)

        do i = 1, NUM_FITTING_POINT
            x_i(i) = GRID%COORD(IND_FITTING_NODE(i), 1)
            y_i(i) = GRID%COORD(IND_FITTING_NODE(i), 2)
            RHS_FITTING(i) = NODAL_VALS(IND_FITTING_NODE(i))
        end do

        LONGEST_LENGTH = CALC_LONGEST_LENGTH_OF_NODE_CONNECTED_TO_A_NODE(x_i, y_i, DUM_IND)
        ! print *, "longest length", LONGEST_LENGTHprint *

        do i = 1, NUM_FITTING_POINT
            xi(i) = (x_i(i) - x_i(1))/LONGEST_LENGTH
            eta(i) = (y_i(i) - y_i(1))/LONGEST_LENGTH
        end do

        do i = 1, NUM_FITTING_POINT
            FITTING_MATRIX(i, 1) = 1.0
            FITTING_MATRIX(i, 2) = xi(i)
            FITTING_MATRIX(i, 3) = eta(i)
            FITTING_MATRIX(i, 4) = (xi(i))**2.0
            FITTING_MATRIX(i, 5) = xi(i)*eta(i)
            FITTING_MATRIX(i, 6) = (eta(i))**2.0
        end do

        NORMED_FITTING_MATRIX = MATMUL(transpose(FITTING_MATRIX), FITTING_MATRIX)
        NORMED_RHS_FITTING = MATMUL(transpose(FITTING_MATRIX), RHS_FITTING)
        FITTING_COEFF = SOLVE_DENSE_MATRIX_SYSTEM(NORMED_FITTING_MATRIX, NORMED_RHS_FITTING)
        ! print *, FITTING_COEFF(:)
        X_DEV = FITTING_COEFF(2)/LONGEST_LENGTH
        Y_DEV = FITTING_COEFF(3)/LONGEST_LENGTH
        XX_DEV = 2.0 * FITTING_COEFF(4)/(LONGEST_LENGTH**2.0)
        XY_DEV = FITTING_COEFF(5)/(LONGEST_LENGTH**2.0)
        YY_DEV = 2.0 * FITTING_COEFF(6)/(LONGEST_LENGTH**2.0)

        NODE_GRAD_RECOVERY(BOUND_NODE_IND, 1) = X_DEV
        NODE_GRAD_RECOVERY(BOUND_NODE_IND, 2) = Y_DEV
        NODE_GRAD_RECOVERY(BOUND_NODE_IND, 3) = XX_DEV
        NODE_GRAD_RECOVERY(BOUND_NODE_IND, 4) = XY_DEV
        NODE_GRAD_RECOVERY(BOUND_NODE_IND, 5) = YY_DEV
    end subroutine CALC_BOUND_NODE_GRAD_RECOVERY

    subroutine CALC_3_CORNER_NODE_GRAD_RECOVERY(CORNEr3_IND)
        implicit none
        integer, intent(in) :: CORNEr3_IND
        integer :: IND_INTERNAL_NODE                                                                ! 2 IS A SPECIFIC NUMBER AS WE USE STRUCTURED GRID
        integer :: INTERNAL_NODE_COUNT, NUM_FITTING_IND, IND_FITTING_NODE(7)                        ! 7 IS A SPECIFIC NUMBER AS WE USE STRUCTURED GRID
        integer :: i, j, DUM_IND
        real(8) :: x_i(7), y_i(7), RHS_FITTING(7), xi(7), eta(7)
        real(8) :: FITTING_MATRIX(7, 6), NORMED_FITTING_MATRIX(6,6), NORMED_RHS_FITTING(6), FITTING_COEFF(6)
        real(8) :: X_DEV, Y_DEV, XX_DEV, XY_DEV, YY_DEV
        real(8) :: LONGEST_LENGTH
        logical :: IS_REPEATED
        
        IND_FITTING_NODE = 0
        IND_FITTING_NODE(1) = CORNEr3_IND
        do i = 2, 4
            if (GRID%NODE_CONN(GRID%NODE_CONN(CORNEr3_IND, i), 8) == 7) then
                IND_INTERNAL_NODE = GRID%NODE_CONN(CORNEr3_IND, i)
                INTERNAL_NODE_COUNT = INTERNAL_NODE_COUNT + 1
            end if
            IND_FITTING_NODE(i) = GRID%NODE_CONN(CORNEr3_IND, i)
        end do
        do i = 2, 7 ! CONNECTED INTERNAL NODE INDICES
            DUM_IND = GRID%NODE_CONN(IND_INTERNAL_NODE, i)
            do j = 1, size(IND_FITTING_NODE)
                if (IND_FITTING_NODE(j) .eq. DUM_IND) then
                    IS_REPEATED = .true.
                    exit
                end if
            end do
            if (.not. IS_REPEATED) then
                IND_FITTING_NODE(NUM_FITTING_IND) = DUM_IND
                NUM_FITTING_IND = NUM_FITTING_IND + 1
            end if
        end do

        do i = 1, 7
            x_i(i) = GRID%COORD(IND_FITTING_NODE(i), 1)
            y_i(i) = GRID%COORD(IND_FITTING_NODE(i), 2)
            RHS_FITTING(i) = NODAL_VALS(IND_FITTING_NODE(i))
        end do
        LONGEST_LENGTH = CALC_LONGEST_LENGTH_OF_NODE_CONNECTED_TO_A_NODE(x_i, y_i, DUM_IND)

        do i = 1, 7
            xi(i) = (x_i(i) - x_i(1))/LONGEST_LENGTH
            eta(i) = (y_i(i) - y_i(1))/LONGEST_LENGTH
        end do

        do i = 1, 7
            FITTING_MATRIX(i, 1) = 1.0
            FITTING_MATRIX(i, 2) = xi(i)
            FITTING_MATRIX(i, 3) = eta(i)
            FITTING_MATRIX(i, 4) = (xi(i))**2.0
            FITTING_MATRIX(i, 5) = xi(i)*eta(i)
            FITTING_MATRIX(i, 6) = (eta(i))**2.0
        end do

        NORMED_FITTING_MATRIX = MATMUL(transpose(FITTING_MATRIX), FITTING_MATRIX)
        NORMED_RHS_FITTING = MATMUL(transpose(FITTING_MATRIX), RHS_FITTING)
        FITTING_COEFF = SOLVE_DENSE_MATRIX_SYSTEM(NORMED_FITTING_MATRIX, NORMED_RHS_FITTING)

        X_DEV = FITTING_COEFF(2)/LONGEST_LENGTH
        Y_DEV = FITTING_COEFF(3)/LONGEST_LENGTH
        XX_DEV = 2.0 * FITTING_COEFF(4)/(LONGEST_LENGTH**2.0)
        XY_DEV = FITTING_COEFF(5)/(LONGEST_LENGTH**2.0)
        YY_DEV = 2.0 * FITTING_COEFF(6)/(LONGEST_LENGTH**2.0)

        NODE_GRAD_RECOVERY(CORNEr3_IND, 1) = X_DEV
        NODE_GRAD_RECOVERY(CORNEr3_IND, 2) = Y_DEV
        NODE_GRAD_RECOVERY(CORNEr3_IND, 3) = XX_DEV
        NODE_GRAD_RECOVERY(CORNEr3_IND, 4) = XY_DEV
        NODE_GRAD_RECOVERY(CORNEr3_IND, 5) = YY_DEV
    end subroutine CALC_3_CORNER_NODE_GRAD_RECOVERY

    subroutine CALC_2_CORNER_NODE_GRAD_RECOVERY(CORNER_2_IND)
        implicit none
        integer, intent(in) :: CORNER_2_IND
        integer :: IND_INTERNAL_NODE                                                                ! 2 IS A SPECIFIC NUMBER AS WE USE STRUCTURED GRID
        integer :: INTERNAL_NODE_COUNT, NUM_FITTING_IND, IND_FITTING_NODE(8)                        ! 8 IS A SPECIFIC NUMBER AS WE USE STRUCTURED GRID
        integer :: i, j, DUM_IND
        real(8) :: x_i(8), y_i(8), RHS_FITTING(8), xi(8), eta(8)
        real(8) :: FITTING_MATRIX(8, 6), NORMED_FITTING_MATRIX(6,6), NORMED_RHS_FITTING(6), FITTING_COEFF(6)
        real(8) :: X_DEV, Y_DEV, XX_DEV, XY_DEV, YY_DEV
        real(8) :: LONGEST_LENGTH
        logical :: IS_REPEATED
        integer :: IND_AUX

        IND_FITTING_NODE = 0
        IND_FITTING_NODE(1) = CORNER_2_IND
        IND_AUX = GRID%NODE_CONN(CORNER_2_IND, 2)
        do i = 2, 4
            if (GRID%NODE_CONN(GRID%NODE_CONN(IND_AUX, i), 8) == 7) then
                IND_INTERNAL_NODE = GRID%NODE_CONN(IND_AUX, i)
                INTERNAL_NODE_COUNT = INTERNAL_NODE_COUNT + 1
            end if
            IND_FITTING_NODE(i) = GRID%NODE_CONN(CORNER_2_IND, i)
        end do
        NUM_FITTING_IND = 3
        do i = 2, 7 ! CONNECTED INTERNAL NODE INDICES
            DUM_IND = GRID%NODE_CONN(IND_INTERNAL_NODE, i)
            do j = 2, NUM_FITTING_IND
                if (IND_FITTING_NODE(j) .eq. DUM_IND) then
                    IS_REPEATED = .true.
                    exit
                end if
            end do
            if (.not. IS_REPEATED) then
                IND_FITTING_NODE(NUM_FITTING_IND) = DUM_IND
                NUM_FITTING_IND = NUM_FITTING_IND + 1
            end if
        end do

        do i = 1, 8
            x_i(i) = GRID%COORD(IND_FITTING_NODE(i), 1)
            y_i(i) = GRID%COORD(IND_FITTING_NODE(i), 2)
            RHS_FITTING(i) = NODAL_VALS(IND_FITTING_NODE(i))
        end do
        LONGEST_LENGTH = CALC_LONGEST_LENGTH_OF_NODE_CONNECTED_TO_A_NODE(x_i, y_i, DUM_IND)

        do i = 1, 8
            xi(i) = (x_i(i) - x_i(1))/LONGEST_LENGTH
            eta(i) = (y_i(i) - y_i(1))/LONGEST_LENGTH
        end do

        do i = 1, 8
            FITTING_MATRIX(i, 1) = 1.0
            FITTING_MATRIX(i, 2) = xi(i)
            FITTING_MATRIX(i, 3) = eta(i)
            FITTING_MATRIX(i, 4) = (xi(i))**2.0
            FITTING_MATRIX(i, 5) = xi(i)*eta(i)
            FITTING_MATRIX(i, 6) = (eta(i))**2.0
        end do

        NORMED_FITTING_MATRIX = MATMUL(transpose(FITTING_MATRIX), FITTING_MATRIX)
        NORMED_RHS_FITTING = MATMUL(transpose(FITTING_MATRIX), RHS_FITTING)
        FITTING_COEFF = SOLVE_DENSE_MATRIX_SYSTEM(NORMED_FITTING_MATRIX, NORMED_RHS_FITTING)

        X_DEV = FITTING_COEFF(2)/LONGEST_LENGTH
        Y_DEV = FITTING_COEFF(3)/LONGEST_LENGTH
        XX_DEV = 2.0 * FITTING_COEFF(4)/(LONGEST_LENGTH**2.0)
        XY_DEV = FITTING_COEFF(5)/(LONGEST_LENGTH**2.0)
        YY_DEV = 2.0 * FITTING_COEFF(6)/(LONGEST_LENGTH**2.0)

        NODE_GRAD_RECOVERY(CORNER_2_IND, 1) = X_DEV
        NODE_GRAD_RECOVERY(CORNER_2_IND, 2) = Y_DEV
        NODE_GRAD_RECOVERY(CORNER_2_IND, 3) = XX_DEV
        NODE_GRAD_RECOVERY(CORNER_2_IND, 4) = XY_DEV
        NODE_GRAD_RECOVERY(CORNER_2_IND, 5) = YY_DEV
    end subroutine CALC_2_CORNER_NODE_GRAD_RECOVERY

    function CALC_LONGEST_LENGTH_OF_NODE_CONNECTED_TO_A_NODE(x_i, y_i, IND) result(MAX_DIST)
        implicit none
        real(8), intent(in) :: x_i(:), y_i(:)
        integer, intent(out) :: IND
        integer :: i, NUMBErOF_ADJACENT_POINTS
        real(8) :: DIST, MAX_DIST

        NUMBErOF_ADJACENT_POINTS = size(x_i) - 1
        MAX_DIST = sqrt((x_i(2) - x_i(1))**2.0 + (y_i(2) - y_i(1))**2.0)
        IND = 2

        do i = 3, NUMBErOF_ADJACENT_POINTS
            DIST = sqrt((x_i(i) - x_i(1))**2.0 + (y_i(i) - y_i(1))**2.0)
            if (DIST > MAX_DIST) then
                MAX_DIST = DIST
                IND = i
            end if
        end do
    end function CALC_LONGEST_LENGTH_OF_NODE_CONNECTED_TO_A_NODE
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
    function INTERPOLATE_DERIVATIVE_AXISYM(z, r) result(INTERPOLATED_DERIVATIVE)
        implicit none
        real(8), intent(in) :: r, z                            ! Coordinates of the point
        real(8) :: INTERPOLATED_DERIVATIVE(5)                 ! Output interpolated [DX, DY, DXX, DXY, DYY]
    
        integer :: ind_elem                                    ! Index of the containing triangle
        integer :: ind_node_1, ind_node_2, ind_node_3          ! Nodes of the triangle
        real(8) :: r1, z1, r2, z2, r3, z3                      ! Coordinates of the triangle nodes
        real(8) :: N_1, N_2, N_3                               ! Shape functions
        real(8) :: DX1, DY1, DXX1, DXY1, DYY1                  ! Nodal derivatives at node 1
        real(8) :: DX2, DY2, DXX2, DXY2, DYY2                  ! Nodal derivatives at node 2
        real(8) :: DX3, DY3, DXX3, DXY3, DYY3                  ! Nodal derivatives at node 3
        real(8) :: b1, b2, b3, c1, c2, c3, d1, d2, d3          ! Coefficients for shape functions
        real(8) :: AREA                                        ! Area of the triangle
    
        ! Step 1: Identify the triangle containing the point (r, z)
        ind_elem = CALC_BELONGING_ELEM(z, r)
        if (ind_elem .eq. -1) then
            print *, "ERROR: POINT LIES OUTSIDE THE MESH."
            INTERPOLATED_DERIVATIVE = -10000.0
            return
        end if
    
        ! Step 2: Get the indices and coordinates of the triangle's nodes
        ind_node_1 = GRID%ELEM_CONN(ind_elem, 1)
        ind_node_2 = GRID%ELEM_CONN(ind_elem, 2)
        ind_node_3 = GRID%ELEM_CONN(ind_elem, 3)
    
        z1 = GRID%COORD(ind_node_1, 1)
        r1 = GRID%COORD(ind_node_1, 2)
        z2 = GRID%COORD(ind_node_2, 1)
        r2 = GRID%COORD(ind_node_2, 2)
        z3 = GRID%COORD(ind_node_3, 1)
        r3 = GRID%COORD(ind_node_3, 2)
        
        print *, "z, r", z1, r1, z2, r2, z3, r3
        ! Step 3: Get the nodal derivatives
        DX1 = NODE_GRAD_RECOVERY(ind_node_1, 1)
        DY1 = NODE_GRAD_RECOVERY(ind_node_1, 2)
        DXX1 = NODE_GRAD_RECOVERY(ind_node_1, 3)
        DXY1 = NODE_GRAD_RECOVERY(ind_node_1, 4)
        DYY1 = NODE_GRAD_RECOVERY(ind_node_1, 5)
    
        DX2 = NODE_GRAD_RECOVERY(ind_node_2, 1)
        DY2 = NODE_GRAD_RECOVERY(ind_node_2, 2)
        DXX2 = NODE_GRAD_RECOVERY(ind_node_2, 3)
        DXY2 = NODE_GRAD_RECOVERY(ind_node_2, 4)
        DYY2 = NODE_GRAD_RECOVERY(ind_node_2, 5)
    
        DX3 = NODE_GRAD_RECOVERY(ind_node_3, 1)
        DY3 = NODE_GRAD_RECOVERY(ind_node_3, 2)
        DXX3 = NODE_GRAD_RECOVERY(ind_node_3, 3)
        DXY3 = NODE_GRAD_RECOVERY(ind_node_3, 4)
        DYY3 = NODE_GRAD_RECOVERY(ind_node_3, 5)
    
        ! Step 4: Compute the triangle's area
        AREA = CALC_TRIANGLE_AREA(z1, r1, z2, r2, z3, r3)
    
        ! Step 5: Compute the coefficients for shape functions
        b1 = r2 * z3 - r3 * z2
        b2 = r3 * z1 - r1 * z3
        b3 = r1 * z2 - r2 * z1
    
        c1 = z2 - z3
        c2 = z3 - z1
        c3 = z1 - z2
    
        d1 = r3 - r2
        d2 = r1 - r3
        d3 = r2 - r1
    
        ! Step 6: Compute the shape functions at (r, z)
        N_1 = (1 / (2 * AREA)) * (b1 + c1 * r + d1 * z)
        N_2 = (1 / (2 * AREA)) * (b2 + c2 * r + d2 * z)
        N_3 = (1 / (2 * AREA)) * (b3 + c3 * r + d3 * z)
    
        ! Step 7: Interpolate the derivatives at the point (r, z)
        INTERPOLATED_DERIVATIVE(1) = N_1 * DX1 + N_2 * DX2 + N_3 * DX3       ! DX (∂ψ/∂Z)
        INTERPOLATED_DERIVATIVE(2) = N_1 * DY1 + N_2 * DY2 + N_3 * DY3       ! DY (∂ψ/∂R)
        INTERPOLATED_DERIVATIVE(3) = N_1 * DXX1 + N_2 * DXX2 + N_3 * DXX3    ! DXX (∂²ψ/∂Z²)
        INTERPOLATED_DERIVATIVE(4) = N_1 * DXY1 + N_2 * DXY2 + N_3 * DXY3    ! DXY (∂²ψ/∂R∂Z)
        INTERPOLATED_DERIVATIVE(5) = N_1 * DYY1 + N_2 * DYY2 + N_3 * DYY3    ! DYY (∂²ψ/∂R²)
    end function INTERPOLATE_DERIVATIVE_AXISYM
!***************************************************************************
!***************************************************************************
!***************************************************************************
    subroutine CALCULATE_VELOCITY_AT_POINT(z, r, VELOCITY)
        implicit none
        real(8), intent(in) :: z, r                            ! Point coordinates (X = Z, Y = R)
        real(8), intent(out) :: VELOCITY(2)                    ! Output velocity [u_r, u_z]
        real(8) :: interpolated_derivatives(5)                 ! Interpolated [Z_DEV, R_DEV, ZZ_DEV, ZR_DEV, RR_DEV]
    
        ! Interpolate derivatives at the given point (X, Y)
        interpolated_derivatives = INTERPOLATE_DERIVATIVE_AXISYM(z, r)
        ! print*, "interpolated_derivatives", interpolated_derivatives(:)

        ! if (r < 1.0E-6) then
        !     VELOCITY(1) = interpolated_derivatives(5)
        !     VELOCITY(2) = 0.0
        ! else
        !     ! VELOCITY(1) =  (interpolated_derivatives(2) &
        !     !                 + z * interpolated_derivatives(4))/r &
        !     !                 + interpolated_derivatives(5)
        !     ! VELOCITY(2) = - interpolated_derivatives(1) &
        !     !                 - z * interpolated_derivatives(3) &
        !     !                 - r * interpolated_derivatives(4)
        !     VELOCITY(1) = - interpolated_derivatives(2)/r
        !     VELOCITY(2) = interpolated_derivatives(1)/r
        ! end if
        VELOCITY(1) = interpolated_derivatives(1) + r * interpolated_derivatives(4) + z * interpolated_derivatives(3)
        VELOCITY(2) = interpolated_derivatives(2) + r * interpolated_derivatives(5) + z * interpolated_derivatives(4)
    end subroutine CALCULATE_VELOCITY_AT_POINT
    
end module streamline_solution_triangular