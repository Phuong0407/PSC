module grid_generation
    implicit none
    private
    real(8), parameter :: PI = 3.141592653589793
    public :: INIT_GRID, GEN_GRID_COORD, GEN_ELEM_CONN, &
              GEN_NODE_CONN, ADD_NODE_CONNECTION
    type :: GRID_STRUCT
        real(8) :: H_1, H_2, H_3                        ! Horizontal dimension of each region
        real(8) :: L_1, L_2                             ! Vertical dimension of each region
        integer :: eH_1, eH_2, eH_3                     ! Horizontal NUMBER OF ELEMENT for each region
        integer :: indH_1, indH_2, indH_3               ! Horizontal INDEX OF NODE for each region
        integer :: eV, eH, n_elem                       ! Total number of VERTICAL and HORIZONTAL and TOTAL elements
        integer :: n_tri_elem                           ! Total number of TRIANGLULAR element
        integer :: nV, nH, nNode                        ! Total number of VERTICAL and HORIZONTAL and TOTAL nodes
        real(8), allocatable :: COORD(:,:)              ! Coordinate matrices, NODE_IND[X, Y]
        integer, allocatable :: ELEM_CONN(:,:)          ! Connection Data for Assemble: ELEMENT[NODE1_IND, NODE_2_IND, NODE_3_IND, NODE_4_IND]
        ! integer, allocatable :: IND_DIRICHLET_BOUND(:)  ! Connection Data for Boundary Prescription: ELEMENT[NODE1_IND, NODE_2_IND, NODE_3_IND, NODE_4_IND]
        ! integer, allocatable :: EDGE_CONN(:,:)          ! Connection Data for Element Boundary Flux: BOUND[NODE1_IND, NODE_2_IND]
        ! integer, allocatable :: NODE_CONN(:,:)          ! Connection Data for Element Boundary Flux: BOUND[NODE1_IND, NODE_2_IND]
        real(8) :: dx1, dx2, dx3                        ! Finite-Different horizontal dimension of each region
        real(8) :: dy1, dy2                             ! Finite-Different vertical dimension of each region
    end type GRID_STRUCT

    type(grid_struct) :: GRID
contains
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

        GRID%n_elem = GRID%eH * GRID%eV
        ! UPDATE HERE AS WE ADAPT THE QUADRILATERAL MESH TO THE TRIANGULAR MESH,
        ! THEN THE NUMBER OF NODE IS UNCHANGED, BUT THE NUMBER OF ELEMENT INCREASE DOUBLY
        GRID%n_Tri_Elem = 2 * GRID%n_elem
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

end modue grid_generation