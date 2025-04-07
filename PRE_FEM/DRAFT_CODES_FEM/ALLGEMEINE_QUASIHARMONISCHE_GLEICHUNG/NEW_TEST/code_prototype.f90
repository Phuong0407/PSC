module ICE_CUBE_AERODYNAMICS
    implicit none
    ! private :: REMOVE_NODE_INSIDE_OBSTACLE, REMAP_NODE_IND, GEN_NODE_CONN, CALC_STIFF_MAT_GLO
    public :: IND_3D_TO_1D
    real(8) :: V_INF, ALPHA, BETA, P_INF
    real(8) :: LENGTH_DOM_1, LENGTH_DOM_2, LENGTH_DOM_3
    real(8) :: OBS_L, OBS_W, OBS_H, CEN_X, CEN_Y, CEN_Z
    integer :: NUM_NODE_X, NUM_NODE_Y, NUM_NODE_Z, NUM_NODE_TOT
    real(8) :: STEP_X, STEP_Y, STEP_Z, STEP_T
    integer, allocatable :: NODE_CONN(:,:)
    integer, allocatable :: NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES(:)
    integer, allocatable :: IND_OBS_BOUND_NODE(:,:), IND_FREE_STREAM_BOUND_NODE(:)
    integer, allocatable :: IND_INSIDE_OBS(:)
    real(8), allocatable :: GRID_COORDS(:,:) ! IND
    real(8), allocatable :: P(:,:,:), U(:,:,:), V(:,:,:), W(:,:,:), RHO(:,:,:)
contains
    pure function IND_3D_TO_1D(i, j, k, N_X, N_Y) result(IND)
        integer, intent(in) :: i, j, k
        integer, intent(in) :: N_X, N_Y
        integer :: IND
        IND = (k - 1) * (N_X * N_Y) + (j - 1) * N_X + (i - 1) + 1
    end function IND_3D_TO_1D

    function IS_ON_CUBE_FACE (x, y, z) result(FACE)
        implicit none
        real(8), intent(in) :: x, y, z
        integer :: FACE
        
        FACE = 0
    
        if (DABS(x - CEN_X) == OBS_L / 2.0 .and. DABS(y - CEN_Y) <= OBS_W / 2.0 .and. DABS(z - CEN_Z) <= OBS_H / 2.0) then
            if (x - CEN_X < 0.0) then
                FACE = 1 ! Left FACE (-X)
            else
                FACE = 2 ! Right FACE (+X)
            end if
        else if (DABS(y - CEN_Y) == OBS_W / 2.0 .and. DABS(x - CEN_X) <= OBS_L / 2.0 .and. DABS(z - CEN_Z) <= OBS_H / 2.0) then
            if (y - CEN_Y < 0.0) then
                FACE = 3 ! Bottom FACE (-Y)
            else
                FACE = 4 ! Top FACE (+Y)
            end if
        else if (DABS(z - CEN_Z) == OBS_H / 2.0 .and. DABS(x - CEN_X) <= OBS_L / 2.0 .and. DABS(y - CEN_Y) <= OBS_W / 2.0) then
            if (z - CEN_Z < 0.0) then
                FACE = 5 ! Back FACE (-Z)
            else
                FACE = 6 ! Front FACE (+Z)
            end if
        end if
    end function IS_ON_CUBE_FACE
    
    function IS_ON_FREE_STREAM_BOUND(x, y, z) result(FACE)
        implicit none
        real(8), intent(in) :: x, y, z
        real(8), parameter :: TOL = 1.0E-8
        integer :: FACE
    
        if (DABS(x) <= TOL) then
            FACE = 1
        else if (DABS(x - LENGTH_DOM_1) <= TOL) then
            FACE = 2
        else if (DABS(y) <= TOL) then
            FACE = 3
        else if (DABS(y - LENGTH_DOM_2) <= TOL) then
            FACE = 4
        else if (DABS(z) <= TOL) then
            FACE = 5
        else if (DABS(z - LENGTH_DOM_3) <= TOL) then
            FACE = 6
        else
            FACE = 0
        end if
    end function IS_ON_FREE_STREAM_BOUND
    
    
    function IS_STRICT_INSIDE_CUBE(x, y, z) result(IS_INSIDE)
        implicit none
        real(8) :: x, y, z
        logical :: IS_INSIDE
        
        IS_INSIDE = & 
            (DABS(x - CEN_X) < OBS_L / 2.0) .AND. &
            (DABS(y - CEN_Y) < OBS_W / 2.0) .AND. &
            (DABS(z - CEN_Z) < OBS_H / 2.0)
    end function IS_STRICT_INSIDE_CUBE
    
    subroutine INIT_SOLUTION(N_X, N_Y, N_Z)
        implicit none
        integer, intent(in) :: N_X, N_Y, N_Z
        NUM_NODE_TOT = N_X * N_Y * N_Z
        allocate(GRID_COORDS(NUM_NODE_TOT, 3))
    end subroutine INIT_SOLUTION
    
    subroutine GEN_GRID_CONFIG()
        use omp_lib
        implicit none
        integer :: IND, i, j, k, FACE_IND, FREE_FACE_IND, TOTAL_INSIDE, TOTAL_OBS_BOUND, TOTAL_FREE_BOUND
        integer, allocatable :: TEMP_ARR(:), TEMP_ARR_FACE(:,:)
        real(8) :: x, y, z

        STEP_X = LENGTH_DOM_1 / (NUM_NODE_X - 1)
        STEP_Y = LENGTH_DOM_2 / (NUM_NODE_Y - 1)
        STEP_Z = LENGTH_DOM_3 / (NUM_NODE_Z - 1)
    
        TOTAL_INSIDE = 0
        TOTAL_OBS_BOUND = 0
        TOTAL_FREE_BOUND = 0
    
        allocate(GRID_COORDS(NUM_NODE_TOT, 3))
        allocate(IND_INSIDE_OBS(NUM_NODE_TOT))
        allocate(IND_OBS_BOUND_NODE(NUM_NODE_TOT, 2))
        allocate(IND_FREE_STREAM_BOUND_NODE(NUM_NODE_TOT))
    
        !$OMP PARALLEL DO PRIVATE(i, j, k, x, y, z, IND, FACE_IND, FREE_FACE_IND) REDUCTION(+:TOTAL_INSIDE, TOTAL_OBS_BOUND, TOTAL_FREE_BOUND)
        do k = 1, NUM_NODE_Z
            do j = 1, NUM_NODE_Y
                do i = 1, NUM_NODE_X
                    IND = IND_3D_TO_1D(i, j, k, NUM_NODE_X, NUM_NODE_Y)
                    x = (i - 1) * STEP_X
                    y = (j - 1) * STEP_Y
                    z = (k - 1) * STEP_Z
    
                    GRID_COORDS(IND, :) = (/x, y, z/)
    
                    if (IS_STRICT_INSIDE_CUBE(x, y, z)) then
                        TOTAL_INSIDE = TOTAL_INSIDE + 1
                        IND_INSIDE_OBS(TOTAL_INSIDE) = IND
                    else
                        FACE_IND = IS_ON_CUBE_FACE(x, y, z)
                        if (FACE_IND /= 0) then
                            TOTAL_OBS_BOUND = TOTAL_OBS_BOUND + 1
                            IND_OBS_BOUND_NODE(TOTAL_OBS_BOUND, 1) = IND
                            IND_OBS_BOUND_NODE(TOTAL_OBS_BOUND, 2) = FACE_IND
                        end if
                        FREE_FACE_IND = IS_ON_FREE_STREAM_BOUND(x, y, z)
                        if (FREE_FACE_IND /= 0) then
                            TOTAL_FREE_BOUND = TOTAL_FREE_BOUND + 1
                            IND_FREE_STREAM_BOUND_NODE(TOTAL_FREE_BOUND) = IND
                        end if
                    end if
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    
        if (TOTAL_INSIDE > 0) then
            allocate(TEMP_ARR(TOTAL_INSIDE))
            TEMP_ARR(:) = IND_INSIDE_OBS(1:TOTAL_INSIDE)
            deallocate(IND_INSIDE_OBS)
            call move_alloc(TEMP_ARR, IND_INSIDE_OBS)
        else
            if (allocated(IND_INSIDE_OBS)) then
                deallocate(IND_INSIDE_OBS)
            end if
        end if
        
        if (TOTAL_OBS_BOUND > 0) then
            allocate(TEMP_ARR_FACE(TOTAL_OBS_BOUND, 2))
            TEMP_ARR_FACE(:,:) = IND_OBS_BOUND_NODE(1:TOTAL_OBS_BOUND, :)
            deallocate(IND_OBS_BOUND_NODE)
            call move_alloc(TEMP_ARR_FACE, IND_OBS_BOUND_NODE)
        else
            if (allocated(IND_OBS_BOUND_NODE)) then
                deallocate(IND_OBS_BOUND_NODE)
            end if
        end if
        if (TOTAL_FREE_BOUND > 0) then
            allocate(TEMP_ARR(TOTAL_FREE_BOUND))
            TEMP_ARR(:) = IND_FREE_STREAM_BOUND_NODE(1:TOTAL_FREE_BOUND)
            deallocate(IND_FREE_STREAM_BOUND_NODE)
            call move_alloc(TEMP_ARR, IND_FREE_STREAM_BOUND_NODE)
        else
            if (allocated(IND_FREE_STREAM_BOUND_NODE)) then
                deallocate(IND_FREE_STREAM_BOUND_NODE)
            end if
        end if
    end subroutine GEN_GRID_CONFIG

    subroutine GEN_NEW_IND()
        implicit none
        integer :: NUM_NODE_EXCLUDE_INSIDE_NODES
        integer :: i, NEW_NODE_IND_COUNTER
    
        NUM_NODE_EXCLUDE_INSIDE_NODES = NUM_NODE_TOT - size(IND_INSIDE_OBS)
        allocate(NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES(NUM_NODE_EXCLUDE_INSIDE_NODES))
        NEW_NODE_IND_COUNTER = 1
    
        do i = 1, NUM_NODE_TOT
            if (.not. IS_IND_INSIDE_OBS(i)) then
                NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES(NEW_NODE_IND_COUNTER) = i
                NEW_NODE_IND_COUNTER = NEW_NODE_IND_COUNTER + 1
            end if
        end do
    
        if (NEW_NODE_IND_COUNTER /= NUM_NODE_EXCLUDE_INSIDE_NODES + 1) then
            deallocate(NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES)
            return
        end if
    end subroutine GEN_NEW_IND

    function IS_IND_INSIDE_OBS(IND) result(IS_INSIDE)
        implicit none
        integer, intent(in) :: IND
        logical :: IS_INSIDE
        IS_INSIDE = any(IND_INSIDE_OBS == IND)
    end function IS_IND_INSIDE_OBS

    subroutine GEN_OBS_BOUND_COND()
        implicit none
        integer :: i, IND_NODE
        real(8) :: x, y, z
        real(8) :: V_INF1, V_INF2, V_INF3

        do i = 1, size(NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES)
            IND_NODE = NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES(i)
            x = GRID_COORDS(IND_NODE, 1)
            y = GRID_COORDS(IND_NODE, 2)
            z = GRID_COORDS(IND_NODE, 3)

            if (IS_ON_CUBE_FACE(x, y, z) /= 0) then
                VEL_VALS(IND_NODE, :) = 0.0D0
            else if (IS_ON_FREE_STREAM_BOUND(x, y, z) /= 0) then
                PRESS_VALS(IND_NODE, :) = P_INF
                VEL_VALS(IND_NODE, 1) = V_INF1
                VEL_VALS(IND_NODE, 2) = V_INF2
                VEL_VALS(IND_NODE, 3) = V_INF3
            end if
        end do
    end subroutine GEN_OBS_BOUND_COND

    function IS_ACTIVE_NODES(IND_NODE) result(IS_ACTIVE)
        implicit none
        integer, intent(in) :: IND_NODE
        integer :: i
        real(8) :: x, y, z
        logical :: IS_ACTIVE

        x = GRID_COORDS(IND_NODE, 1)
        y = GRID_COORDS(IND_NODE, 2)
        z = GRID_COORDS(IND_NODE, 3)

        if (IS_ON_FREE_STREAM_BOUND(x, y, z) /= 0) then
            IS_ACTIVE = .false.
        else if (IS_ON_CUBE_FACE(x, y, z) /= 0) then
            IS_ACTIVE = .false.
        else if (.not. IS_STRICT_INSIDE_CUBE(x, y, z)) then
            IS_ACTIVE = .false.
        else
            IS_ACTIVE = .true.
        end if

    end function IS_ACTIVE_NODES

    subroutine CALC_SOLUTION_NAVIER_STOKES()
        implicit none
        integer :: i, IND_NODE
        real(8) :: x, y, z
        do i = 1, size(NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES)
            IND_NODE = NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES(i)
            x = GRID_COORDS(IND_NODE, 1)
            y = GRID_COORDS(IND_NODE, 2)
            z = GRID_COORDS(IND_NODE, 3)

            if (IS_ACTIVE_NODES(i)) then

            end if
        end do
    End subroutine CALC_SOLUTION_NAVIER_STOKES
end module ICE_CUBE_AERODYNAMICS
