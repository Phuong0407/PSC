module ICE_CUBE_AERODYNAMICS
    implicit none
    real(8) :: U_INF, ALPHA, BETA, P_INF, RHO_INF, ARTI_COMPRE
    real(8) :: LENGTH_DOM_1, LENGTH_DOM_2, LENGTH_DOM_3
    real(8) :: OBS_L, OBS_W, OBS_H, X_CEN, Y_CEN, Z_CEN
    integer :: NUM_NODE_X, NUM_NODE_Y, NUM_NODE_Z, NUM_NODE_TOT
    ! integer, allocatable :: IND_OBS_BOUND_NODES(:,:), IND_FREE_STREAM_BOUND_NODES(:), IND_INSIDE_OBS(:)
    real(8) :: STEP_X, STEP_Y, STEP_Z, STEP_T
    real(8), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
    real(8), allocatable :: P(:,:,:), U(:,:,:), V(:,:,:), W(:,:,:), RHO(:,:,:)
contains
    function IS_ON_CUBE_FACE(X_IND, Y_IND, Z_IND) result(FACE)
        implicit none
        integer, intent(in) :: X_IND, Y_IND, Z_IND
        real(8), parameter :: TOL = 1.0E-8
        real(8) :: x_loc, y_loc, z_loc
        real(8) :: x_len, y_len, z_len
        integer :: FACE

        x_loc = x(X_IND, Y_IND, Z_IND)
        y_loc = y(X_IND, Y_IND, Z_IND)
        z_loc = z(X_IND, Y_IND, Z_IND)

        x_len = DABS(x_loc - X_CEN)
        y_len = DABS(y_loc - Y_CEN)
        z_len = DABS(z_loc - Z_CEN)

        if (DABS(x_len - OBS_L / 2.0) <= TOL .and. y_len <= OBS_W / 2.0 .and. z_len <= OBS_H / 2.0) then
            if (x_loc - X_CEN < 0.0) then
                FACE = 1
            else
                FACE = 2
            end if
        else if (DABS(y_len - OBS_W / 2.0) <= TOL .and. x_len <= OBS_L / 2.0 .and. z_len <= OBS_H / 2.0) then
            if (y_loc - Y_CEN < 0.0) then
                FACE = 3
            else
                FACE = 4
            end if
        else if (DABS(z_len - OBS_H / 2.0) <= TOL .and. x_len <= OBS_L / 2.0 .and. y_len <= OBS_W / 2.0) then
            if (z_loc - Z_CEN < 0.0) then
                FACE = 5
            else
                FACE = 6
            end if
        else
            FACE = 0
        end if
    end function IS_ON_CUBE_FACE
    
    function IS_ON_FREE_STREAM_BOUND(X_IND, Y_IND, Z_IND) result(FACE)
        implicit none
        integer, intent(in) :: X_IND, Y_IND, Z_IND
        real(8), parameter :: TOL = 1.0E-8
        real(8) :: x_loc, y_loc, z_loc
        integer :: FACE

        x_loc = x(X_IND, Y_IND, Z_IND)
        y_loc = y(X_IND, Y_IND, Z_IND)
        z_loc = z(X_IND, Y_IND, Z_IND)

        if (DABS(x_loc) <= TOL) then
            FACE = 1
        else if (DABS(x_loc - LENGTH_DOM_1) <= TOL) then
            FACE = 2
        else if (DABS(y_loc) <= TOL) then
            FACE = 3
        else if (DABS(y_loc - LENGTH_DOM_2) <= TOL) then
            FACE = 4
        else if (DABS(z_loc) <= TOL) then
            FACE = 5
        else if (DABS(z_loc - LENGTH_DOM_3) <= TOL) then
            FACE = 6
        else
            FACE = 0
        end if
    end function IS_ON_FREE_STREAM_BOUND
    
    
    function IS_STRICT_INSIDE_CUBE(X_IND, Y_IND, Z_IND) result(IS_INSIDE)
        implicit none
        integer, intent(in) :: X_IND, Y_IND, Z_IND
        real(8) :: x_loc, y_loc, z_loc
        real(8) :: x_len, y_len, z_len
        logical :: IS_INSIDE

        x_loc = x(X_IND, Y_IND, Z_IND)
        y_loc = y(X_IND, Y_IND, Z_IND)
        z_loc = z(X_IND, Y_IND, Z_IND)
        x_len = DABS(x_loc - X_CEN)
        y_len = DABS(y_loc - Y_CEN)
        z_len = DABS(z_loc - Z_CEN)

        IS_INSIDE = & 
            (DABS(x_len - OBS_L / 2.0) < 0) .AND. &
            (DABS(y_len - OBS_W / 2.0) < 0) .AND. &
            (DABS(z_len - OBS_H / 2.0) < 0)
    end function IS_STRICT_INSIDE_CUBE
    
    subroutine INIT_SOLUTION(N_X, N_Y, N_Z)
        implicit none
        integer, intent(in) :: N_X, N_Y, N_Z

        NUM_NODE_TOT = N_X * N_Y * N_Z

        STEP_X = LENGTH_DOM_1 / (NUM_NODE_X - 1)
        STEP_Y = LENGTH_DOM_2 / (NUM_NODE_Y - 1)
        STEP_Z = LENGTH_DOM_3 / (NUM_NODE_Z - 1)
    end subroutine INIT_SOLUTION
    
    subroutine GEN_GRID_COORDS()
        implicit none
        integer :: i, j, k

        allocate(x(NUM_NODE_X, NUM_NODE_Y, NUM_NODE_Z))
        allocate(y(NUM_NODE_X, NUM_NODE_Y, NUM_NODE_Z))
        allocate(z(NUM_NODE_X, NUM_NODE_Y, NUM_NODE_Z))
    
        do k = 1, NUM_NODE_Z
            do j = 1, NUM_NODE_Y
                do i = 1, NUM_NODE_X
                    x(i,j,k) = (i - 1) * STEP_X
                    y(i,j,k) = (j - 1) * STEP_Y
                    z(i,j,k) = (k - 1) * STEP_Z
                end do
            end do
        end do
    end subroutine GEN_GRID_COORDS

    ! function IS_IND_INSIDE_OBS(IND) result(IS_INSIDE)
    !     implicit none
    !     integer, intent(in) :: IND
    !     logical :: IS_INSIDE
    !     IS_INSIDE = any(IND_INSIDE_OBS == IND)
    ! end function IS_IND_INSIDE_OBS

    ! subroutine GEN_OBS_BOUND_COND()
    !     implicit none
    !     integer :: i, IND_NODE
    !     real(8) :: x, y, z
    !     real(8) :: V_INF1, V_INF2, V_INF3

    !     do i = 1, size(NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES)
    !         IND_NODE = NEW_IND_TO_OLD_IND_EXCLUDE_INSIDE_NODES(i)
    !         x = GRID_COORDS(IND_NODE, 1)
    !         y = GRID_COORDS(IND_NODE, 2)
    !         z = GRID_COORDS(IND_NODE, 3)

    !         if (IS_ON_CUBE_FACE(i, j, k) /= 0) then
    !             VEL_VALS(IND_NODE, :) = 0.0D0
    !         else if (IS_ON_FREE_STREAM_BOUND(i, j, k) /= 0) then
    !             PRESS_VALS(IND_NODE, :) = P_INF
    !             VEL_VALS(IND_NODE, 1) = V_INF1
    !             VEL_VALS(IND_NODE, 2) = V_INF2
    !             VEL_VALS(IND_NODE, 3) = V_INF3
    !         end if
    !     end do
    ! end subroutine GEN_OBS_BOUND_COND

    ! function IS_ACTIVE_NODES(IND_NODE) result(IS_ACTIVE)
    !     implicit none
    !     integer, intent(in) :: IND_NODE
    !     integer :: i
    !     real(8) :: x, y, z
    !     logical :: IS_ACTIVE

    !     x = GRID_COORDS(IND_NODE, 1)
    !     y = GRID_COORDS(IND_NODE, 2)
    !     z = GRID_COORDS(IND_NODE, 3)

    !     if (IS_ON_FREE_STREAM_BOUND(i, j, k) /= 0) then
    !         IS_ACTIVE = .false.
    !     else if (IS_ON_CUBE_FACE(i, j, k) /= 0) then
    !         IS_ACTIVE = .false.
    !     else if (.not. IS_STRICT_INSIDE_CUBE(i, j, k)) then
    !         IS_ACTIVE = .false.
    !     else
    !         IS_ACTIVE = .true.
    !     end if

    ! end function IS_ACTIVE_NODES

    subroutine GEN_OBS_BOUND_COND()
        implicit none
        integer :: i, j, k
        do i = 1, NUM_NODE_X
            do j = 1, NUM_NODE_Y
                do k = 1, NUM_NODE_Z
                    if(IS_ON_CUBE_FACE(i, j, k) /= 0) then
                        U(i, j, k) = 0.0
                        V(i, j, k) = 0.0
                        W(i, j, k) = 0.0
                    end if
                end do
            end do
        end do
    end subroutine GEN_OBS_BOUND_COND

    subroutine GEN_FREE_STREAM_COND()
        implicit none
        real(8) :: U_INF1, U_INF2, U_INF3
    
        U_INF1 = U_INF * cos(ALPHA) * cos(BETA)
        U_INF2 = U_INF * sin(ALPHA) * cos(BETA)
        U_INF3 = U_INF * sin(BETA)
    
        U(1, :, :) = U_INF1
        V(1, :, :) = U_INF2
        W(1, :, :) = U_INF3
        P(1, :, :) = P_INF
        RHO(1, :, :) = P_INF / ARTI_COMPRE
    
        U(NUM_NODE_X, :, :) = U_INF1
        V(NUM_NODE_X, :, :) = U_INF2
        W(NUM_NODE_X, :, :) = U_INF3
        P(NUM_NODE_X, :, :) = P_INF
        RHO(NUM_NODE_X, :, :) = P_INF / ARTI_COMPRE
    
        U(:, 1, :) = U_INF1
        V(:, 1, :) = U_INF2
        W(:, 1, :) = U_INF3
        P(:, 1, :) = P_INF
        RHO(:, 1, :) = P_INF / ARTI_COMPRE
    
        U(:, NUM_NODE_Y, :) = U_INF1
        V(:, NUM_NODE_Y, :) = U_INF2
        W(:, NUM_NODE_Y, :) = U_INF3
        P(:, NUM_NODE_Y, :) = P_INF
        RHO(:, NUM_NODE_Y, :) = P_INF / ARTI_COMPRE
    
        U(:, :, 1) = U_INF1
        V(:, :, 1) = U_INF2
        W(:, :, 1) = U_INF3
        P(:, :, 1) = P_INF
        RHO(:, :, 1) = P_INF / ARTI_COMPRE
    
        U(:, :, NUM_NODE_Z) = U_INF1
        V(:, :, NUM_NODE_Z) = U_INF2
        W(:, :, NUM_NODE_Z) = U_INF3
        P(:, :, NUM_NODE_Z) = P_INF
        RHO(:, :, NUM_NODE_Z) = P_INF / ARTI_COMPRE    
    end subroutine GEN_FREE_STREAM_COND
    

    subroutine CALC_SOLUTION_NAVIER_STOKES()
        implicit none
        
    End subroutine CALC_SOLUTION_NAVIER_STOKES
end module ICE_CUBE_AERODYNAMICS
