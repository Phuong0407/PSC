module csr_sparse_matrix
    implicit none
    private
    real(8) :: ZERO_COMPARED_TOLERANCE, MATRIX_NORM    ! AXILIARY VARIABLE to check for NULL-INDICES in sparse matrix
    integer :: SIZE_MAT, NON_ZERO_ELEM              ! Size of matrix and number of non-zero elements
    real(8), allocatable :: VAL(:)                     ! Non-zero values
    integer, allocatable :: COL_IND(:)              ! Column indices of non-zero elements
    integer, allocatable :: ROW_PTR(:)              ! Row pointer array
    real(8) :: CONV_TOL = 1.0E-6                       ! Tolerance for convergence
    real(8) :: MAX_ITER = 100000                       ! Maximum number of iterations

    public :: INIT_SPARSE_SOLVER, MULT_SPARSE_MAT, SOLVE_SPARSE_SYSTEM

contains
    subroutine INIT_SPARSE_SOLVER(SPARSE_MAT)
        implicit none
        real(8), intent(in), dimension(:,:) :: SPARSE_MAT
        integer :: i, j, index

        NON_ZERO_ELEM = 0
        SIZE_MAT = size(SPARSE_MAT, 1)

        MATRIX_NORM = 0.0
        do i = 1, SIZE_MAT
            MATRIX_NORM = max(MATRIX_NORM, sum(abs(SPARSE_MAT(i, :))))
        end do
        ZERO_COMPARED_TOLERANCE = EPSILON(1.0) * MATRIX_NORM

        do i = 1, SIZE_MAT
            do j = 1, SIZE_MAT
                if (abs(SPARSE_MAT(i, j)) > ZERO_COMPARED_TOLERANCE) then
                    NON_ZERO_ELEM = NON_ZERO_ELEM + 1
                end if
            end do
        end do

        allocate(VAL(NON_ZERO_ELEM))
        allocate(COL_IND(NON_ZERO_ELEM))
        allocate(ROW_PTR(SIZE_MAT + 1))

        ROW_PTR = 0
        index = 1
        ROW_PTR(1) = 1
        do i = 1, SIZE_MAT
            do j = 1, SIZE_MAT
                if (abs(SPARSE_MAT(i, j)) > ZERO_COMPARED_TOLERANCE) then
                    VAL(index) = SPARSE_MAT(i, j)
                    COL_IND(index) = j
                    index = index + 1
                end if
            end do
            ROW_PTR(i + 1) = index
        end do
    end subroutine INIT_SPARSE_SOLVER

    function MULT_SPARSE_MAT(X) result(Y)
        implicit none
        real(8), intent(in) :: X(:)
        real(8), allocatable :: Y(:)
        integer :: i, j
    
        if (size(X) /= SIZE_MAT) then
            stop
        end if

        allocate(Y(SIZE_MAT))
        Y = 0.0
        do i = 1, SIZE_MAT
            do j = ROW_PTR(i), ROW_PTR(i + 1) - 1
                Y(i) = Y(i) + VAL(j) * X(COL_IND(j))
            end do
        end do
    end function MULT_SPARSE_MAT
    
    function SOLVE_SPARSE_SYSTEM(RHS_VECT) result(X)
        implicit none
        real(8), intent(in) :: RHS_VECT(:)
        real(8), allocatable :: r(:), p(:), Ap(:), X(:)
        real(8) :: alpha, beta, rnorm, bnorm, tol
        integer :: ITER
    
        allocate(r(SIZE_MAT), p(SIZE_MAT), Ap(SIZE_MAT))
        tol = CONV_TOL
    
        allocate(X(SIZE_MAT))
        X = 0.0
        bnorm = sqrt(dot_product(RHS_VECT, RHS_VECT))
        r = RHS_VECT - MULT_SPARSE_MAT(X)
        p = r
        rnorm = dot_product(r, r)
        
        ITER = 0
        do while (sqrt(rnorm) / bnorm > tol .and. ITER < MAX_ITER)
            ITER = ITER + 1
            Ap = MULT_SPARSE_MAT(p)
            alpha = rnorm / dot_product(p, Ap)
            if (alpha == 0.0) exit
            
            X = X + alpha * p
            r = r - alpha * Ap
            
            if (sqrt(dot_product(r, r)) / bnorm < tol) exit
            beta = dot_product(r, r) / rnorm
            if (beta == 0.0) exit
        
            p = r + beta * p
            rnorm = dot_product(r, r)
        end do
    
        if (ITER == MAX_ITER) then
            print *, "Failed to converge in ", MAX_ITER, " iterations."
        end if
        deallocate(r, p, Ap)
    end function SOLVE_SPARSE_SYSTEM
end module csr_sparse_matrix