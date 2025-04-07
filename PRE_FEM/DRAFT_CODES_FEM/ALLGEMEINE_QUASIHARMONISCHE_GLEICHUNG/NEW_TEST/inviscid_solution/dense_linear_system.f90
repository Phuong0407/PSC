!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
module DENSE_MATRIX_SOLVER
contains
    !  ***************************************************************
    !  * Given an N x N matrix A, this routine replaces it by the LU *
    !  * decomposition of a rowwise permutation of itself. A and N   *
    !  * are input. IND is an output vector which records the row   *
    !  * permutation effected by the partial pivoting; D is output   *
    !  * as -1 or 1, depending on whether the number of row inter-   *
    !  * changes was even or odd, respectively. This routine is used *
    !  * in combination with LUBKSB to solve linear equations or to  *
    !  * invert a matrix. Return code is 1, if matrix is singular.   *
    !  ***************************************************************
    subroutine DECOMPOSE_LU(A, IND, PARITY, IS_SINGULAR)
        implicit none
        real(8), parameter :: TINY = 1.5D-16

        real(8), intent(inout), dimension(:,:) :: A
        integer, intent(out) :: PARITY, IS_SINGULAR
        integer, intent(out), dimension(:) :: IND
        integer :: MAT_SIZE
        real(8)  :: SUM, ROW_MAX_VAL, PIVOT, PIVOT_COMPARATOR, ROW_SWAP, PIVOT_DIVISOR
        real(8), allocatable :: IMPLICIT_SCALING_FACTOR(:)
        INTEGER :: i, j, k, i_MAX

        MAT_SIZE = size(A,1)

        allocate(IMPLICIT_SCALING_FACTOR(MAT_SIZE))

        PARITY = 1
        IS_SINGULAR = 0
    
        do i = 1, MAT_SIZE
            ROW_MAX_VAL = 0.d0
            do j = 1, MAT_SIZE
                if (DABS(A(i,j)) .GT. ROW_MAX_VAL) ROW_MAX_VAL = DABS(A(i,j))
            end do
            if (ROW_MAX_VAL .LT. TINY) then
                IS_SINGULAR = 1
                RETURN
            end if
            IMPLICIT_SCALING_FACTOR(i) = 1.d0 / ROW_MAX_VAL
        end do

        do j =1, MAT_SIZE
            do i = 1, j - 1
                SUM = A(i,j)
                do K=1,I-1
                    SUM = SUM - A(I,K)*A(K,J) 
                end do
                A(i,j) = SUM
            end do

            PIVOT = 0.d0
            do i = j, MAT_SIZE
                SUM = A(i,j)
                do k=1, j - 1
                    SUM = SUM - A(I,K)*A(K,J) 
                end do
                A(i,j) = SUM
                PIVOT_COMPARATOR = IMPLICIT_SCALING_FACTOR(i)*DABS(SUM)
                if (PIVOT_COMPARATOR .GE. PIVOT) then
                    i_MAX = i
                    PIVOT = PIVOT_COMPARATOR
                end if
            end do
            
            if (j .NE. i_MAX) then
                do k = 1, MAT_SIZE
                    ROW_SWAP = A(i_MAX, k)
                    A(i_MAX,K) = A(j, k)
                    A(j, k) = ROW_SWAP
                end do
                PARITY = - PARITY
                IMPLICIT_SCALING_FACTOR(i_MAX) = IMPLICIT_SCALING_FACTOR(j)
            end if

            IND(j) = i_MAX
            if (DABS(A(j, j)) < TINY) A(j, j) = TINY

            if (j .NE. MAT_SIZE) then
                PIVOT_DIVISOR = 1.d0 / A(j, j)
                do i = j + 1, MAT_SIZE
                    A(i,j) = A(i,j) * PIVOT_DIVISOR
                end do
            end if 
        end do

        RETURN
     end subroutine DECOMPOSE_LU
    
    
    !  ******************************************************************
    !  * Solves the set of N linear equations A . X = B.  Here A is     *
    !  * input, not as the matrix A but rather as its LU decomposition, *
    !  * determined by the routine LUDCMP. IND is input as the permuta-*
    !  * tion vector returned by LUDCMP. B is input as the right-hand   *
    !  * side vector B, and returns with the solution vector X. A, N and*
    !  * IND are not modified by this routine and can be used for suc- *
    !  * cessive calls with different right-hand sides. This routine is *
    !  * also efficient for plain matrix inversion.                     *
    !  ******************************************************************
    subroutine SUBSTITUTE_BACKWARD_LU(A, IND, B)
        implicit none
        real(8), intent(in), dimension(:,:) :: A
        integer, intent(in), dimension(:) :: IND
        real(8), intent(inout), dimension(:) :: B
        real(8) :: SUM
        integer :: MAT_SIZE, FIRST_NON_NULL_IND, ROW_PERM_IND, i, j
        
        MAT_SIZE = size(A,1)
        
        FIRST_NON_NULL_IND = 0
        
        do i = 1, MAT_SIZE
            ROW_PERM_IND = IND(i)
            SUM = B(ROW_PERM_IND)
            B(ROW_PERM_IND) = B(i)
            if (FIRST_NON_NULL_IND .NE. 0) then
                do j = FIRST_NON_NULL_IND, i - 1
                    SUM = SUM - A(i,j)*B(j)
                end do
            else if (SUM.NE.0.d0) then
                FIRST_NON_NULL_IND = i
            end if
            B(i) = SUM
        end do
        
        do i = MAT_SIZE , 1, -1
            SUM = B(i)
            if (i < MAT_SIZE) then
                do j = i + 1, MAT_SIZE
                    SUM = SUM - A(i,j)*B(j)
                end do
            end if
            B(i) = SUM / A(I,I)
        end do
    end subroutine SUBSTITUTE_BACKWARD_LU

    function SOLVE_DENSE_MATRIX_SYSTEM(A, B) result(SOLUTION)
        implicit none
        real(8), allocatable :: SOLUTION(:)
        real(8), dimension(:, :), intent(inout) :: A(:,:)
        real(8), dimension(:), intent(in) :: B(:)
        integer,pointer ::  IND(:)
        integer :: IS_SINGULAR, MAT_SIZE, PARITY
        MAT_SIZE = size(A,1)
        
        allocate(IND(MAT_SIZE))
        allocate(SOLUTION(MAT_SIZE))
        SOLUTION = B

        call DECOMPOSE_LU(A, IND, PARITY, IS_SINGULAR)
        if (IS_SINGULAR .EQ. 0) then
            call SUBSTITUTE_BACKWARD_LU(A, IND, SOLUTION)
        end if
    end function SOLVE_DENSE_MATRIX_SYSTEM
end module DENSE_MATRIX_SOLVER