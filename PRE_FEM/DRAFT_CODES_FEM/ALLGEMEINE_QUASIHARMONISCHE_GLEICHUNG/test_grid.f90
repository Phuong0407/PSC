program test_grid
    use streamline_solution
    ! use csr_sparse_matrix
    implicit none
    integer :: i
    integer :: j
    integer :: k
    integer :: file_unit
    ! logical :: is_positive_definite
    ! real, allocatable :: L(:,:)
    ! real :: sum
    ! logical :: is_diagonally_dominant
    ! integer :: i, j
    ! real :: row_sum
    ! real, allocatable :: AFD(:)
    ! real, allocatable :: rAFD(:)
    ! integer :: i
    call INIT_SOLUTION(3.0, 3.0, 3.0, 4.0, 2.0, 2, 2, 2, 2, 1.0)
    ! print *, GRID%nElem

    ! do i = 1, GRID%nNode
    !     write(*, '(f10.4, 1x)', advance="no") &
    !         "[", GRID%COORD(i, 1), ", ", &
    !             GRID%COORD(i, 2), ", ", &
    !         "] "
    !     print *  ! End the current row and move to the next line
    ! end do

    ! do i = 1, GRID%nElem
    !         write(*, '(A, I6, A, I6, A, I6, A, I6, A)', advance="no") &
    !             "[", GRID%CONN(i, 1), ", ", &
    !                  GRID%CONN(i, 2), ", ", &
    !                  GRID%CONN(i, 3), ", ", &
    !                  GRID%CONN(i, 4), "] "
    !     print *  ! End of row
    ! end do
    !===========================
    ! PRINT TO FILE
    ! file_unit = 20  ! Logical unit number
    ! open(unit=file_unit, file='stiffness_matrix.txt', status='replace', action='write')

    ! ! Write the matrix to the file
    ! do i = 1, GRID%nNode
    !     do j = 1, GRID%nNode
    !         write(file_unit, '(f10.4, 1x)', advance="no") STIFF_MAT_GLO(i, j)
    !     end do
    !     write(file_unit, *)  ! End the current row and move to the next line
    ! end do
    !===========================

    !===========================
    ! PRINT TO SCREEN
    do i = 1, GRID%nNode
        do j = 1, GRID%nNode
            write(*, '(f10.4, 1x)', advance="no") STIFF_MAT_GLO(i, j)
        end do
        write(*, *)
    end do

    do i = 1, GRID%nNode
        write(*, '(f10.4, 1x)', advance="no") FORCE_MAT_GLO(i)
        write(*, *)
    end do
    !===========================


    ! Close the file
    ! close(file_unit)
    
    ! do i = 1, GRID%nNode
    !     do j = 1, GRID%nNode
    !         if (STIFF_MAT_GLO(i, j) /= STIFF_MAT_GLO(j, i)) then
    !             print *, "STIFFNESS MATRIX IS NOT SYMMETRIC"
    !             exit
    !         end if
    !     end do
    ! end do
    ! print *, "STIFFNESS MATRIX IS SYMMETRIC"
    
    ! do i = 1, GRID%nNode
    !     print *, DIRICHLET_BOUND(i, 1)
    ! end do

    ! ! Check if the matrix is positive definite using Cholesky decomposition
    ! allocate(L(GRID%nNode, GRID%nNode))
    ! L = 0.0d0
    ! is_positive_definite = .true.
    ! do i = 1, GRID%nNode
    !     do j = 1, i
    !         sum = STIFF_MAT_GLO(i, j)
    !             do k = 1, j-1
    !                 sum = sum - L(i, k) * L(j, k)
    !             end do
    !             if (i == j) then
    !                 if (sum <= 0.0d0) then
    !                     is_positive_definite = .false.
    !                     exit
    !                 end if
    !                 L(i, j) = sqrt(sum)
    !             else
    !             L(i, j) = sum / L(j, j)
    !         end if
    !     end do
    !         if (.not. is_positive_definite) exit
    ! end do

    ! ! Output results
    ! if (is_positive_definite) then
    !     print *, "The matrix is positive definite."
    ! else
    !     print *, "The matrix is not positive definite."
    ! end if

    ! do i = 1, GRID%nNode
    !     row_sum = 0.0d0
    !     do j = 1, GRID%nNode
    !        if (j /= i) then
    !           row_sum = row_sum + abs(STIFF_MAT_GLO(i, j))
    !        end if
    !     end do
    !     if (abs(STIFF_MAT_GLO(i, i)) <= row_sum) then
    !        is_diagonally_dominant = .false.
    !        exit
    !     end if
    ! end do
   
    !  ! Output result
    ! if (is_diagonally_dominant) then
    !     print *, "The matrix is diagonally dominant."
    ! else
    !     print *, "The matrix is not diagonally dominant."
    ! end if

    ! allocate(AFD(size(STIFF_MAT_GLO, 1)))
    ! allocate(rAFD(size(STIFF_MAT_GLO, 1)))
    ! AFD = 1.0
    ! ! call INIT_SPARSE_SOLVER(STIFF_MAT_GLO)
    ! ! call MULT_SPARSE_MAT()
    ! ! call PRINT_MAT()
    ! rAFD = matmul(STIFF_MAT_GLO, AFD)
    ! do i = 1, size(STIFF_MAT_GLO, 1)
    !     print *, rAFD(i)
    ! end do
end program test_grid