program test_csr_solver
    use csr_sparse_matrix
    implicit none

    real, dimension(10, 10) :: A
    real, dimension(10) :: RHS_VECT, SOLUTION
    integer :: i

    A = reshape([&
    10.75925616,  0.15334631,  0.0       ,  0.0       ,  0.0       , &
     0.93835208,  0.0       ,  0.0       ,  0.77274546,  0.58800104, &
     0.15334631, 10.0       ,  0.0       ,  0.0       ,  0.0       , &
     0.0       ,  0.73700786,  0.0       ,  0.27559238,  0.69211309, &
     0.0       ,  0.0       , 10.0       ,  0.45147133,  0.52801772, &
     0.0       ,  0.0       ,  0.07429782,  0.44942516,  0.0       , &
     0.0       ,  0.0       ,  0.45147133, 10.0       ,  0.0       , &
     0.0       ,  0.0       ,  0.85894124,  0.0       ,  0.0       , &
     0.0       ,  0.0       ,  0.52801772,  0.0       , 10.0       , &
     0.0       ,  0.0       ,  0.0       ,  0.0       ,  0.0       , &
     0.93835208,  0.0       ,  0.0       ,  0.0       ,  0.0       , &
    10.0       ,  0.49383478,  0.0       ,  0.0       ,  0.0       , &
     0.0       ,  0.73700786,  0.0       ,  0.0       ,  0.0       , &
     0.49383478, 10.26045892,  1.34211025,  0.0       ,  0.0       , &
     0.0       ,  0.0       ,  0.07429782,  0.85894124,  0.0       , &
     0.0       ,  1.34211025, 10.0       ,  0.0       ,  0.39783897, &
     0.77274546,  0.27559238,  0.44942516,  0.0       ,  0.0       , &
     0.0       ,  0.0       ,  0.0       , 11.11748848,  0.0       , &
     0.58800104,  0.69211309,  0.0       ,  0.0       ,  0.0       , &
     0.0       ,  0.0       ,  0.39783897,  0.0       , 11.02878898], &
    shape=[10, 10])


    RHS_VECT = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    SOLUTION = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    call INIT_SPARSE_SOLVER(A)

    print *, "Solving the system A * X = RHS_VECT..."
    SOLUTION = SOLVE_SPARSE_SYSTEM(RHS_VECT)

    print *, "Solution vector:"
    do i = 1, size(SOLUTION)
        print *, "X(", i, ") = ", SOLUTION(i)
    end do

    print *, "Test completed."
end program test_csr_solver