module grid_generation
    implicit none
    private
    public :: grid_struct

    type :: grid_struct
        real :: H_1, H_2, H_3       ! Widths of the regions
        real :: L_1, L_2            ! Heights of the regions
        integer :: N                ! Vertical number of element divisions
        integer :: M_1, M_2, M_3    ! Horizontal element divisions for each region
        real, allocatable :: x(:,:), y(:,:)  ! Coordinate matrices
    contains
        procedure :: init_grid, gen_grid
    end type grid_struct

contains

    subroutine init_grid(this, H_1, H_2, H_3, L_1, L_2, N, M_1, M_2, M_3)
        class(grid_struct), intent(inout) :: this
        real, intent(in) :: H_1, H_2, H_3, L_1, L_2
        integer, intent(in) :: N, M_1, M_2, M_3

        this%N = N
        this%M_1 = M_1
        this%M_2 = M_2
        this%M_3 = M_3

        this%H_1 = H_1
        this%H_2 = H_2
        this%H_3 = H_3
        this%L_1 = L_1
        this%L_2 = L_2

        allocate(this%x(N+1, M_1 + M_2 + M_3 + 1))
        allocate(this%y(N+1, M_1 + M_2 + M_3 + 1))
    end subroutine init_grid

    subroutine gen_grid(this)
        class(grid_struct), intent(inout) :: this
        integer :: i, j
        real :: x_2

        do i = 0, this%N
            ! Region 1: First Rectangular Region
            do j = 0, this%M_1
                this%x(i+1, j+1) = this%H_1 / this%M_1 * j
                this%y(i+1, j+1) = this%L_1 / this%N * i
            end do

            ! Region 2: Oblique Region
            do j = this%M_1 + 1, this%M_1 + this%M_2
                x_2 = this%H_1 + (this%H_2 / this%M_2) * (j - this%M_1)
                this%x(i+1, j+1) = x_2
                this%y(i+1, j+1) = ((this%L_2 - this%L_1) / (this%H_2 * this%N)) * i * x_2 + &
                                   (this%L_1 / this%N) * i - ((this%L_2 - this%L_1) / (this%H_2 * this%N)) * i * this%H_1
            end do

            ! Region 3: Second Rectangular Region
            do j = this%M_1 + this%M_2 + 1, this%M_1 + this%M_2 + this%M_3
                this%x(i+1, j+1) = this%H_1 + this%H_2 + (this%H_3 / this%M_3) * (j - this%M_1 - this%M_2)
                this%y(i+1, j+1) = this%L_2 / this%N * i
            end do
        end do
    end subroutine gen_grid

end module grid_generation