module ssa_model_m
  implicit none

  private

  public :: rk
  public :: ssa_model_t

  integer, parameter :: rk = selected_real_kind(12)

  type ssa_model_t
     integer :: n_sites
     real(kind=rk) :: time
     real(kind=rk) :: k_move
     real(kind=rk) :: k_flip
     real(kind=rk) :: k_1
     real(kind=rk) :: k_2
     integer, allocatable :: n_left(:)
     integer, allocatable :: n_right(:)
     real(kind=rk), allocatable :: k(:)
   contains
     procedure :: init => ssa_model_init
     procedure :: step => ssa_step
  end type ssa_model_t

contains

  subroutine ssa_model_init(this, n_sites)
    class(ssa_model_t), intent(out) :: this
    integer, intent(in) :: n_sites

    this%n_sites = n_sites
    this%time = 0
    allocate(this%n_left(n_sites))
    this%n_left = 0
    allocate(this%n_right(n_sites))
    this%n_right = 0
    allocate(this%k(n_sites))
    this%k = 0

  end subroutine ssa_model_init

  subroutine ssa_step(this)
    class(ssa_model_t), intent(inout) :: this

    integer :: i, n_l, n_r, n
    real(kind=rk) :: k_move, k_flip, k_1, k_2
    real(kind=rk) :: k, k_total, xi, cumulative_k, alpha

    k_move = this%k_move
    k_flip = this%k_flip
    k_1 = this%k_1
    k_2 = this%k_2

    k_total = 0
    do i = 1, this%n_sites
       n_l = this%n_left(i)
       n_r = this%n_right(i)
       n = n_l + n_r

       k = n * (k_flip + k_move) + n*k_1 + n*(n-1)*k_2

       k_total = k_total + k

       this%k(i) = k
    end do

    call random_number(xi)
    this%time = this%time - log(xi)/k_total

    call random_number(xi)
    xi = xi*k_total

    cumulative_k = this%k(1)
    do i = 2, this%n_sites
       cumulative_k = cumulative_k + this%k(i)
       ! react cell i
       if (xi < cumulative_k) exit
    end do

    n_l = this%n_left(i)
    n_r = this%n_right(i)
    n = n_l + n_r

    k = n * (k_flip + k_move + k_1) + n*(n-1)*k_2

    call random_number(xi)
    xi = xi*k

    if (xi < n_l*k_flip) then
       ! flip a left-going particle
       this%n_left(i) = this%n_left(i) - 1
       this%n_right(i) = this%n_right(i) + 1
    else if (xi < n*k_flip) then
       ! flip a right-going particle
       this%n_right(i) = this%n_right(i) - 1
       this%n_left(i) = this%n_left(i) + 1
    else if (xi < n*k_flip + n_l*k_move) then
       ! move a left-going particle
       this%n_left(i) = this%n_left(i) - 1
       if (i-1 >0) this%n_left(i-1) = this%n_left(i-1) + 1
    else if (xi < n * k_flip + n*k_move) then
       ! move a right-going particle
       this%n_right(i) = this%n_right(i) - 1
       if (i+1 <= this%n_sites) this%n_right(i+1) = this%n_right(i+1) + 1
    else if (xi < n * (k_flip + k_move + k_1)) then
       ! add a particle
       call random_number(alpha)
       if (alpha < 0.5_rk) then
          this%n_right(i) = this%n_right(i) + 1
       else
          this%n_left(i) = this%n_left(i) + 1
       end if
    else if (xi < n * (k_flip + k_move + k_1) + n*(n-1)*k_2) then
       ! remove a particle
       call random_number(alpha)
       if (alpha < 0.5_rk) then
          this%n_right(i) = this%n_right(i) - 1
       else
          this%n_left(i) = this%n_left(i) - 1
       end if
    else
       write(*,*) 'error in ssa_step'
    end if

  end subroutine ssa_step

end module ssa_model_m
