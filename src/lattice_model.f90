module lattice_model
  implicit none

  integer, parameter :: rk = selected_real_kind(12)

  type lattice_t
     integer, allocatable :: n(:)
     integer, allocatable :: n_buffer(:)
     integer, allocatable :: v(:,:)
     integer, allocatable :: v_buffer(:,:)
     integer :: rho_0
     logical :: do_diffusion = .false.
     logical :: do_drive = .false.
     logical :: do_reaction = .false.
     logical :: cst_bc = .false.
     real(kind=rk) :: p_flip, p_drive
     real(kind=rk) :: dt
     real(kind=rk) :: k1, k2
   contains
     procedure :: init
     procedure :: step
     procedure :: move
     procedure :: react
     procedure :: reaction_step
  end type lattice_t

contains

  subroutine init(l, n_bins, n_per_bin, do_diffusion, do_drive, do_reaction, p_flip, p_drive, dt)
    class(lattice_t), intent(out) :: l
    integer, intent(in) :: n_bins
    integer, intent(in), optional :: n_per_bin
    logical, intent(in), optional :: do_diffusion
    logical, intent(in), optional :: do_drive
    logical, intent(in), optional :: do_reaction
    real(kind=rk), intent(in), optional :: p_flip
    real(kind=rk), intent(in), optional :: p_drive
    real(kind=rk), intent(in) :: dt
 
    allocate(l%n(n_bins))
    l%n = 0
    allocate(l%n_buffer(n_bins))
    l%n_buffer = 0

    l%dt = dt

    if (present(do_drive)) then
       if (do_drive) then
          l%do_drive = do_drive
          if (.not. present(n_per_bin)) then
             stop 'n_per_bin absent in init'
          end if
          if (.not. present(p_flip)) then
             stop 'p_flip absent in init'
          end if
          if (.not. present(p_drive)) then
             stop 'p_drive absent in init'
          end if
          allocate(l%v(n_per_bin, n_bins))
          allocate(l%v_buffer(n_per_bin, n_bins))
          l%p_flip = p_flip
          l%p_drive = p_drive
       end if
    end if

    if (present(do_diffusion)) then
       if (do_diffusion) then
          l%do_diffusion = do_diffusion
       end if
    end if

    if (present(do_reaction)) then
       if (do_reaction) then
          l%do_reaction = do_reaction
          ! set the rates
       end if
    end if

  end subroutine init

  !> Perform one step of joint diffusive-driven motion
  subroutine step(l)
    class(lattice_t), intent(inout) :: l

    integer :: i, j
    integer :: n_bins, n_per_bin
    integer :: jump, old_n1
    real(kind=rk) :: xi

    n_bins = size(l%n)
    if (l%do_drive) then
       n_per_bin = size(l%v, dim=1)
    else
       n_per_bin = 0
    end if

    do i = 1, n_bins
       do j = 1, l%n(i)
          if (l%do_drive) then
             call random_number(xi)
             if (xi < l%p_drive) then
                l%v(j,i) = l%v(j,i) * random_flip(proba=l%p_flip)
                jump = l%v(j,i)
             else
                jump = 0
             end if
          else
             jump = 0
          end if
          if (l%do_diffusion) then
             jump = jump + random_step(dt=l%dt)
          end if

          call l%move(i, j, i + jump)

       end do
    end do

    ! move buffer data to actual arrays

    l%n = l%n_buffer
    if (l%do_drive) l%v = l%v_buffer

    l%n_buffer = 0

    if (l%cst_bc) then
       old_n1 = l%n(1)
       l%n(1) = l%rho_0

       ! Generate -1, 1 initial velocities for the new particles
       if (l%do_drive) then
          if (old_n1 < l%rho_0) then
             do i = old_n1+1, l%rho_0
                l%v(i, 1) = random_flip(proba=0.5_rk)
             end do
          end if
       end if
    end if

  end subroutine step

  !> Move a single particle to a new bin (in the buffer bin)
  subroutine move(l, from_bin, from_idx, to_bin)
    class(lattice_t), intent(inout) :: l
    integer, intent(in) :: from_bin
    integer, intent(in) :: from_idx
    integer, intent(in) :: to_bin

    integer :: n_bins, n_per_bin, idx, to_bin_mod

    n_bins = size(l%n)

    ! If the boundary condition is constant on the left, drop
    ! outgoing particles
    if (l%cst_bc .and. &
         ((to_bin <=0 ) .or. (to_bin > n_bins))) return

    to_bin_mod = modulo(to_bin - 1, n_bins) + 1

    if (l%do_drive) then
       ! check if there is room to store the velocity
       n_per_bin = size(l%v_buffer, dim=1)

       if (l%n_buffer(to_bin_mod) <= n_per_bin) then
          idx = l%n_buffer(to_bin_mod) + 1
          l%v_buffer(idx, to_bin_mod) = l%v(from_idx, from_bin)
       else
          stop 'size exceeded in move'
       end if
    end if

    l%n_buffer(to_bin_mod) = l%n_buffer(to_bin_mod) + 1

  end subroutine move

  subroutine reaction_step(l)
    class(lattice_t), intent(inout) :: l

    integer :: i, n_bins, r, n, n_per_bin

    n_bins = size(l%n)

    if (l%do_drive) then
       n_per_bin = size(l%v, dim=1)
    else
       n_per_bin = 0
    end if

    do i = 1, n_bins
       r = l%react(i)

       if (r==0) cycle

       if (r==1) then
          ! create_particle
          n = l%n(i) + 1
          if (l%do_drive) then
             if (n>n_per_bin) then
                write(*,*) 'number of particles exceeded in reaction_step'
                stop
             end if
             l%v(n, i) = random_flip(0.5_rk)
          end if
          l%n(i) = n
       else if (r==2) then
          ! kill 2
          l%n(i) = l%n(i) - 2
       end if

    end do

  end subroutine reaction_step

  function react(l, i) result(r)
    class(lattice_t), intent(in) :: l
    integer, intent(in) :: i
    integer :: r

    integer :: n
    real(kind=rk) :: xi, proba_1, proba_2

    n = l%n(i)

    if (n == 0) then
       r = 0
       return
    end if

    proba_1 = n*l%k1*l%dt

    if (n>1) then
       proba_2 = n*(n-1)*l%k2*l%dt/2
    else
       proba_2 = 0
    end if

    call random_number(xi)

    r = 0
    if (xi < proba_1 + proba_2) then
       call random_number(xi)
       xi = xi*(proba_1+proba_2)
       if (xi < proba_1) then
          r = 1
       else
          r = 2
       end if
    end if

  end function react

  !> Return 1 with probability 1-proba and -1 with probability proba
  function random_flip(proba) result(r)
    real(kind=rk), intent(in) :: proba
    integer :: r

    real(kind=rk) :: xi

    call random_number(xi)

    if (xi < proba) then
       r = -1
    else
       r = 1
    end if

  end function random_flip

  !> Return +/- 1 with probability proba and 0 with probability 1-proba
  function random_step(dt) result(r)
    real(kind=rk), intent(in) :: dt
    integer :: r

    real(kind=rk) :: xi

    call random_number(xi)

    if (xi < dt) then
       r = 1
    else if (xi < 2*dt) then
       r = -1
    else
       r = 0
    end if

  end function random_step

end module lattice_model
