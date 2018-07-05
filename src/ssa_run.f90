program ssa_run
  use config_m
  use ssa_model_m
  implicit none

  type(config_t) :: conf
  type(ssa_model_t) :: l

  real(kind=rk) :: xi
  integer :: funit

  logical :: do_drive, do_diffusion, do_reaction
  real(kind=rk) :: k_move, k_flip
  integer :: i, j, n_loops, n_inner
  character(len=:), allocatable :: output_filename
  character(len=:), allocatable :: seed_argument

  call conf%init(filename=get_character_argument(1))
  output_filename = get_character_argument(2)


  call l%init(conf%get_i('n_sites'))
  l%k_move = conf%get_d('k_move')
  l%k_flip = conf%get_d('k_flip')
  l%k_1 = conf%get_d('k_1')
  l%k_2 = conf%get_d('k_2') / conf%get_i('u_0')

  n_loops = conf%get_i('n_loops')
  n_inner = conf%get_i('n_inner')

  l%n_left(l%n_sites/2) = conf%get_i('u_0') / 2
  l%n_right(l%n_sites/2) = conf%get_i('u_0') / 2

  do i = 1, n_loops
     do j = 1, n_inner
        call l%step()
     end do
     write(31,*) l%time
     write(32,*) l%n_left
     write(33,*) l%n_right
  end do


contains

  !> Return the position-th command-line argument as an integer
  function get_character_argument(position) result(arg)
    integer, intent(in) :: position
    character(len=:), allocatable :: arg
    character(len=128) :: r

    if (command_argument_count() < position) then
       stop 'missing argument for parameter file'
    end if

    call get_command_argument(position, r)

    arg = trim(r)

  end function get_character_argument

end program ssa_run
