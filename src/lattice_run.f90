program use_lm
  use config_m
  use lattice_model
  implicit none

  type(config_t) :: conf
  type(lattice_t) :: l

  integer :: i, j
  real(kind=rk) :: xi

  logical :: do_drive, do_diffusion, do_reaction
  real(kind=rk) :: p_flip, p_move
  integer :: n_bins, n_per_bin
  integer :: width, n_particles
  integer :: left, right, middle
  integer :: n_loops, n_inner

  call conf%init(filename=get_character_argument(1))

  do_diffusion = conf%get_l('do_diffusion')
  do_drive = conf%get_l('do_drive')
  do_reaction = conf%get_l('do_reaction')

  p_flip = conf%get_d('p_flip')
  p_move = conf%get_d('p_move')

  n_bins = conf%get_i('n_bins')
  n_per_bin = conf%get_i('n_per_bin')

  n_particles = conf%get_i('n_particles')
  width = conf%get_i('width')

  n_loops = conf%get_i('n_loops')
  n_inner = conf%get_i('n_inner')

  call l%init(n_bins=n_bins, n_per_bin=n_per_bin, do_drive=do_drive, &
       do_diffusion=do_diffusion, p_move=p_move, p_flip=p_flip)

  middle = n_bins / 2
  right = width / 2
  left = width - 1 - right
  l%n(middle-left:middle+right) = n_particles/width

  if (do_drive) then
     do i = middle-left, middle+right
        do j = 1, l%n(i)
           call random_number(xi)
           if (xi > 0.5_rk) then
              l%v(j, i) = 1
           else
              l%v(j, i) = -1
           end if
        end do
     end do
  end if

  do i = 1, n_loops
     do j = 1, n_inner
        call l%step()
     end do
     write(31, *) l%n
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

end program use_lm
