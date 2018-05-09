program use_lm
  use config_m
  use lattice_model
  implicit none

  type(config_t) :: conf
  type(lattice_t) :: l

  integer :: i, j
  real(kind=rk) :: xi
  integer :: funit

  logical :: do_drive, do_diffusion, do_reaction
  real(kind=rk) :: p_flip, p_drive
  real(kind=rk) :: k1, k2
  integer :: n_bins, n_per_bin
  integer :: width, n_particles
  integer :: left, right, middle
  integer :: n_loops, n_inner
  character(len=:), allocatable :: output_filename

  call conf%init(filename=get_character_argument(1))
  output_filename = get_character_argument(2)

  do_diffusion = conf%get_l('do_diffusion')
  do_drive = conf%get_l('do_drive')
  do_reaction = conf%get_l('do_reaction')

  p_flip = conf%get_d('p_flip')
  p_drive = conf%get_d('p_drive')

  n_bins = conf%get_i('n_bins')
  n_per_bin = conf%get_i('n_per_bin')

  n_particles = conf%get_i('n_particles')
  width = conf%get_i('width')

  n_loops = conf%get_i('n_loops')
  n_inner = conf%get_i('n_inner')

  call l%init(n_bins=n_bins, n_per_bin=n_per_bin, do_drive=do_drive, &
       do_diffusion=do_diffusion, p_flip=p_flip, &
       p_drive=p_drive, &
       dt=conf%get_d('dt'), &
       do_reaction=do_reaction)

  k1 = conf%get_d('k1')
  k2 = conf%get_d('k2')
  l%rho_0 = conf%get_d('rho_0')
  l%cst_bc = conf%get_l('cst_bc')

  l%k1 = k1
  l%k2 = k2 / l%rho_0

  if (l%cst_bc) then
     left = 1
     right = width
     l%n(1:width) = n_particles/width
  else
     middle = n_bins / 2
     right = middle + width / 2
     left = middle - (width - 1 - width/2)
  end if

  l%n(left:right) = n_particles/width

  if (do_drive) then
     do i = left, right
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

  open(newunit=funit, file=output_filename)

  do i = 1, n_loops
     do j = 1, n_inner
        call l%step()
        if (l%do_reaction) call l%reaction_step()
     end do
     write(funit, '(*(i5))') l%n
  end do

  close(funit)

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
