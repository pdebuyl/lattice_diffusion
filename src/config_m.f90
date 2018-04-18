module config_m
  implicit none
  private

  public :: config_t

  integer, parameter :: dp = selected_real_kind(12)

  !> Single entry in the dictionary
  type entry_t
     character(len=:), allocatable :: key
     character(len=:), allocatable :: value
  end type entry_t

  !> A bucket contains several entries
  type bucket_t
     type(entry_t), allocatable :: entries(:)
     integer :: current_size = 0
     integer :: current_idx = 0
   contains
     procedure :: find
  end type bucket_t

  !> The dictionary contains dict_size buckets (defined at run time)
  type config_t
     type(bucket_t), allocatable :: buckets(:)
     integer :: dict_size = 0
   contains
     procedure :: djb2
     procedure :: set
     procedure :: init
     procedure :: show
     procedure :: get_s
     procedure :: get_i
     procedure :: get_d
     procedure :: get_l
  end type config_t

  integer, parameter :: BUCKET_EMPTY = -2
  integer, parameter :: BUCKET_ENTRY_NOT_FOUND = -4

contains

  !> djb2 hash function
  !!
  !! \param this the config_t object
  !! \param s a string
  !!
  !! \return the hash value between 0 and dict_size-1
  function djb2(this, s) result(r)
    class(config_t), intent(in) :: this
    character(len=*), intent(in) :: s
    integer :: r

    integer :: i, l

    l = len(s)

    r = 5381

    do i = 1, l
       r = r*33 + ichar(s(i:i))
    end do

    r = modulo(r, this%dict_size)

  end function djb2

  !> Add or replace an entry in the dictionary
  !!
  !! \param this the config_t object
  !! \param k the key
  !! \param v the value
  subroutine set(this, k, v)
    class(config_t), intent(inout) :: this
    character(len=*), intent(in) :: k
    character(len=*), intent(in) :: v

    type(bucket_t) :: tmp_bucket

    integer :: h, i, b_idx

    h = this%djb2(k) + 1

    b_idx = this%buckets(h)%find(k)

    if (b_idx == BUCKET_EMPTY) then
       ! allocate bucket for 1 entry
       ! also, means we can take the first entry
       allocate(this%buckets(h)%entries(1))
       this%buckets(h)%current_size = 1
       this%buckets(h)%current_idx = 1
       b_idx = 1
       this%buckets(h)%entries(1)%key = trim(k)
       this%buckets(h)%entries(1)%value = trim(v)
       ! the values are registered, exit
       return
    end if

    if (b_idx == BUCKET_ENTRY_NOT_FOUND) then
       ! copy and grow bucket entries
       
       allocate(tmp_bucket%entries(this%buckets(h)%current_size + 1))
       tmp_bucket%current_size = this%buckets(h)%current_size + 1
       tmp_bucket%current_idx = this%buckets(h)%current_idx + 1

       do i = 1, this%buckets(h)%current_size
          tmp_bucket%entries(i)%key = this%buckets(h)%entries(i)%key
          tmp_bucket%entries(i)%value = this%buckets(h)%entries(i)%value
       end do

       deallocate(this%buckets(h)%entries)
       allocate(this%buckets(h)%entries, source=tmp_bucket%entries)
       deallocate(tmp_bucket%entries)

       this%buckets(h)%current_size = tmp_bucket%current_size
       this%buckets(h)%current_idx = tmp_bucket%current_idx
       b_idx = this%buckets(h)%current_idx
    end if

    if (b_idx > 0) then
       this%buckets(h)%entries(b_idx)%key = k
       this%buckets(h)%entries(b_idx)%value = v
    end if

  end subroutine set

  !> Initialize a dictionary object
  !!
  !! \param this the config_t object
  !! \param dict_size the size of the hash table
  subroutine init(this, filename, dict_size)
    class(config_t), intent(out) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: dict_size

    integer :: funit, iostat, eq_location
    character(len=256) :: line
    character(len=:), allocatable :: k, v

    if (present(dict_size)) then
       allocate(this%buckets(dict_size))
       this%dict_size = dict_size
    else
       allocate(this%buckets(1024))
       this%dict_size = 1024
    end if

    open(file=filename, newunit=funit)
    read_loop: do
       read(funit, '(a)', iostat=iostat) line
       if (iostat < 0) exit read_loop
       if (iostat > 0) cycle read_loop
       eq_location = index(line, '=')
       if (eq_location == 0) cycle read_loop

       k = trim(adjustl(line(1:eq_location-1)))
       v = trim(adjustl(line(eq_location+1:)))

       if (len(k) ==0) then
          write(*,*) 'empty key while reading data'
          stop
       end if
       if (len(v) ==0) then
          write(*,*) 'empty value while reading data'
          stop
       end if

       call this%set(k, v)

    end do read_loop

  end subroutine init

  !> Display the content of a dictionary
  !!
  !! \param this the config_t object
  subroutine show(this)
    class(config_t), intent(in) :: this

    integer :: i, j, s
    integer :: n

    n = 0
    do i = 1, this%dict_size
       s = this%buckets(i)%current_idx
       if (s > 0) then
             write(*,*) 'bucket   : ', i, ' size ', s
          do j = 1, s
             write(*,*) 'key      : ', this%buckets(i)%entries(j)%key
             write(*,*) 'value    : ', this%buckets(i)%entries(j)%value
          end do
       end if
    end do

  end subroutine show

  !> Find the "in-bucket" index for a given key
  !!
  !! Negative return values correspond to module-defined return codes.
  !!
  !! \param this the bucket_t object
  !! \param k the key
  !!
  !! \return the index (1-based) of the key in the bucket or a return code
  function find(this, k) result(r)
    class(bucket_t), intent(in) :: this
    character(len=*), intent(in) :: k
    integer :: r

    integer :: i

    if (this%current_size == 0) then
       r = BUCKET_EMPTY
       return
    end if

    r = BUCKET_ENTRY_NOT_FOUND
    do i = 1, this%current_size
       if (this%entries(i)%key == trim(k)) then
          r = i
          exit
       end if
    end do

  end function find

  !> Fetch an entry in the dictionary.
  !!
  !! \param this the config_t object
  !! \param k the key
  !!
  !! \return the value if found, an empty string else
  function get_s(this, k) result(r)
    class(config_t), intent(in) :: this
    character(len=*), intent(in) :: k

    character(len=:), allocatable :: r

    integer :: h, b_idx

    h = this%djb2(k) + 1

    b_idx = this%buckets(h)%find(k)

    if ( (b_idx == BUCKET_EMPTY) .or. &
         (b_idx == BUCKET_ENTRY_NOT_FOUND) ) then
       write(*,*) 'key ', k, ' not found'
       stop
    end if

    if (b_idx>0) then
       r = this%buckets(h)%entries(b_idx)%value
    end if

  end function get_s

  !> Fetch an entry in the dictionary.
  !!
  !! \param this the config_t object
  !! \param k the key
  !!
  !! \return the value if found, an empty string else
  function get_i(this, k) result(r)
    class(config_t), intent(in) :: this
    character(len=*), intent(in) :: k
    integer :: r

    character(len=:), allocatable :: v

    v = this%get_s(k)

    read(v, *) r

  end function get_i

  !> Fetch an entry in the dictionary.
  !!
  !! \param this the config_t object
  !! \param k the key
  !!
  !! \return the value if found, an empty string else
  function get_d(this, k) result(r)
    class(config_t), intent(in) :: this
    character(len=*), intent(in) :: k
    real(kind=dp) :: r

    character(len=:), allocatable :: v

    v = this%get_s(k)

    read(v, *) r

  end function get_d

  !> Fetch an entry in the dictionary.
  !!
  !! \param this the config_t object
  !! \param k the key
  !!
  !! \return the value if found, an empty string else
  function get_l(this, k) result(r)
    class(config_t), intent(in) :: this
    character(len=*), intent(in) :: k
    logical :: r

    character(len=:), allocatable :: v

    v = this%get_s(k)

    read(v, *) r

  end function get_l

end module config_m
