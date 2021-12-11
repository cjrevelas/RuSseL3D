!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

#ifdef __DO_NOT_PREPROCESS_DOC__
! Hash table implementation imitating to GCC STL (with singly linked list).
! DO NOT COMPILE THIS TEMPLATE FILE DIRECTLY.
! Use a wrapper module and include this file instead, e.g. fhash_modules.f90.
! Remove is not implemented since not needed currently.
!
! #define                         | meaning
! --------------------------------+-----------------------------------------------------
! SHORTNAME <Name>                | (optional) The name of the type this FHASH table is
!                                 | for. If set, it overrides all settings that have
!                                 | have possibly been made for FHASH_MODULE_NAME,
!                                 | FHASH_TYPE_NAME and FHASH_TYPE_ITERATOR_NAME.
!                                 |
! FHASH_MODULE_NAME <Name>        | The name of the module that encapsulates the FHASH
!                                 | types and functionality
! FHASH_TYPE_NAME <Name>          | The name of the actual FHASH type
! FHASH_TYPE_ITERATOR_NAME <Name> | The name of the FHASH type that can iterate through
!                                 | the whole FHASH
!                                 |
! KEY_USE <use stmt>              | (optional) A use statement that is required to use
!                                 | a specific type as a key for the FHASH
! KEY_TYPE <typename>             | The type of the keys. May require KEY_USE to be
!                                 | accessible.
!                                 |
! VALUE_USE <use stmt>            | (optional) A use statement that is required to use
!                                 | a specific type as a value for the FHASH
! VALUE_TYPE <typename>           | The type of the values. May require VALUE_USE to be
!                                 | accessible.
!                                 |
! VALUE_VALUE                     | Flag indicating that the values in FHASH are value
!                                 | values. This is the default. (see VALUE_POINTER)
! VALUE_POINTER                   | Flag indicating that the values in FHASH are value
!                                 | pointers.
! VALUE_ASSIGNMENT                | (internal) The assignment operator, do not set it
!                                 | anywhere, it is configured based on VALUE_VALUE or
!                                 | VALUE_POINTER
#endif

#ifdef SHORTNAME
#undef FHASH_MODULE_NAME
#undef FHASH_TYPE_NAME
#undef FHASH_TYPE_ITERATOR_NAME

#ifdef __GFORTRAN__
#define PASTE(a) a
#define CONCAT(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define CONCAT(a,b) PASTE(a,b)
#endif

#define FHASH_MODULE_NAME CONCAT(fhash_module__,SHORTNAME)
#define FHASH_TYPE_NAME CONCAT(fhash_type__,SHORTNAME)
#define FHASH_TYPE_ITERATOR_NAME CONCAT(fhash_type_iterator__,SHORTNAME)
#endif

#undef VALUE_ASSIGNMENT
#ifndef VALUE_VALUE
#ifndef VALUE_POINTER
#define VALUE_VALUE
#endif
#endif

#ifdef VALUE_POINTER
#define VALUE_ASSIGNMENT =>
#else
#define VALUE_ASSIGNMENT =
#endif

module FHASH_MODULE_NAME !The module starts here (after all the preproc directives)
#undef FHASH_MODULE_NAME

#ifdef KEY_USE
  KEY_USE ! USE statement regarding the key type
#undef KEY_USE
#endif

#ifdef VALUE_USE
  VALUE_USE ! USE statement regarding the value type
#undef VALUE_USE
#endif

  implicit none

  private

  public :: FHASH_TYPE_NAME
  public :: FHASH_TYPE_ITERATOR_NAME

    ! Define (key, value) pair of type (KEY_TYPE, VALUE_TYPE)
    type kv_type
        KEY_TYPE   :: key
        VALUE_TYPE :: value
    end type kv_type

    ! Define node_type as the type of bucket elements
    type node_type
        type(kv_type), allocatable :: kv
        type(node_type), pointer   :: next => NULL()

      contains
        procedure :: node_set
        procedure :: node_get
        procedure :: node_remove
        procedure :: node_clear
        procedure :: node_depth
    end type node_type

    type FHASH_TYPE_NAME
      private
        integer :: n_buckets = 0
        integer :: n_keys    = 0
        type(node_type), allocatable :: buckets(:)

      contains
        procedure, public :: bucket_count
        procedure, public :: n_collisions
        procedure, public :: reserve
        procedure, public :: key_count
        procedure, public :: set
        procedure, public :: get
        procedure, public :: remove
        procedure, public :: clear
    end type FHASH_TYPE_NAME

    type FHASH_TYPE_ITERATOR_NAME
      private
        integer :: bucket_id
        type(node_type),       pointer :: node_ptr  => NULL()
        type(FHASH_TYPE_NAME), pointer :: fhash_ptr => NULL()

      contains
        procedure, public :: begin ! Set the iterator to the beginning of a hash table (i.e., the first bucket element)
        procedure, public :: next  ! Get the (key, value) of the next bucket element and advance the iterator
    end type FHASH_TYPE_ITERATOR_NAME

  contains ! Module procedures are defined below

    function bucket_count(this)
        class(FHASH_TYPE_NAME), intent(inout) :: this
        integer :: bucket_count

        bucket_count = this%n_buckets
    end function bucket_count

    function n_collisions(this)
        class(FHASH_TYPE_NAME), intent(inout) :: this
        integer :: n_collisions
        integer :: ii

        n_collisions = 0
        do ii = 1, this%n_buckets
            n_collisions = n_collisions + node_depth(this%buckets(ii)) - 1
        enddo
    end function n_collisions

    ! Return the length of the linked list starting from the current node
    recursive function node_depth(this) result(depth)
        class(node_type), intent(inout) :: this
        integer :: depth

        if (.not. ASSOCIATED(this%next)) then
            depth = 1
        else
            depth = 1 + node_depth(this%next)
        endif
    end function node_depth

    subroutine reserve(this, n_buckets)
        class(FHASH_TYPE_NAME), intent(inout) :: this
        integer, intent(in)    :: n_buckets
        integer, dimension(29) :: sizes
        integer                :: ii

        if (this%key_count() > 0) STOP "Cannot reserve when fhash is not empty."

        sizes = (/5, 11, 23, 47, 97, 199, 409, 823, 1741, 3469, 6949, 14033,     &
                & 28411, 57557, 116731, 236897, 480881, 976369,1982627, 4026031, &
                & 8175383, 16601593, 33712729, 68460391, 139022417, 282312799,   &
                & 573292817, 1164186217, 2147483647/)

        do ii = 1, SIZE(sizes)
            if (sizes(ii) >= n_buckets) then
                this%n_buckets = sizes(ii)
                allocate(this%buckets(this%n_buckets))
                return
            endif
        enddo
    end subroutine reserve

    function key_count(this)
        class(FHASH_TYPE_NAME), intent(inout) :: this
        integer :: key_count

        key_count = this%n_keys
    end function key_count

    subroutine set(this, key, value)
        class(FHASH_TYPE_NAME), intent(inout) :: this
        KEY_TYPE, intent(in)   :: key
        VALUE_TYPE, intent(in) :: value
        integer                :: bucket_id
        logical                :: is_new

        bucket_id = MODULO(hash_value(key), this%n_buckets) + 1 ! This is the actual hash function

        call this%buckets(bucket_id)%node_set(key, value, is_new)

        if (is_new) this%n_keys = this%n_keys + 1
    end subroutine set

    ! If kv is not allocated, allocate and set to the (key, value) passed in
    ! If key is present and the same as the key passed in, overwrite the value
    ! If key is present and different than the key passed in, defer to the next node (allocate if not allocated)
    recursive subroutine node_set(this, key, value, is_new)
        class(node_type), intent(inout) :: this
        KEY_TYPE, intent(in)            :: key
        VALUE_TYPE, intent(in)          :: value
        logical, optional, intent(out)  :: is_new

        if (.not. ALLOCATED(this%kv)) then
            allocate(this%kv)
            this%kv%key = key
            this%kv%value VALUE_ASSIGNMENT value
            if (PRESENT(is_new)) is_new = .true.
        else if (this%kv%key == key) then
            this%kv%value VALUE_ASSIGNMENT value
            if (PRESENT(is_new)) is_new = .false.
        else ! A collision has occured: different (key, value) pair into the same bucket id
            if (.not. ASSOCIATED(this%next)) allocate(this%next)
            call this%next%node_set(key, value, is_new)
        endif
    end subroutine node_set

    subroutine get(this, key, value, success)
        class(FHASH_TYPE_NAME), intent(inout) :: this
        KEY_TYPE, intent(in)                  :: key
        VALUE_TYPE, intent(out)               :: value
        logical, optional, intent(out)        :: success
        integer :: bucket_id

        bucket_id = MODULO(hash_value(key), this%n_buckets) + 1

        call this%buckets(bucket_id)%node_get(key, value, success)
    end subroutine get

    ! If kv is not allocated, fail and return 0.
    ! If key is present and the same as the key passed in, return the value in kv.
    ! If next pointer is associated, delegate to it.
    ! Otherwise, fail and return 0.
    recursive subroutine node_get(this, key, value, success)
        class(node_type), intent(inout) :: this
        KEY_TYPE, intent(in)            :: key
        VALUE_TYPE, intent(out)         :: value
        logical, optional, intent(out)  :: success

        if (.not. ALLOCATED(this%kv)) then
            ! Not found. (Initial node in the bucket not set)
            if (PRESENT(success)) success = .false.
        else if (this%kv%key == key) then
            value VALUE_ASSIGNMENT this%kv%value
            if (PRESENT(success)) success = .true.
        else if (ASSOCIATED(this%next)) then
            call this%next%node_get(key, value, success)
        else
            if (PRESENT(success)) success = .false.
        endif
    end subroutine node_get

    subroutine remove(this, key, success)
        class(FHASH_TYPE_NAME), intent(inout) :: this
        KEY_TYPE, intent(in)           :: key
        logical, optional, intent(out) :: success
        integer                        :: bucket_id
        type(node_type)                :: first
        logical                        :: locSuccess

        bucket_id = MODULO(hash_value(key), this%n_buckets) + 1

        first = this%buckets(bucket_id)

        if (ALLOCATED(first%kv)) then
            if (first%kv%key == key) then
                if (ASSOCIATED(first%next)) then
                    this%buckets(bucket_id)%kv%key = this%buckets(bucket_id)%next%kv%key
                    this%buckets(bucket_id)%kv%value VALUE_ASSIGNMENT this%buckets(bucket_id)%next%kv%value
                    deallocate(first%next%kv)
                    this%buckets(bucket_id)%next => this%buckets(bucket_id)%next%next
                else
                    deallocate(this%buckets(bucket_id)%kv)
                endif
                locSuccess = .true.
            else
                call node_remove(first%next, key, locSuccess, first)
            endif
        else
            locSuccess = .false.
        endif

        if (locSuccess) this%n_keys = this%n_keys - 1

        if (PRESENT(success)) success = locSuccess
    end subroutine remove

    ! If kv is not allocated, fail and return
    ! If key is present and node is first in bucket, set first node in bucket to
    !   the next node of first. Return success
    ! If key is present and the node is another member of the linked list, link the
    !   previous node's next node to this node's next node, deallocate this node,
    !   return success
    ! Otherwise, fail and return 0
    recursive subroutine node_remove(this, key, success, last)
        class(node_type), intent(inout) :: this, last
        KEY_TYPE, intent(in) :: key
        logical, intent(out) :: success

        if (.not. ALLOCATED(this%kv)) then
            ! Not found. (Initial node in the bucket not set)
            success = .false.
        else if (this%kv%key == key) then
            last%next => this%next
            nullify(this%next)
            deallocate(this%kv)
            success = .true.
        else if (ASSOCIATED(this%next)) then
            call this%next%node_remove(key, success, this)
        else
            success = .false.
        endif
    end subroutine node_remove

    subroutine clear(this)
        class(FHASH_TYPE_NAME), intent(inout) :: this
        integer :: ii

        if (.not. ALLOCATED(this%buckets)) return

        do ii = 1, SIZE(this%buckets)
            if (ASSOCIATED(this%buckets(ii)%next)) then
                call this%buckets(ii)%next%node_clear()
                deallocate(this%buckets(ii)%next)
            endif
        enddo
        deallocate(this%buckets)
        this%n_keys    = 0
        this%n_buckets = 0
    end subroutine clear

    ! Deallocate kv if allocated.
    ! Call the clear method of the next node if the next pointer associated.
    ! Deallocate and nullify the next pointer.
    recursive subroutine node_clear(this)
        class(node_type), intent(inout) :: this

        if (ASSOCIATED(this%next)) then
            call this%next%node_clear()
            deallocate(this%next)
            nullify(this%next)
        endif
    end subroutine node_clear

    subroutine begin(this, fhash_target)
        class(FHASH_TYPE_ITERATOR_NAME), intent(inout) :: this
        type(FHASH_TYPE_NAME), target, intent(in) :: fhash_target

        this%bucket_id = 1
        this%node_ptr  => fhash_target%buckets(1)
        this%fhash_ptr => fhash_target
    end subroutine begin

    subroutine next(this, key, value, status)
        class(FHASH_TYPE_ITERATOR_NAME), intent(inout) :: this
        KEY_TYPE, intent(out)          :: key
        VALUE_TYPE, intent(out)        :: value
        integer, optional, intent(out) :: status

        do while (.not. ASSOCIATED(this%node_ptr) .or. .not. ALLOCATED(this%node_ptr%kv))
            if (this%bucket_id < this%fhash_ptr%n_buckets) then
                this%bucket_id = this%bucket_id + 1
                this%node_ptr => this%fhash_ptr%buckets(this%bucket_id)
            else
                if (PRESENT(status)) status = -1
#ifdef VALUE_TYPE_INIT
                value VALUE_ASSIGNMENT INT(VALUE_TYPE_INIT)
#endif
                return
            endif
         enddo

         key = this%node_ptr%kv%key
         value VALUE_ASSIGNMENT this%node_ptr%kv%value
         if (PRESENT(status)) status = 0
         this%node_ptr => this%node_ptr%next
    end subroutine next
end module

#undef KEY_TYPE
#undef VALUE_TYPE
#undef VALUE_TYPE_INIT
#undef VALUE_ASSIGNMENT
#undef FHASH_TYPE_NAME
#undef FHASH_TYPE_ITERATOR_NAME
#undef SHORTNAME
#undef CONCAT
#undef PASTE
