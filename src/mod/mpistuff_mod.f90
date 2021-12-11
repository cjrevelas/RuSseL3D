!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module mpistuff_mod
    integer(kind=4) :: my_id
    integer(kind=4) :: n_proc
    integer(kind=4) :: ierr
    integer(kind=4) :: world_group_id

    logical :: root
    logical :: flag_continue
end module mpistuff_mod
