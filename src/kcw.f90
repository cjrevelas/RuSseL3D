module kcw
    
    double precision,allocatable :: c(:,:)
    double precision, allocatable :: k(:,:)
    double precision, allocatable :: w(:,:)
    double precision, allocatable :: g(:,:)
    double precision, allocatable :: rh(:,:)
    type mum_matrix
        sequence 
        DOUBLE PRECISION, DIMENSION(:), POINTER :: value
        INTEGER, DIMENSION(:), POINTER :: row, col
        !integer::col
        !integer::row
        !double precision:: m
    end type 
    integer non_zero
    
    
    
    type( mum_matrix)::c_m
    end 
