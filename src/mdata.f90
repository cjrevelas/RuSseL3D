 module mdata
    integer,allocatable::ix(:,:)
    real(8), allocatable ::xc(:,:)
    integer::dmnum, dmel
    integer numel, nel, ndm,numnp  
    integer:: fcnum,fcel
    integer, allocatable::fcelement(:,:)
    integer,allocatable::fcentity(:)
    
end module