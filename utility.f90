

module utility
  implicit none

contains
  integer(8) function hash(arr, size) 
    implicit none
    real, intent(in) :: arr(size)
    integer, intent(in) :: size
    integer(8) :: ret
    external hash_double, hash_float
    select case(kind(arr))
    case(4)
       call hash_float(arr, size, ret)
    case(8)
       call hash_double(arr, size, ret)
    case default
       stop "Error: hash(): unknown real type"
    end select
    hash = ret
  end function hash
end module utility
  
