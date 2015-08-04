! engtraj format ver 0.1.0 (HDF5):
!   attribute
!     /version
!     /sltspec        ! solute species index of this engtraj
!     /sltfirsttag    ! tag of the 1st solute in this engtraj
!     /numslt
!     /nummol
!     /numpair
!     /numframe
!   dataset
!     /engergy        ! when described in column-major: energy(1:numframe, 1:numpair)
!                     ! where numpair is the number of unordered pairs among
!                     ! molecules of tagpt(sltfirsttag:nummol)

module engtraj
  use HDF5
  use mpiproc, only: mpi_comm_activeprocs, MPI_INFO_NULL, myrank
  implicit none
  character(len=*), parameter :: ENGTRAJ_VER = "0.1.0"
  character(len=*), parameter :: ENGTRAJ_FILENAME = "engtraj.h5"
  character(len=*), parameter :: ENGDSET_NAME = "energy"
  integer(hid_t) :: file_id       ! File identifier 
  integer(hid_t) :: dset_id       ! Dataset identifier 
  integer(hid_t) :: filespace     ! Dataspace identifier in file 
  integer :: sltfirsttag, num_relevant_mol

contains
  subroutine engtraj_init(sltspec, slttag1, numslt, nummol, maxcnf, skpcnf, engdiv)
    use H5LT
    implicit none
    integer(hid_t) :: plist_id      ! Property list identifier 
    integer, intent(in) :: sltspec, slttag1, numslt, nummol, maxcnf, skpcnf, engdiv
    integer :: err, numframe, num_record_per_frame, frame_remainder
    logical, save :: is_first_init = .true.

    sltfirsttag = slttag1
    numframe = maxcnf / skpcnf
    frame_remainder = mod(numframe, engdiv)
    if (is_first_init) then
      if (myrank == 0) then
        if (frame_remainder /= 0) then
          write(*, *) "Warning: ", frame_remainder, " frames at the end are not processed due to numframe not divisible by engdiv"
          write(*, *) "Set engdiv to 1 or a number dividing ", numframe, " if you want all the frames to be processed"
        end if
      end if
    end if
    numframe = numframe - frame_remainder

    num_relevant_mol = nummol - sltfirsttag + 1
    num_record_per_frame = numslt * (numslt - 1) / 2 + numslt * (num_relevant_mol - numslt)

    ! Initialize FORTRAN predefined datatypes
    call H5open_f(err)

    ! Setup file access property list with parallel I/O access.
    call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
    call H5Pset_fapl_mpio_f(plist_id, mpi_comm_activeprocs, MPI_INFO_NULL, err)

    if (is_first_init) then
      ! Create the file collectively.
      call H5Fcreate_f(ENGTRAJ_FILENAME, H5F_ACC_EXCL_F, file_id, err, access_prp = plist_id)
    else
      ! Open existing file collectively.
      call H5Fopen_f(ENGTRAJ_FILENAME, H5F_ACC_RDWR_F, file_id, err, access_prp = plist_id)
    end if

    ! Close property list
    call h5pclose_f(plist_id, err)

    if (is_first_init) then
      call H5LTset_attribute_string_f(file_id, "/", "version", ENGTRAJ_VER, err)
      call H5LTset_attribute_int_f(file_id, "/", "sltspec", [sltspec], 1_hsize_t, err)
      call H5LTset_attribute_int_f(file_id, "/", "sltfirsttag", [sltfirsttag], 1_hsize_t, err)
      call H5LTset_attribute_int_f(file_id, "/", "numslt", [numslt], 1_hsize_t, err)
      call H5LTset_attribute_int_f(file_id, "/", "nummol", [nummol], 1_hsize_t, err)
      call H5LTset_attribute_int_f(file_id, "/", "numpair", [num_record_per_frame], 1_hsize_t, err)
      call H5LTset_attribute_int_f(file_id, "/", "numframe", [numframe], 1_hsize_t, err)
    end if

    ! create filespace
    call H5Screate_simple_f(2, [int(numframe, kind=hsize_t), int(num_record_per_frame, kind=hsize_t)], filespace, err)

    if (is_first_init) then
      ! create dataset
      call H5Dcreate_f(file_id, ENGDSET_NAME, H5T_NATIVE_DOUBLE, filespace, dset_id, err)
    else
      ! open existing dataset
      call H5Dopen_f(file_id, ENGDSET_NAME, dset_id, err)
    end if

    if (is_first_init) is_first_init = .false.
  end subroutine engtraj_init

  subroutine engtraj_write(stnum, slttag, energy)
    implicit none
    integer, intent(in) :: stnum, slttag
    real, intent(in) :: energy(:)
    integer(hsize_t) :: offset(2), count(2), stride(2), block(2), size_energy(1)
    integer(hid_t) :: memspace      ! Dataspace identifier in memory
    integer(hid_t) :: plist_id      ! Property list identifier 
    integer :: err, tag, begin_index

    size_energy = [size(energy)]
    tag = slttag - sltfirsttag + 1
    begin_index = (tag - 1) * num_relevant_mol + (tag + 1) - (tag + 1 ) * tag / 2
    !          c
    !    | 1  2  3  4
    !  --+------------
    !  1 |    1  2  3
    !    |
    !  2 |       4  5
    !r   |
    !  3 |          6
    !    |
    !  4 |
    !
    !  index(r, c) = (r - 1) * n + c - (r + 1) * r / 2
    !  begin_index = index(r, r + 1)
    !  where r is equivalent to tag

    offset = [stnum - 1, begin_index - 1]
    block = [1_hsize_t, size_energy]
    stride = [1, 1]
    count = [1, 1]
    ! select data slab in filespace
    call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, err, stride, block) 

    ! create memory space
    call H5Screate_simple_f(1, [int(size_energy, kind=hsize_t)], memspace, err)

    ! Create property list for collective dataset write
    call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
    call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)

    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, energy, size_energy, err, &
                    file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    call H5Pclose_f(plist_id, err)
    call H5Sclose_f(memspace, err)
  end subroutine engtraj_write

  subroutine engtraj_finish()
    implicit none
    integer :: err
    
    call H5Sclose_f(filespace, err)
    call H5Dclose_f(dset_id, err)
    call H5Fclose_f(file_id, err)

    ! this call should be put at the end of the program when HDF5 is no longer needed
    ! putting it here results in a significant performance impact
    !call H5close_f(err)
  end subroutine engtraj_finish
end module engtraj
