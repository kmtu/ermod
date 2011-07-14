! -*- F90 -*-
! interfacing module for vmdfio
module vmdfio_interface
  implicit none
contains
  
  subroutine init_vmdplugins
    external vmdfio_init_traj
    call vmdfio_init_traj
  end subroutine init_vmdplugins

  

  subroutine finish_vmdplugins
    external vmdfio_fini_traj
    call vmdfio_fini_traj
  end subroutine finish_vmdplugins

end module vmdfio_interface



