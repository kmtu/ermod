! -*- F90 -*-
! DO NOT rewrite realcal.F90,
! realcal.F90 is generated from realcal.F90.php.
module realcal
  implicit none
  integer :: nsolu_atom, nsolv_atom
  integer, allocatable :: block_solu(:, :), block_solv(:, :)
  integer, allocatable :: belong_solu(:), belong_solv(:)
  integer, allocatable :: atomno_solu(:), atomno_solv(:)
  integer, allocatable :: counts_solu(:, :, :), counts_solv(:, :, :)
  integer, allocatable :: psum_solu(:), psum_solv(:)
  real, allocatable :: sitepos_solu(:, :), sitepos_solv(:, :)
  integer, allocatable :: ljtype_solu(:, :), ljtype_solv(:, :)

  integer :: max_solu_block, max_solv_block

  real, allocatable :: ljeps_lowlj(:, :), ljsgm2_lowlj(:, :), dist_lowlj(:, :)
  integer, allocatable :: belong_lowlj(:, :)
  real, allocatable :: ljeps_switch(:, :), ljsgm2_switch(:, :), dist_switch(:, :)
  integer, allocatable :: belong_switch(:, :)
  real, allocatable :: charge_el(:, :), dist_el(:, :)
  integer, allocatable :: belong_el(:, :)
  real, allocatable :: e_t(:, :)

  ! subcell_neighbour only stores the neighbour list on y-z plane. x direction is stored on subcell_xlen.
  ! (looks like sub"pillar" rather than cell)
  ! each subcell spans x=(-subcell_xlen(ix) .. 1+subcell_xlen(ix)),  y=subcell_neighbour(2, ix), z=subcell_neighbour(3, ix)
  integer, allocatable :: subcell_neighbour(:, :) ! of (2:3, subcell_num_neighbour)
  integer, allocatable :: subcell_xlen(:)
  integer :: subcell_num_neighbour

  integer :: block_size(3)
  
  ! "straight" coordinate system
  real, allocatable :: sitepos_normal(:, :)
  real :: cell_normal(3, 3), invcell_normal(3), cell_len_normal(3)

contains
  subroutine realcal_proc(target_solu, tagpt, slvmax, uvengy)
    use engmain, only: numsite
    !$ use omp_lib, only: omp_get_num_threads
    integer, intent(in) :: target_solu, tagpt(:), slvmax
    real, intent(out) :: uvengy(0:slvmax)
    real, allocatable :: eng(:, :)
    integer :: lsize, i, j
    integer :: npar

    npar = 1
    !$ npar = omp_get_num_threads()
    
    ! print *, "DEBUG: relcal_proc called"
    ! FIXME: fix calling convention & upstream call tree
    ! to calculate several solutes at once
    nsolu_atom = numsite(target_solu)
    nsolv_atom = count_solv(target_solu, tagpt, slvmax)

    call set_block_info()

    allocate(block_solu(3, nsolu_atom), block_solv(3, nsolv_atom))
    allocate(belong_solu(nsolu_atom), belong_solv(nsolv_atom))
    allocate(atomno_solu(nsolu_atom), atomno_solv(nsolv_atom))
    allocate(counts_solu(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))
    allocate(counts_solv(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))
    allocate(psum_solu(0:block_size(1) * block_size(2) *  block_size(3)))
    allocate(psum_solv(0:block_size(1) * block_size(2) *  block_size(3)))

    call set_solv_atoms(target_solu, tagpt, slvmax)
    call set_solu_atoms(target_solu)

    ! assertion
    ! if (.not. all(belong_solu(:) == target_solu)) stop "realcal_blk: target_solu bugged"

    call blockify(nsolu_atom, atomno_solu, block_solu)
    call blockify(nsolv_atom, atomno_solv, block_solv)

    call sort_block(block_solu, nsolu_atom, belong_solu, atomno_solu, counts_solu, psum_solu)
    call sort_block(block_solv, nsolv_atom, belong_solv, atomno_solv, counts_solv, psum_solv)

    allocate(sitepos_solv(3, nsolv_atom + 1))
    sitepos_solv(1:3, 1:nsolv_atom) = sitepos_normal(1:3, atomno_solv(1:nsolv_atom))
    sitepos_solv(1:3, nsolv_atom+1) = 0

    max_solu_block = maxval(counts_solu)
    max_solv_block = 0
    do i = 0, block_size(2) - 1
       do j = 0, block_size(3) - 1
          max_solv_block = max(max_solv_block, sum(counts_solv(:, i, j)))
       end do
    end do
    lsize = max_solu_block * max_solv_block
    
    allocate(ljeps_lowlj(lsize, npar), ljsgm2_lowlj(lsize, npar), dist_lowlj(lsize, npar), belong_lowlj(lsize, npar))
    allocate(ljeps_switch(lsize, npar), ljsgm2_switch(lsize, npar), dist_switch(lsize, npar), belong_switch(lsize, npar))
    allocate(charge_el(lsize, npar), dist_el(lsize, npar), belong_el(lsize, npar))
    allocate(e_t(lsize, npar))
    ! assertion
    ! if (.not. all(belong_solu(:) == target_solu)) stop "realcal_blk: target_solu bugged after sorting"

    allocate(eng(1:slvmax, npar))
    eng(:, :) = 0
    call get_pair_energy(eng)

    uvengy(1:slvmax) = sum(eng(1:slvmax, 1:npar), 2)

    deallocate(ljeps_lowlj, ljsgm2_lowlj, dist_lowlj, belong_lowlj)
    deallocate(ljeps_switch, ljsgm2_switch, dist_switch, belong_switch)
    deallocate(charge_el, dist_el, belong_el)
    deallocate(e_t)

    deallocate(eng)
    deallocate(block_solu, belong_solu, atomno_solu, counts_solu, psum_solu)
    deallocate(block_solv, belong_solv, atomno_solv, counts_solv, psum_solv)
    deallocate(sitepos_solv)
    deallocate(subcell_neighbour, subcell_xlen)
  end subroutine realcal_proc

  subroutine realcal_prepare
    use engmain, only: numatm, sitepos, boxshp, SYS_PERIODIC
    implicit none
    
    allocate(sitepos_normal(3, numatm))
    sitepos_normal(:, :) = sitepos(:, :)

    ! "Straighten" box, and normalize coordinate system
    if(boxshp == SYS_PERIODIC) call normalize_periodic
  end subroutine realcal_prepare

  subroutine realcal_cleanup
    implicit none

    deallocate(sitepos_normal)
  end subroutine realcal_cleanup


  ! Calculate i-j interaction energy.
  ! This routine is called as a dispatcher to realcal_self or bare coulomb interaction
  subroutine realcal_bare(i,j,pairep)
    use engmain, only:  nummol,maxsite,numatm,boxshp,numsite,&
         elecut,lwljcut,upljcut,cmbrule,cltype,screen,&
         charge,specatm,sitepos,&
         cell,invcl,volume,pi,&
         ljtype, ljtype_max, ljene_mat, ljlensq_mat,&
         SYS_NONPERIODIC, SYS_PERIODIC
    implicit none
    integer i,j,is,js,ismax,jsmax,ati,atj,m,k
    real reelcut,pairep,ljeps,ljsgm,chr2,rst,dis2,rtp1,rtp2
    real :: eplj,epcl,xst(3),clm(3),swth, half_cell(3)
    integer :: ljtype_i, ljtype_j
    real, parameter :: infty=1.0e50      ! essentially equal to infinity
    !
    if(i.eq.j) then
       stop
       call realcal_self(i, pairep)
       return
    endif

    if(cltype /= 0) stop "cannot happen: realcal() is called only when cltype is 'bare coulomb'."

    if(boxshp == SYS_NONPERIODIC) reelcut=infty
    if(boxshp == SYS_PERIODIC) then
       reelcut=elecut
       half_cell(:) = 0.5 * cell_len_normal(:)
    endif

    pairep=0.0e0
    ismax=numsite(i)
    jsmax=numsite(j)

    do is=1,ismax
       do js=1,jsmax
          ati=specatm(is,i)
          atj=specatm(js,j)
          ljtype_i = ljtype(ati)
          ljtype_j = ljtype(atj)
          xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
          if(boxshp == SYS_PERIODIC) then              ! when the system is periodic
             xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))
          endif
          dis2=xst(1)*xst(1)+xst(2)*xst(2)+xst(3)*xst(3)
          rst=sqrt(dis2)
          if(rst > upljcut) then
             eplj=0.0e0
          else
             ljeps=ljene_mat(ljtype_i, ljtype_j)
             ljsgm=ljlensq_mat(ljtype_i, ljtype_j)

             rtp1=ljsgm/dis2
             rtp2=rtp1*rtp1*rtp1
             eplj=4.0e0*ljeps*rtp2*(rtp2-1.0e0)
             if(rst > lwljcut) then    ! CHARMM form of switching function
                rtp1=lwljcut*lwljcut
                rtp2=upljcut*upljcut
                swth=(2.0e0*dis2+rtp2-3.0e0*rtp1)*(dis2-rtp2)*(dis2-rtp2)&
                     /(rtp2-rtp1)/(rtp2-rtp1)/(rtp2-rtp1)
                eplj=swth*eplj
             endif
          endif
          if(rst >= reelcut) then
             epcl=0.0e0
          else
             chr2=charge(ati)*charge(atj)

             epcl=chr2 / rst
          endif
          pairep = pairep + (eplj + epcl)
       end do
    end do
    !
    return
  end subroutine realcal_bare

  subroutine realcal_self(i, pairep)
    use engmain, only:  nummol,maxsite,numatm,boxshp,numsite,&
         screen,cltype,&
         charge,specatm,sitepos,&
         cell,invcl,EL_COULOMB, pi
    implicit none
    integer, intent(in) :: i
    real, intent(inout) :: pairep
    integer :: is,js,ismax,ati,atj,m,k
    real :: reelcut,chr2,rst,dis2,rtp1
    real :: epcl,xst(3),clm(3),swth, half_cell(3)

    pairep=0.0e0
    if(cltype == EL_COULOMB) return

    half_cell(:) = 0.5 * cell_len_normal(:) 

    ismax=numsite(i)

    do is=1,ismax
       ati=specatm(is,i)

       ! Atom residual
       chr2=charge(ati)*charge(ati)
       epcl=-chr2*screen/sqrt(pi)

       pairep = pairep + epcl

       do js=is+1,ismax
          atj=specatm(js,i)
 
          xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
          xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))

          dis2=xst(1)*xst(1)+xst(2)*xst(2)+xst(3)*xst(3)

          rst=sqrt(dis2)
          chr2=charge(ati)*charge(atj)
          epcl=-chr2*derf(screen*rst)/rst

          pairep=pairep+epcl
       enddo
    enddo

  end subroutine realcal_self


  integer function count_solv(solu, tagpt, slvmax)
    use engmain, only: numsite
    integer, intent(in) :: solu, tagpt(:), slvmax
     integer :: i, j, cnt
    cnt = 0
    do i = 1, slvmax
       j = tagpt(i)
       if(j == solu) cycle
       cnt = cnt + numsite(j)
    end do
    count_solv = cnt
  end function count_solv

  subroutine set_solu_atoms(solu)
    use engmain, only: numsite, mol_begin_index
    integer, intent(in) :: solu
    integer :: i
    do i = 1, numsite(solu)
       atomno_solu(i) = mol_begin_index(solu) + (i - 1)
       belong_solu(i) = solu
    end do
  end subroutine set_solu_atoms

  subroutine set_solv_atoms(solu, tagpt, slvmax)
    use engmain, only: numsite, mol_begin_index
    integer, intent(in) :: solu, tagpt(:), slvmax
    integer :: i, j, k, cnt
    cnt = 1
    do i = 1, slvmax
       j = tagpt(i)
       if(j == solu) cycle
       do k = 1, numsite(j)
          atomno_solv(cnt) = mol_begin_index(j) + (k - 1)
          belong_solv(cnt) = i
          cnt = cnt + 1
       end do
    end do
  end subroutine set_solv_atoms

  subroutine set_block_info()
    use engmain, only: block_threshold, upljcut, elecut
    real :: unit_axes(3), cut2, l
    integer :: i, j, k, bmax, ix, xlen
    real, allocatable :: grid_dist(:, :)
    real, allocatable :: box_dist(:, :)
    
    ! get the length of axes
    ! assumes cell's 1st axis being x-axis, 2nd axis on x-y plane

    ! set block size
    block_size(:) = ceiling(cell_len_normal(:) / block_threshold)

    ! pre-calculate grid-grid distance and box-box distance
    bmax = maxval(block_size(:))
    allocate(grid_dist(0:bmax-1, 3))
    allocate(box_dist(0:bmax-1, 3))
    unit_axes(:) = cell_len_normal(:) / block_size(:)

    ! only have to calculate axis-wise distance (because of orthogonality)
    do j = 1, 3
       do i = 0, block_size(j) - 1
          grid_dist(i, j) = min(i, modulo(-i, block_size(j))) * unit_axes(j)
       end do
    end do
    do j = 1, 3
       do i = 0, block_size(j) - 1
          box_dist(i, j) = min(grid_dist(i, j),&
               grid_dist(modulo(i - 1, block_size(j)), j), &
               grid_dist(modulo(i + 1, block_size(j)), j))
       end do
    end do

    deallocate(grid_dist)
    
    ! make a subcell list
    cut2 = max(upljcut, elecut) ** 2
    ! count...
    subcell_num_neighbour = 0
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          l = box_dist(j, 2) ** 2 + box_dist(k, 3) ** 2
          if(l < cut2) subcell_num_neighbour = subcell_num_neighbour + 1
       end do
    end do
    allocate(subcell_neighbour(2:3, subcell_num_neighbour))
    allocate(subcell_xlen(subcell_num_neighbour))

    ! then generate the list
    ix = 1
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          l = box_dist(j, 2) ** 2 + box_dist(k, 3) ** 2
          if(l < cut2) then
             xlen = ceiling(sqrt(cut2 - l) / unit_axes(1))
             subcell_neighbour(2, ix) = j
             subcell_neighbour(3, ix) = k
             subcell_xlen(ix) = xlen ! this subpillar contains subcell [-xlen .. xlen]; note the end of region is inclusive
             ix = ix + 1
          endif
       end do
    end do
    deallocate(box_dist)
  end subroutine set_block_info

  subroutine blockify(natom, atomlist, blk)
    integer, intent(in) :: natom, atomlist(:)
    integer, intent(out) :: blk(:, :)
    integer :: i, j, a, blktmp(3)

    do i = 1, natom
       a = atomlist(i)
       blktmp(:) = int(floor(sitepos_normal(:, a) * invcell_normal(:) * block_size(:)))
       do j = 1, 3
          blk(j, i) = modulo(blktmp(j), block_size(j))
          if(blk(j,i) < 0) then
             print *, cell_len_normal(:), blktmp(:), block_size(:)
             STOP "INVLBLK"
          endif
       end do
    end do
  end subroutine blockify

  subroutine sort_block(blk, nmol, belong, atomno, counts, psum)
    integer, intent(inout) :: blk(:, :)
    integer, intent(in) :: nmol
    integer, intent(inout) :: belong(:)
    integer, intent(inout) :: atomno(:)
    integer, intent(inout) :: counts(0:block_size(1) - 1, 0:block_size(2) - 1, 0:block_size(3) - 1)
    integer, intent(out) :: psum(0:block_size(1) * block_size(2) * block_size(3))
    integer, allocatable :: buffer(:, :) ! FIXME: ugly!
    integer, allocatable :: pnum(:, :, :)
    integer :: a, b, c
    integer :: i, j, k, partialsum, pos
    
    counts(:, :, :) = 0

    do i = 1, nmol
       a = blk(1, i)
       b = blk(2, i)
       c = blk(3, i)
       if(a < 0 .or. b < 0 .or. c < 0) STOP "INVL"
       counts(a, b, c) = counts(a, b, c) + 1
    end do

    allocate(pnum(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))

    partialsum = 0
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          do i = 0, block_size(1) - 1
             pnum(i, j, k) = partialsum
             partialsum = partialsum + counts(i, j, k)
          end do
       end do
    end do

    allocate(buffer(5, nmol))
    do i = 1, nmol
       a = blk(1, i)
       b = blk(2, i)
       c = blk(3, i)
       pos = pnum(a, b, c) + 1  ! buffer (and output) index starts from 1
       pnum(a, b, c) = pos
       buffer(1:3, pos) = blk(1:3, i)
       buffer(4, pos) = belong(i)
       buffer(5, pos) = atomno(i)
    end do
    deallocate(pnum)

    blk(1:3, :) = buffer(1:3, :)
    belong(:) = buffer(4, :)
    atomno(:) = buffer(5, :)

    deallocate(buffer)

    partialsum = 0
    pos = 0
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          do i = 0, block_size(1) - 1
             psum(pos) = partialsum + 1 ! output index starts from 1
             pos = pos + 1
             partialsum = partialsum + counts(i, j, k)
          end do
       end do
    end do
    psum(pos) = partialsum + 1
    ! print *, psum
  end subroutine sort_block

  subroutine get_pair_energy(energy_vec)
    ! calculate for each subcell
    ! cut-off by subcell distance
    real, intent(out) :: energy_vec(:, :)
    integer :: u1, u2, u3
    integer :: vbs(3)
    integer :: i, upos, vpos_base, vpos_line_end, vpos_begin, vpos_end
    integer :: xlen

    do u3 = 0, block_size(3) - 1
       do u2 = 0, block_size(2) - 1
          do u1 = 0, block_size(1) - 1
             !$omp task
             upos = u1 + block_size(1) * (u2 + block_size(2) * u3)
             if(psum_solu(upos + 1) /= psum_solu(upos)) then ! if solute have atoms in the block
                do i = 1, subcell_num_neighbour
                   vbs(2) = mod(u2 + subcell_neighbour(2, i) , block_size(2))
                   vbs(3) = mod(u3 + subcell_neighbour(3, i) , block_size(3))
                   vpos_base = block_size(1) * (vbs(2) + block_size(2) * vbs(3))
                   
                   xlen = subcell_xlen(i)
                   if(2 * xlen + 1 >= block_size(1)) then
                      ! spans all over x-axis
                      call get_pair_energy_block(upos, vpos_base, vpos_base + block_size(1), energy_vec)
                   else
                      vpos_begin = vpos_base + u1 - xlen
                      vpos_end = vpos_base + u1 + xlen + 1
                      vpos_line_end = vpos_base + block_size(1)

                      if(vpos_begin < vpos_base) then
                         ! spans over periodic boundary, case 1
                         call get_pair_energy_block(upos, vpos_begin + block_size(1), vpos_line_end, energy_vec)
                         call get_pair_energy_block(upos, vpos_base, vpos_end, energy_vec)
                      elseif(vpos_end > vpos_line_end) then
                         ! spans over periodic boundary, case 2
                         call get_pair_energy_block(upos, vpos_begin, vpos_line_end, energy_vec)
                         call get_pair_energy_block(upos, vpos_base, vpos_end - block_size(1), energy_vec)
                      else
                         ! standard case
                         call get_pair_energy_block(upos, vpos_begin, vpos_end, energy_vec)
                      endif
                   endif
                end do
             end if
             !$omp end task
          end do
       end do
    end do
    !$omp taskwait
  end subroutine get_pair_energy

  ! Computational kernel to calculate distance between particles
  subroutine get_pair_energy_block(upos, vpos_b, vpos_e, energy_vec)
    use engmain, only: cltype, boxshp, upljcut, lwljcut, elecut, screen, charge,&
         ljtype, ljtype_max, ljene_mat, ljlensq_mat
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    integer, intent(in) :: upos, vpos_b, vpos_e
    real, intent(out) :: energy_vec(:, :)
    integer :: ui, vi, ua, va
    integer :: belong_u, belong_v, ljtype_u, ljtype_v
    real :: crdu(3), crdv(3), d(3), dist, r, dist_next
    real :: elj, eel, rtp1, rtp2, chr2, swth, ljeps, ljsgm2
    real :: upljcut2, lwljcut2, elecut2, half_cell(3)
    integer :: n_lowlj, n_switch, n_el
    integer :: i, curp
    
    if(cltype == 0) stop "realcal%get_pair_energy_block: cltype assertion failure"
    if(boxshp == 0) stop "realcal%get_pair_energy_block: boxshp assertion failure"

    curp = 1
    !$ curp = omp_get_thread_num() + 1

    half_cell(:) = 0.5 * cell_len_normal(:)

    n_lowlj = 0
    n_switch = 0
    n_el = 0

    lwljcut2 = lwljcut ** 2
    upljcut2 = upljcut ** 2
    elecut2 = elecut ** 2

    ! TODO optimize:
    ! if you sort / reorder atomno, ljtype etc.
    ! this loop can be vectorized
    ! TODO optimize:
    ! software pipelining
    do ui = psum_solu(upos), psum_solu(upos + 1) - 1
       ua = atomno_solu(ui)
       belong_u = belong_solu(ui) ! FIXME: not used in later calculation
       ljtype_u = ljtype(ua)
       crdu(:) = sitepos_normal(:, ua)

       ! hide latency by calculating distance of next coordinate set
       d(:) = crdu(:) - sitepos_solv(:, psum_solv(vpos_b))
       d(:) = half_cell(:) - abs(half_cell(:) - abs(d(:)))
       dist_next = sum(d(:) ** 2)
       
       do vi = psum_solv(vpos_b), psum_solv(vpos_e) - 1
          va = atomno_solv(vi)
          belong_v = belong_solv(vi)
          ljtype_v = ljtype(va)
          
          crdv(:) = sitepos_solv(:, vi + 1)

          d(:) = crdv(:) - crdu(:)
          d(:) = half_cell(:) - abs(half_cell(:) - abs(d(:))) ! get nearest image
          ! assumes that only a single image matters for both electrostatic and LJ.
          ! if the box is very small and strongly anisotropic,
          ! there is a risk that second nearest image still being inside the cutoff length.
          ! But it's not considered in this case ...
          
          dist = dist_next
          dist_next = sum(d(:) ** 2) ! CHECK: any sane compiler will expand and unroll

          ! lines up all variables, to enable vectorization in 2nd phase
          if(dist <= lwljcut2) then
             n_lowlj = n_lowlj + 1
             ljeps_lowlj (n_lowlj, curp) = ljene_mat(ljtype_v, ljtype_u)
             ljsgm2_lowlj(n_lowlj, curp) = ljlensq_mat(ljtype_v, ljtype_u)
             dist_lowlj(n_lowlj, curp) = dist
             belong_lowlj(n_lowlj, curp) = belong_v
          elseif(dist <= upljcut2) then
             n_switch = n_switch + 1
             ljeps_switch (n_switch, curp) = ljene_mat(ljtype_v, ljtype_u)
             ljsgm2_switch(n_switch, curp) = ljlensq_mat(ljtype_v, ljtype_u)
             dist_switch(n_switch, curp) = dist
             belong_switch(n_switch, curp) = belong_v
          end if
          
          if(dist <= elecut2) then
             n_el = n_el + 1
             charge_el(n_el, curp) = charge(ua) * charge(va)
             dist_el(n_el, curp) = dist
             belong_el(n_el, curp) = belong_v
          end if
       end do
    end do

    ! 2nd phase: calculate actual values with vectorized loop

    ! TODO optimize:
    ! explicit vectorization may increase performance

    ! LJ inside low cutoff
    do i = 1, n_lowlj
       rtp1 = ljsgm2_lowlj(i, curp) / dist_lowlj(i, curp)
       rtp2 = rtp1 * rtp1 * rtp1
       e_t(i, curp) = 4.0e0 * ljeps_lowlj(i, curp) * rtp2 * (rtp2 - 1.0e0)
    end do
    do i = 1, n_lowlj
       energy_vec(belong_lowlj(i, curp), curp) = energy_vec(belong_lowlj(i, curp), curp) + e_t(i, curp)
    end do

    ! LJ switching region
    do i = 1, n_switch
       rtp1 = ljsgm2_switch(i, curp) / dist_switch(i, curp)
       rtp2 = rtp1 * rtp1 * rtp1
       e_t(i, curp) = 4.0e0 * ljeps_switch(i, curp) * rtp2 * (rtp2 - 1.0e0) * &
            (2.0e0 * dist_switch(i, curp) + upljcut2 - 3.0e0 * lwljcut2) * &
            ((dist_switch(i, curp) - upljcut2) ** 2) * (1.0 / (upljcut2 - lwljcut2) ** 3)
    end do
    do i = 1, n_switch
       energy_vec(belong_switch(i, curp), curp) = energy_vec(belong_switch(i, curp), curp) + e_t(i, curp)
    end do

    ! ewald electrostatic
    do i = 1, n_el
       r = sqrt(dist_el(i, curp))
       e_t(i, curp) = charge_el(i, curp) * (1.0e0 - erf(screen * r)) / r
    end do
    do i = 1, n_el
       energy_vec(belong_el(i, curp), curp) = energy_vec(belong_el(i, curp), curp) + e_t(i, curp)
    end do
  end subroutine get_pair_energy_block

  ! Rotate box and coordinate so that the cell(:,:) is upper triangular:
  ! cell(2, 1) = cell(3, 1) = cell(3, 2) = 0
  ! (1st axis to be aligned to x-axis, 2nd axis to be within xy-plane)
  subroutine normalize_periodic
    use engmain, only: cell, sitepos
    implicit none
    integer :: n
    real :: newcell(3, 3), scale(3), qr(3,3)
    integer :: lwork
    real, allocatable :: work(:)
    integer :: info, perm(3)
    integer :: i
    real :: dummy

    n = size(sitepos, 2)
    lwork = max(3 * 3 + 1, n)
    allocate(work(lwork))

    ! QR-factorize box vector
    qr(:, :) = cell
    perm(:) = 0

    if(kind(dummy) == 8) then
       call dgeqp3(3, 3, qr, 3, perm, scale, work, lwork, info)
    else
       call sgeqp3(3, 3, qr, 3, perm, scale, work, lwork, info)
    endif
    if(info /= 0) stop "engproc.f90, normalize_periodic: failed to factorize box vector"

    ! reorganize R
    ! AP = QR    <=>  Q^T AP = R
    newcell(:, :) = cell(:, perm(:))

    if(kind(dummy) == 8) then
       call dormqr('L', 'T', 3, 3, 3, qr, 3, scale, newcell, 3, work, lwork, info)
    else
       call sormqr('L', 'T', 3, 3, 3, qr, 3, scale, newcell, 3, work, lwork, info)
    endif
    if(info /= 0) stop "engproc.f90, normalize_periodic: failed to rotate cell"

    cell_normal(:, :) = newcell(:, :)
    if(abs(newcell(2, 1)) > 1e-8 .or. &
       abs(newcell(3, 1)) > 1e-8 .or. &
       abs(newcell(3, 2)) > 1e-8) then
       print *, newcell
       stop "engproc.f90, normalize_periodic: assertion failed, box rotation is bugged"
    endif

    ! rotate coordinates
    ! Note: sitepos does not need permutation. (the order of cell axis is not important)
    if(kind(dummy) == 8) then
       call dormqr('L', 'T', 3, n, 3, qr, 3, scale, sitepos_normal, 3, work, lwork, info)
    else
       call sormqr('L', 'T', 3, n, 3, qr, 3, scale, sitepos_normal, 3, work, lwork, info)
    endif
    if(info /= 0) stop "engproc.f90, normalize_periodic: failed to rotate coordinate"

    deallocate(work)

    do i = 1, 3
       cell_len_normal(i) = abs(cell_normal(i, i))
    end do
    invcell_normal(:) = 1 / cell_len_normal(:)

    ! normalize coordinate within single periodicity
    ! move all particles inside the cuboid spanned by (0 .. cell(1, 1)), (0 .. cell(2,2)), (0 .. cell(3,3))
    do i = 1, 3
       sitepos_normal(i, :) = sitepos_normal(i, :) - &
            floor(sitepos_normal(i, :) * invcell_normal(i)) * cell_len_normal(i)
    end do
  end subroutine normalize_periodic

end module realcal
