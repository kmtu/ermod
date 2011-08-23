module realcal
  implicit none
  integer :: nsolu_atom, nsolv_atom
  integer, allocatable :: block_solu(:, :), block_solv(:, :)
  integer, allocatable :: belong_solu(:), belong_solv(:)
  integer, allocatable :: atomno_solu(:), atomno_solv(:)
  integer, allocatable :: counts_solu(:, :, :), counts_solv(:, :, :)
  integer, allocatable :: psum_solu(:), psum_solv(:)

  integer, allocatable :: subcell_neighbour(:, :) ! of (3, subcell_num_neighbour)
  integer :: subcell_num_neighbour

  integer :: block_size(3)
  
  ! "straight" coordinate system
  real, allocatable :: sitepos_normal(:, :)
  real :: cell_normal(3, 3), invcell_normal(3), cell_len_normal(3)

contains
  subroutine realcal_proc(target_solu, tagpt, slvmax, uvengy)
    use engmain, only: numsite
    integer, intent(in) :: target_solu, tagpt(:), slvmax
    real, intent(out) :: uvengy(0:slvmax)
    real, allocatable :: eng(:, :)
    
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

    ! assertion
    ! if (.not. all(belong_solu(:) == target_solu)) stop "realcal_blk: target_solu bugged after sorting"

    allocate(eng(1:slvmax, 1))
    eng(:, :) = 0
    call get_pair_energy(eng)

    uvengy(1:slvmax) = eng(1:slvmax, 1)

    deallocate(eng)
    deallocate(block_solu, belong_solu, atomno_solu, counts_solu, psum_solu)
    deallocate(block_solv, belong_solv, atomno_solv, counts_solv, psum_solv)
    deallocate(subcell_neighbour)
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
         charge,ljene,ljlen,specatm,sitepos,&
         cell,invcl,volume,pi,&
         SYS_NONPERIODIC, SYS_PERIODIC
    implicit none
    integer i,j,is,js,ismax,jsmax,ati,atj,m,k
    real reelcut,pairep,ljeps,ljsgm,chr2,rst,dis2,rtp1,rtp2
    real :: eplj,epcl,xst(3),clm(3),swth, half_cell(3)
    real, parameter :: infty=1.0e50      ! essentially equal to infinity
    !
    if(i.eq.j) then
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
          xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
          if(boxshp == SYS_PERIODIC) then              ! when the system is periodic
             xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))
          endif
          dis2=xst(1)*xst(1)+xst(2)*xst(2)+xst(3)*xst(3)
          rst=sqrt(dis2)
          if(rst > upljcut) then
             eplj=0.0e0
          else
             ljeps=ljene(ati) * ljene(atj)

             if(cmbrule.eq.0) ljsgm=(ljlen(ati)+ljlen(atj))/2.0e0
             if(cmbrule.eq.1) ljsgm=sqrt(ljlen(ati)*ljlen(atj))

             rtp1=ljsgm*ljsgm/dis2
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

    call residual_ene(i, i, pairep)
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
    integer :: i, j, k, bmax, ix
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
          do i = 0, block_size(1) - 1
             l = box_dist(i, 1) ** 2 + box_dist(j, 2) ** 2 + box_dist(k, 3) ** 2
             if(l < cut2) subcell_num_neighbour = subcell_num_neighbour + 1
          end do
       end do
    end do
    allocate(subcell_neighbour(3, subcell_num_neighbour))
    ! then generate the list
    ix = 1
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          do i = 0, block_size(1) - 1
             l = box_dist(i, 1) ** 2 + box_dist(j, 2) ** 2 + box_dist(k, 3) ** 2
             if(l < cut2) then
                subcell_neighbour(1, ix) = i
                subcell_neighbour(2, ix) = j
                subcell_neighbour(3, ix) = k
                ix = ix + 1
             endif
          end do
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

  ! FIXME: create pairenergy_single_solu as specilization?
  subroutine get_pair_energy(energy_mat)
    ! calculate for each subcell
    ! cut-off by subcell distance
    real, intent(out) :: energy_mat(:, :)
    integer :: u1, u2, u3
    integer :: vbs(3)
    integer :: i, upos, vpos

    do u3 = 0, block_size(3) - 1
       do u2 = 0, block_size(2) - 1
          do u1 = 0, block_size(1) - 1
             upos = u1 + block_size(1) * (u2 + block_size(2) * u3)
             if(psum_solu(upos + 1) /= psum_solu(upos)) then ! if solute have atoms in the block
                do i = 1, subcell_num_neighbour
                   vbs(1) = mod(u1 + subcell_neighbour(1, i) , block_size(1))
                   vbs(2) = mod(u2 + subcell_neighbour(2, i) , block_size(2))
                   vbs(3) = mod(u3 + subcell_neighbour(3, i) , block_size(3))
                   
                   vpos = vbs(1) + block_size(1) * (vbs(2) + block_size(2) * vbs(3))

                   call get_pair_energy_block(upos, vpos, energy_mat)
                end do
             end if
          end do
       end do
    end do
  end subroutine get_pair_energy

  subroutine get_pair_energy_block(upos, vpos, energy_mat)
    use engmain, only: cltype, boxshp, upljcut, lwljcut, elecut, ljene, ljlen, cmbrule, screen, charge
    implicit none
    integer, intent(in) :: upos, vpos
    real, intent(out) :: energy_mat(:, :)
    integer :: ui, vi, ua, va
    integer :: belong_u, belong_v
    real :: crdu(3), crdv(3), d(3), dist, r
    real :: elj, eel, rtp1, rtp2, chr2, swth, ljeps, ljsgm2
    real :: upljcut2, lwljcut2, elecut2, half_cell(3)
    
    if(cltype == 0) stop "realcal%get_pair_energy_block: cltype assertion failure"
    if(boxshp == 0) stop "realcal%get_pair_energy_block: boxshp assertion failure"

    half_cell(:) = 0.5 * cell_len_normal(:)

    lwljcut2 = lwljcut ** 2
    upljcut2 = upljcut ** 2
    elecut2 = elecut ** 2

    do ui = psum_solu(upos), psum_solu(upos + 1) - 1
       ua = atomno_solu(ui)
       belong_u = belong_solu(ui) ! FIXME: not used in later calculation
       crdu(:) = sitepos_normal(:, ua)
       do vi = psum_solv(vpos), psum_solv(vpos + 1) - 1
          va = atomno_solv(vi)
          belong_v = belong_solv(vi)
          crdv(:) = sitepos_normal(:, va)

          d(:) = crdv(:) - crdu(:)
          d(:) = half_cell(:) - abs(half_cell(:) - abs(d(:))) ! get nearest image
          ! FIXME:
          ! assumes that only a single image matters for both electrostatic and LJ.
          ! if the box is very small and strongly anisotropic,
          ! there is a risk that second nearest image still being inside the cutoff length.
          ! But it's not considered in this case ...
          
          dist = sum(d(:) ** 2) ! CHECK: any sane compiler will expand and unroll
          r = sqrt(dist)
          if(dist >= upljcut2) then
             elj = 0.0
          else
             ljeps = ljene(va) * ljene(ua)
             select case (cmbrule)
                case (0)
                   ljsgm2 = ((ljlen(va) + ljlen(ua)) * 0.5) ** 2
                case (1)
                   ljsgm2 = ljlen(va) * ljlen(ua)
                case default
                   stop "Unknown coulomb type @ realcal%pairenergy"
             end select
             rtp1 = ljsgm2 / dist
             rtp2 = rtp1 * rtp1 * rtp1
             elj = 4.0e0 * ljeps * rtp2 * (rtp2 - 1.0e0)
             if(dist > lwljcut2) then    ! CHARMM form of switching function
                rtp1 = lwljcut2
                rtp2 = upljcut2
                swth = (2.0e0 * dist + rtp2 - 3.0e0 * rtp1) * (dist - rtp2) * (dist - rtp2) &
                     / ((rtp2 - rtp1) * (rtp2 - rtp1) * (rtp2 - rtp1))
                elj = swth * elj
             endif
          end if
          if(dist > elecut2) then
             eel = 0.0
          else
             chr2 = charge(ua) * charge(va)
             eel = chr2 * (1.0e0 - erf(screen * r)) / r 
          end if
          energy_mat(belong_v, 1) = energy_mat(belong_v, 1) + elj + eel
       end do
    end do
  end subroutine get_pair_energy_block

  ! Calculate residual energy in Ewald / PME
  subroutine residual_ene(i, j, pairep)
    use engmain, only: screen, volume, numsite, mol_charge, cltype
    implicit none
    integer, intent(in) :: i, j
    real, intent(inout) :: pairep
    real :: rtp1, rtp2, epcl
    integer :: is, js, ismax, jsmax, ati, atj
    real, parameter :: pi = 3.141592653589793283462
    if(cltype == 0) stop "Error: residual_ene: called when cltype == 0, cannot happen"

    rtp1 = mol_charge(i)
    rtp2 = mol_charge(j)
    epcl=pi*rtp1*rtp2/screen/screen/volume
    if(i.eq.j) epcl=epcl/2.0e0 ! self-interaction
    pairep=pairep-epcl
  end subroutine residual_ene

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
       cell_len_normal(i) = cell_normal(i, i)
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
