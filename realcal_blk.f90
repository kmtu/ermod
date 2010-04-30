module realcal_blk
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
  real :: laxes(3)
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
    use engmain, only: numsite, specatm
    integer, intent(in) :: solu
    integer :: i
    do i = 1, numsite(solu)
       atomno_solu(i) = specatm(i, solu)
       belong_solu(i) = solu
    end do
  end subroutine set_solu_atoms

  subroutine set_solv_atoms(solu, tagpt, slvmax)
    use engmain, only: specatm, numsite
    integer, intent(in) :: solu, tagpt(:), slvmax
    integer :: i, j, k, cnt
    cnt = 1
    do i = 1, slvmax
       j = tagpt(i)
       if(j == solu) cycle
       do k = 1, numsite(j)
          atomno_solv(cnt) = specatm(k, j)
          belong_solv(cnt) = i
          cnt = cnt + 1
       end do
    end do
  end subroutine set_solv_atoms

  subroutine set_block_info()
    use engmain, only: cell, block_threshold, upljcut, elecut
    real :: unit_axes(3), cut2, l
    integer :: i, j, k, bmax, ix
    real, allocatable :: grid_dist(:, :)
    real, allocatable :: box_dist(:, :)
    
    if(abs(cell(2, 1)) > 1e-8 .or. &
       abs(cell(3, 1)) > 1e-8 .or. &
       abs(cell(3, 2)) > 1e-8) then
       print *, cell
       stop "realcal%set_box_info: assertion failed"
    endif

    ! get the length of axes
    ! assumes cell's 1st axis being x-axis, 2nd axis on x-y plane
    do i = 1, 3
       laxes(i) = cell(i, i)
    end do

    ! set block size
    block_size(:) = ceiling(laxes(:) / block_threshold)
    if(any(block_size(:) < 0)) stop "realcal%set_box_info: assertion failed (blocksize)"

    ! pre-calculate grid-grid distance and box-box distance
    bmax = maxval(block_size(:))
    allocate(grid_dist(0:bmax-1, 3))
    allocate(box_dist(0:bmax-1, 3))
    unit_axes(:) = laxes(:) / block_size(:)

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
    use engmain, only: sitepos, specatm
    integer, intent(in) :: natom, atomlist(:)
    integer, intent(out) :: blk(:, :)
    integer :: i, j, a, blktmp(3)

    do i = 1, natom
       a = atomlist(i)
       blktmp(:) = int(floor(sitepos(:, a) / laxes(:) * block_size(:)))
       do j = 1, 3
          blk(j, i) = modulo(blktmp(j), block_size(j))
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
    integer :: ubs(3), vbs(3)
    integer :: i

    do u3 = 0, block_size(3) - 1
       do u2 = 0, block_size(2) - 1
          do u1 = 0, block_size(1) - 1
             do i = 1, subcell_num_neighbour
                vbs(1) = modulo(u1 + subcell_neighbour(1, i) , block_size(1))
                vbs(2) = modulo(u2 + subcell_neighbour(2, i) , block_size(2))
                vbs(3) = modulo(u3 + subcell_neighbour(3, i) , block_size(3))
                
                ubs(:) = (/ u1, u2, u3 /)
                
                call get_pair_energy_block(ubs, vbs, energy_mat)
             end do
          end do
       end do
    end do
  end subroutine get_pair_energy

  subroutine get_pair_energy_block(ubs, vbs, energy_mat)
    use engmain, only: sitepos, cltype, boxshp, upljcut, lwljcut, elecut, ljene, ljlen, cmbrule, screen, charge
    integer, intent(in) :: ubs(3), vbs(3)
    real, intent(out) :: energy_mat(:, :)
    integer :: upos, vpos
    integer :: ui, vi, ua, va
    integer :: belong_u, belong_v
    real :: crdu(3), crdv(3), d(3), invbox(3), dist, r
    real :: elj, eel, rtp1, rtp2, chr2, swth, ljeps, ljsgm
    
    if(cltype == 0) stop "realcal%get_pair_energy_block: cltype assertion failure"
    if(boxshp == 0) stop "realcal%get_pair_energy_block: boxshp assertion failure"

    invbox(:) = 1.0 / laxes(:)

    upos = ubs(1) + block_size(1) * (ubs(2) + block_size(2) * ubs(3))
    vpos = vbs(1) + block_size(1) * (vbs(2) + block_size(2) * vbs(3))
    
    do ui = psum_solu(upos), psum_solu(upos + 1) - 1
       ua = atomno_solu(ui)
       belong_u = belong_solu(ui)
       crdu(:) = sitepos(:, ua)
       ! if(belong_u /= 1) then
       !    print  *, upos, ui, ua, belong_u, crdu
       !    stop "INVALID U"
       ! end if
       do vi = psum_solv(vpos), psum_solv(vpos + 1) - 1
          va = atomno_solv(vi)
          belong_v = belong_solv(vi)
          crdv(:) = sitepos(:, va)
          ! if(belong_v <= 1 .or. belong_v > 14500) then
          !    print  *, vpos, vi, va, belong_v, crdv, psum_solv(vpos) + 1, psum_solv(vpos + 1)
          !    stop "INVALID V"
          ! end if

          d(:) = crdv(:) - crdu(:)
          d(:) = d(:) - laxes(:) * anint(invbox(:) * d(:)) ! get nearest image
          ! FIXME:
          ! assumes that only a single image matters for both electrostatic and LJ.
          ! if the box is very small and strongly anisotropic,
          ! there is a risk that second nearest image still being inside the cutoff length.
          ! But it's not considered in this case ...
          
          dist = sum(d(:) ** 2) ! CHECK: any sane compiler will expand and unroll
          r = sqrt(dist)
          if(r > upljcut) then
             elj = 0.0
          else
             ljeps = sqrt(ljene(va)*ljene(ua))
             select case (cmbrule)
                case (0)
                   ljsgm = (ljlen(va) + ljlen(ua)) * 0.5
                case (1)
                   ljsgm = sqrt(ljlen(va) * ljlen(ua))
                case default
                   stop "Unknown coulomb type @ realcal%pairenergy"
             end select
             rtp1 = ljsgm * ljsgm / dist
             rtp2 = rtp1 * rtp1 * rtp1
             elj = 4.0e0 * ljeps * rtp2 * (rtp2 - 1.0e0)
             if(r > lwljcut) then    ! CHARMM form of switching function
                rtp1 = lwljcut * lwljcut
                rtp2 = upljcut * upljcut
                swth = (2.0e0 * dist + rtp2 - 3.0e0 * rtp1) * (dist - rtp2) * (dist - rtp2) &
                     / ((rtp2 - rtp1) * (rtp2 - rtp1) * (rtp2 - rtp1))
                elj = swth * elj
             endif
          end if
          if(r > elecut) then
             eel = 0.0
          else
             chr2 = charge(ua) * charge(va)
             eel = chr2 * (1.0e0 - erf(screen * r)) / r 
          end if
          energy_mat(belong_v, belong_u) = energy_mat(belong_v, belong_u) + elj + eel
       end do
    end do
  end subroutine get_pair_energy_block
end module realcal_blk
