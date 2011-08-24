
module engproc
  implicit none
  integer :: cntdst, slvmax
  integer :: maxdst
  integer :: tagslt
  integer, allocatable :: tagpt(:)

  ! flceng needs to be output in-order
  logical, allocatable :: flceng_stored(:)
  real, allocatable :: flceng(:, :)

contains
  !
  !  procedure for constructing energy distribution functions
  !
  subroutine enginit
    use engmain, only: numtype,nummol,engdiv,corrcal,slttype,&
         moltype,sluvid,&
         ermax,numslv,uvmax,uvsoft,esmax,uvspec,&
         uvcrd,edens,ecorr,escrd,eself,&
         voffset, &
         aveuv,slnuv,avediv,minuv,maxuv,numslt,sltlist,&
         ene_param, ene_confname, &
         io_flcuv, CAL_SOLN
    use mpiproc, only: halt_with_error, myrank
    implicit none
    real ecdmin,ecfmns,ecmns0,ecdcen,ecpls0,ecfpls,eccore,ecdmax
    real eclbin,ecfbin,ec0bin,finfac,ectmvl
    integer peread,pemax,pesoft,pecore,sltmltp
    character(*), parameter :: ecdfile='EcdInfo'
    integer, parameter :: ecdio=51       ! IO for ecdfile
    real, parameter :: infty=1.0e50      ! essentially equal to infinity
    !
    integer, parameter :: rglmax=5, large=10000
    real, parameter :: tiny=1.0e-30
    integer iduv,i,k,q,pti,regn,prvmx,curmx,uprgcd(rglmax+1)
    real factor,incre,cdrgvl(0:rglmax+1),rgcnt(rglmax),ecpmrd(10)
    integer, dimension(:), allocatable :: tplst
    real, dimension(:,:), allocatable  :: ercrd
    !
    integer, parameter :: paramfile_io=191
    integer :: param_err
    namelist /hist/ ecdmin, ecfmns, ecmns0, ecdcen, ecpls0, ecfpls, eccore, ecdmax,&
         eclbin, ecfbin, ec0bin, finfac, ectmvl,&
         peread, pemax, pesoft, pecore
    !
    allocate( tplst(nummol) )
    numslt=0
    do i=1,nummol
       if(sluvid(i).gt.0) then
          numslt=numslt+1
          tplst(numslt)=i
          sltmltp=moltype(i)
       endif
    end do
    do i=1,nummol
       if((sluvid(i).gt.0).and.(moltype(i).ne.sltmltp)) then
          call halt_with_error('typ')
       endif
    end do
    iduv=0
    if(numslt.le.0) iduv=9
    if((slttype.ge.2).and.(numslt.ne.1)) iduv=9
    if(iduv.ne.0) call halt_with_error('num')
    allocate( sltlist(numslt) )
    do i=1,numslt
       sltlist(i)=tplst(i)
    end do
    if((slttype.ge.2).and.(sltlist(1).ne.nummol)) call halt_with_error('ins')
    deallocate( tplst )
    !
    if(numslt.eq.1) numslv=numtype-1
    if(numslt.gt.1) numslv=numtype
    !
    allocate( uvspec(nummol) )
    do i=1,nummol
       pti=moltype(i)
       if(sluvid(i).eq.0) then            ! solvent
          if(pti.lt.sltmltp) uvspec(i)=pti
          if(pti.eq.sltmltp) call halt_with_error('typ')
          if(pti.gt.sltmltp) uvspec(i)=pti-1
       endif
       if(sluvid(i).ne.0) then            ! solute
          if(pti.ne.sltmltp) call halt_with_error('typ')
          if(numslt.eq.1) uvspec(i)=0
          if(numslt.gt.1) uvspec(i)=numtype
       endif
    enddo
    !
    allocate( uvmax(numslv),uvsoft(numslv),ercrd(large,0:numslv) )
    !
    peread=0
    ermax=0
    do pti=0,numslv
       open(unit = paramfile_io, file = ene_confname, action = "read", iostat = param_err)
       if (param_err == 0) then
          read(paramfile_io, nml = ene_param) ! FIXME: is this necessary?
          read(paramfile_io, nml = hist)
          close(paramfile_io)
       else
          stop "parameter file does not exist"
       end if
       if(peread.eq.1) then   ! read coordinate parameters from separate file
          open(unit=ecdio,file=ecdfile,status='old')
          read(ecdio,*)        ! comment line
          do i=1,large
             read(ecdio,*,end=3109) q
             if(q.eq.pti) then
                backspace(ecdio)
                if(pti.eq.0) then
                   read(ecdio,*) k,(ecpmrd(k),k=1,8)
                   pecore=0       ! no core region for solute self-energy
                endif
                if(pti.gt.0) then
                   read(ecdio,*) k,(ecpmrd(k),k=1,9),pecore
                   ecdmax=ecpmrd(9)
                endif
                eclbin=ecpmrd(1) ; ecfbin=ecpmrd(2) ; ec0bin=ecpmrd(3)
                finfac=ecpmrd(4) ; ecdmin=ecpmrd(5) ; ecfmns=ecpmrd(6)
                ecdcen=ecpmrd(7) ; eccore=ecpmrd(8)
                if(eccore.lt.tiny) pecore=0
                if(pecore.eq.1) call halt_with_error('ecd')
                exit
             endif
          end do
3109      continue
          close(ecdio)
       end if
       ectmvl=finfac*ecfbin ; ecdmin=ecdmin-ectmvl
       ecfmns=ecfmns-ectmvl ; ecmns0=ecdcen-ectmvl
       ecpls0=2.0e0*ecdcen-ecmns0 ; ecfpls=2.0e0*ecdcen-ecfmns
       eccore=eccore+ecfpls-ecdcen

       rgcnt(1)=(ecfmns-ecdmin)/eclbin
       rgcnt(2)=(ecmns0-ecfmns)/ecfbin
       rgcnt(3)=(ecpls0-ecmns0)/ec0bin
       rgcnt(4)=(ecfpls-ecpls0)/ecfbin
       rgcnt(5)=(eccore-ecfpls)/eclbin
       pesoft=0
       do regn=1,rglmax
          factor=rgcnt(regn)
          if(int(factor).lt.1) call halt_with_error('ecd')
          pesoft=pesoft+nint(factor)
       end do

       pemax=pesoft+pecore
       if(pemax.gt.large) call halt_with_error('siz')

       cdrgvl(0)=ecdmin
       cdrgvl(1)=ecfmns
       cdrgvl(2)=ecmns0
       cdrgvl(3)=ecpls0
       cdrgvl(4)=ecfpls
       cdrgvl(5)=eccore
       cdrgvl(6)=ecdmax
       do regn=1,rglmax-1
          if((regn.eq.1).or.(regn.eq.5)) factor=eclbin
          if((regn.eq.2).or.(regn.eq.4)) factor=ecfbin
          if(regn.eq.3) factor=ec0bin
          iduv=nint((cdrgvl(regn)-cdrgvl(regn-1))/factor)
          if(regn.eq.1) uprgcd(regn)=iduv
          if(regn.gt.1) uprgcd(regn)=uprgcd(regn-1)+iduv
       end do
       uprgcd(rglmax)=pesoft
       uprgcd(rglmax+1)=pemax-1

       do regn=0,rglmax+1
          if(regn.eq.0) iduv=0
          if(regn.gt.0) iduv=uprgcd(regn)
          ercrd(iduv+1,pti)=cdrgvl(regn)
       end do

       if(pecore.eq.0) i=rglmax       ! no core region
       if(pecore.gt.0) i=rglmax+1     ! explicit treatment of core region
       curmx=0
       do regn=1,i
          prvmx=curmx
          curmx=uprgcd(regn)
          if((regn.eq.1).or.(regn.eq.5)) factor=eclbin
          if((regn.eq.2).or.(regn.eq.4)) factor=ecfbin
          if(regn.eq.3) factor=ec0bin
          if(regn.eq.(rglmax+1)) then ! this condition is satisfied only if pecore == 0
             incre=log(ercrd(pemax,pti)/ercrd(pesoft+1,pti))
             factor=incre/real(pecore-1)
          endif
          do iduv=prvmx+1,curmx
             incre=factor*real(iduv-prvmx-1)
             if(regn.le.rglmax) ercrd(iduv,pti)=ercrd(prvmx+1,pti)+incre
             if(regn.eq.(rglmax+1)) then
                ercrd(iduv,pti)=ercrd(prvmx+1,pti)*exp(incre)
             end if
          end do
       end do

       if(pti.eq.0) esmax=pesoft      ! solute self-energy
       if(pti.gt.0) then              ! solute-solvent interaction energy
          uvmax(pti)=pemax
          uvsoft(pti)=pesoft
          ermax=ermax+pemax
       endif
    end do
    !
    allocate( uvcrd(ermax),edens(ermax) )
    if(corrcal.eq.1) allocate( ecorr(ermax,ermax) )
    allocate( escrd(esmax),eself(esmax) )
    if(slttype.eq.1) allocate( aveuv(engdiv,numslv),slnuv(numslv) )
    allocate( avediv(engdiv,2) )
    allocate( minuv(0:numslv),maxuv(0:numslv) )
    !
    i=0
    do pti=1,numslv
       pemax=uvmax(pti)
       do iduv=1,pemax
          i=i+1
          uvcrd(i)=ercrd(iduv,pti)
       end do
    end do
    do iduv=1,esmax
       escrd(iduv)=ercrd(iduv,0)
    end do
    deallocate( ercrd )

    do pti=0,numslv
       minuv(pti)=infty
       maxuv(pti)=-infty
    end do
    voffset = -infty

    call engclear

    ! Output for energy fluctuation
    if((slttype == CAL_SOLN).and.(myrank.eq.0)) then
       open(unit=io_flcuv,file='flcuv.tt',status='new')   ! open flcuv file
    endif

    return
  end subroutine enginit

  subroutine engclear
    use engmain, only: corrcal,slttype,ermax,numslv,esmax,&
         edens,ecorr,eself,slnuv,avslf,engnorm,engsmpl, &
         CAL_SOLN
    implicit none
    integer iduv,iduvp,pti
    edens(1:ermax)=0.0e0
    if(corrcal.eq.1) then
       ecorr(1:ermax,1:ermax)=0.0e0
    endif

    eself(1:esmax)=0.0e0

    if(slttype == CAL_SOLN) then
       slnuv(1:numslv)=0.0e0
    endif
    avslf=0.0e0
    engnorm=0.0e0
    engsmpl=0.0e0
    return
  end subroutine engclear

  subroutine engproc_cleanup
    use engmain, only: slttype, CAL_SOLN, io_flcuv
    use mpiproc
    implicit none
    if((slttype == CAL_SOLN).and.(myrank.eq.0)) then
       endfile(io_flcuv)
       close(io_flcuv)
    endif
  end subroutine engproc_cleanup
  !
  !
  subroutine engconst(stnum, nactiveproc)
    use engmain, only: nummol,maxcnf,skpcnf,corrcal,slttype,wgtslf,&
         estype,sluvid,temp,volume,plmode,&
         maxins,ermax,numslv,esmax,uvspec,&
         edens,ecorr,eself,&
         slnuv,avslf,minuv,maxuv,numslt,sltlist,&
         engnorm,engsmpl,voffset,&
         boxshp, cltype, cell, &
         io_flcuv, &
         SYS_NONPERIODIC, SYS_PERIODIC, &
         EL_COULOMB, EL_PME, &
         CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX, &
         ES_NVT, ES_NPT
    use ptinsrt, only: instslt
    use realcal, only: realcal_proc
    use reciprocal, only: recpcal_init, &
         recpcal_prepare_solute, recpcal_prepare_solvent, recpcal_energy, recpcal_spline_greenfunc, &
         recpcal_self_energy
    use mpiproc                                                      ! MPI
    implicit none
    integer, intent(in) :: stnum, nactiveproc
    integer i,pti,iduv,iduvp,k,q
    integer :: irank
    real engnmfc,pairep,wgtslcf,factor
    integer, dimension(:), allocatable :: insdst,engdst,tplst
    real, dimension(:),    allocatable :: uvengy,svfl
    logical, allocatable :: flceng_stored_g(:,:)
    real, allocatable :: flceng_g(:,:,:)
    real, save :: prevcl(3, 3)
    real, parameter :: tiny = 1.0e-20
    logical :: skipcond
    logical, save :: voffset_initialized = .false.
    logical, save :: pme_initialized = .false.
    call mpi_info                                                    ! MPI
    call mpi_init_active_group(nactiveproc)

    call sanity_check_sluvid()

    ! for soln: maxdst is number of solutes (multiple solute)
    ! for refs: maxdst is number of insertions
    select case(slttype)
    case(CAL_SOLN)
       maxdst=numslt
    case(CAL_REFS_RIGID, CAL_REFS_FLEX)
       maxdst=maxins
    end select

    allocate( tplst(nummol) )
    slvmax=0
    do i=1,nummol
       if(sluvid(i) == 0) then ! 0 == exists in trajectory
          slvmax=slvmax+1
          tplst(slvmax)=i ! which particle is treated by this node?
       end if
    end do
    allocate( tagpt(slvmax) )
    allocate( uvengy(0:slvmax) )
    do k=1,slvmax
       tagpt(k)=tplst(k) ! and copied from tplst
    end do
    deallocate( tplst )

    allocate(flceng_stored(maxdst))
    allocate(flceng(numslv, maxdst))
    flceng_stored(:) = .false.
    flceng(:, :) = 0

    if(myrank < nactiveproc) then
       ! Initialize reciprocal space - grid and charges
       call get_inverted_cell
       if(cltype == EL_PME) then
          if(.not. pme_initialized) then
             call recpcal_init(slvmax,tagpt)
          endif
          
          ! check whether cell size changes
          ! recpcal is called only when cell size differ
          if((.not. pme_initialized) .or. &
               (any(prevcl(:,:) /= cell(:,:)))) call recpcal_spline_greenfunc()
          prevcl(:, :) = cell(:, :)
          
          pme_initialized = .true.
          
          do k=1,slvmax
             i=tagpt(k)
             call recpcal_prepare_solvent(i)
          end do
       endif

       ! cntdst is the loop to select solute MOLECULE from multiple solutes (soln)
       ! cntdst is the iteration no. of insertion (refs)
       do cntdst=1,maxdst
          call get_uv_energy(stnum, wgtslcf, uvengy(0:slvmax), skipcond)
          if(skipcond) cycle
          
          call update_histogram(stnum, wgtslcf, uvengy(0:slvmax))
       end do

       ! for soln only: need to output flceng
       if(slttype == CAL_SOLN) then
          
          allocate(flceng_g(numslv, maxdst, nactiveproc))
          allocate(flceng_stored_g(maxdst, nactiveproc))
          
#ifndef noMPI
          ! gather flceng values to rank 0
          call mpi_gather(flceng_stored, maxdst, mpi_logical, &
               flceng_stored_g, maxdst, mpi_logical, &
               0, mpi_comm_activeprocs, ierror)
          call mpi_gather(flceng, numslv * maxdst, mpi_double_precision, &
               flceng_g, numslv * maxdst, mpi_double_precision, &
               0, mpi_comm_activeprocs, ierror)
#endif
          
          if(myrank == 0) then
             do irank = 1, nactiveproc
                do i = 1, maxdst
                   if(flceng_stored_g(maxdst, irank)) then
                      if(maxdst.eq.1) then
                         write(io_flcuv, 911) (stnum + irank - 1) * skpcnf, (flceng_g(pti, i, irank), pti=1,numslv)
                      else
                         write(io_flcuv, 912) cntdst, (stnum + irank - 1) * skpcnf, (flceng_g(pti, i, irank), pti=1,numslv)
                      endif
911                   format(i9,999f15.5)
912                   format(2i9,999f15.5)
                   endif
                enddo
             enddo
          endif
          deallocate(flceng_g, flceng_stored_g)
       endif
    endif

    deallocate( tagpt,uvengy )
    deallocate(flceng, flceng_stored)

    call mpi_finish_active_group()
    return
  end subroutine engconst
  !
  !
  subroutine engstore(stnum)
    !
    use engmain, only: nummol,maxcnf,engdiv,corrcal,slttype,wgtslf,&
         plmode,ermax,numslv,esmax,temp,&
         edens,ecorr,eself,&
         aveuv,slnuv,avediv,avslf,minuv,maxuv,&
         engnorm,engsmpl,voffset, &
         CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX
    use mpiproc                                                      ! MPI
    implicit none
    integer stnum,i,pti,j,iduv,iduvp,k,q,cntdst, division
    character*10, parameter :: numbers='0123456789'
    character*9 engfile
    character*3 suffeng
    real :: voffset_local, voffset_scale
    real :: factor
    real, parameter :: tiny=1.0e-30
    real, dimension(:), allocatable :: sve1,sve2
    call mpi_info                                                    ! MPI
    !

    ! synchronize voffset
    if(wgtslf == 1) then
       voffset_local = voffset
#ifndef noMPI
       call mpi_allreduce(voffset_local, voffset, 1,&
            mpi_double_precision, mpi_max, mpi_comm_world, ierror)   ! MPI
       ! scale histograms accoording to the maximum voffset
       select case(slttype)
       case (CAL_SOLN)
          voffset_scale = exp((voffset_local - voffset)/temp)
       case (CAL_REFS_RIGID, CAL_REFS_FLEX)
          voffset_scale = exp(-(voffset_local - voffset)/temp)
       end select

       engnorm = engnorm * voffset_scale
       eself(:) = eself(:) * voffset_scale
       if(slttype == CAL_SOLN) slnuv(:) = slnuv(:) * voffset_scale
       edens(:) = edens(:) * voffset_scale
       if(corrcal == 1) ecorr(:, :) = ecorr(:, :) * voffset_scale
#endif
    endif

    ! Gather all information to Master node
#ifndef noMPI
    if(plmode == 2) then                                              ! MPI
       call mpi_reduce(avslf,factor,1,&
            mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
       avslf=factor                                                   ! MPI
       call mpi_reduce(engnorm,factor,1,&
            mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
       engnorm=factor                                                 ! MPI
       call mpi_reduce(engsmpl,factor,1,&
            mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
       engsmpl=factor                                                 ! MPI
       allocate( sve1(esmax) )                                        ! MPI
       do iduv=1,esmax                                                ! MPI
          sve1(iduv)=eself(iduv)                                      ! MPI
       end do
       call mpi_reduce(sve1,eself,esmax,&
            mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
       deallocate( sve1 )                                             ! MPI
       if(slttype == CAL_SOLN) call mympi_reduce_real(slnuv, numslv, mpi_sum, 0)
    endif                                                             ! MPI
    allocate( sve1(0:numslv),sve2(0:numslv) )                         ! MPI
    do pti=0,numslv                                                   ! MPI
       sve1(pti)=minuv(pti)                                           ! MPI
       sve2(pti)=maxuv(pti)                                           ! MPI
    end do
    call mpi_reduce(sve1,minuv,numslv+1,&                             ! MPI
         mpi_double_precision,mpi_min,0,mpi_comm_world,ierror)        ! MPI
    call mpi_reduce(sve2,maxuv,numslv+1,&                             ! MPI
         mpi_double_precision,mpi_max,0,mpi_comm_world,ierror)        ! MPI
    deallocate( sve1,sve2 )                                           ! MPI
    allocate( sve1(ermax) )                                           ! MPI
    do iduv=1,ermax                                                   ! MPI
       sve1(iduv)=edens(iduv)                                         ! MPI
    end do
    call mpi_reduce(sve1,edens,ermax,&                                ! MPI
         mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)        ! MPI
    deallocate( sve1 )                                                ! MPI
#endif
    do iduv=1,ermax
       edens(iduv)=edens(iduv)/engnorm
    end do
    if(corrcal.eq.1) then
       do iduv=1,ermax
#ifndef noMPI
          allocate( sve1(ermax),sve2(ermax) )                          ! MPI
          do iduvp=1,ermax                                        ! MPI
             sve1(iduvp)=ecorr(iduvp,iduv)                              ! MPI
          end do
          call mpi_reduce(sve1,sve2,ermax,&                             ! MPI
               mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)   ! MPI
          do iduvp=1,ermax                                        ! MPI
             ecorr(iduvp,iduv)=sve2(iduvp)                              ! MPI
          end do
          deallocate( sve1,sve2 )                                      ! MPI
#endif
          do iduvp=1,ermax
             ecorr(iduvp,iduv)=ecorr(iduvp,iduv)/engnorm
          end do
       end do
    endif
    do iduv=1,esmax
       eself(iduv)=eself(iduv)/engnorm
    end do
    avslf=avslf/engnorm
    !
    if(myrank.ne.0) go to 7999                                       ! MPI
    division = stnum / (maxcnf / engdiv)
    if(slttype.eq.1) then
       do pti=1,numslv
          aveuv(division,pti)=slnuv(pti)/engnorm
       end do
    endif
    avediv(division,1)=engnorm/engsmpl
    if(slttype.eq.1) avediv(division,2)=voffset-temp*log(avslf)
    if(slttype.ge.2) avediv(division,2)=voffset+temp*log(avslf)
    !
    if(division == engdiv) then
       if(slttype.eq.1) then
          open(unit=75,file='aveuv.tt',status='new')
          do k=1,engdiv
             write(75,751) k,(aveuv(k,pti), pti=1,numslv)
          end do
          endfile(75)
          close(75)
751       format(i5,99999f15.5)
       endif
       if(slttype.eq.1) open(unit=73,file='weight_soln',status='new')
       if(slttype.ge.2) open(unit=73,file='weight_refs',status='new')
       do k=1,engdiv
          if(wgtslf.eq.0) write(73,731) k,avediv(k,1)
          if(wgtslf.eq.1) write(73,732) k,avediv(k,1),avediv(k,2)
       end do
       endfile(73)
       close(73)
731    format(i5,f15.3)
732    format(i5,f15.3,g17.5)
       open(77,file='uvrange.tt',status='new')
       write(77,771)
       do pti=0,numslv
          factor=maxuv(pti)
          if(factor.lt.1.0e5) write(77,772) pti,minuv(pti),factor
          if(factor.ge.1.0e5) write(77,773) pti,minuv(pti),factor
       end do
       endfile(77)
       close(77)
771    format(' species     minimum        maximum')
772    format(i5,2f15.5)
773    format(i5,f15.5,g18.5)
    endif
    !
    j=division/10
    k=division-10*j
    if(engdiv.eq.1) suffeng='.tt'
    if(engdiv.gt.1) suffeng='.'//numbers(j+1:j+1)//numbers(k+1:k+1)
    do cntdst=1,3
       if(cntdst.eq.1) then
          if(slttype.eq.1) engfile='engsln'//suffeng
          if(slttype.ge.2) engfile='engref'//suffeng
       endif
       if((cntdst.eq.2).and.(corrcal.ne.1)) goto 7199
       if((cntdst.eq.2).and.(corrcal.eq.1)) then
          if(slttype.eq.1) engfile='corsln'//suffeng
          if(slttype.ge.2) engfile='corref'//suffeng
       endif
       if(cntdst.eq.3) engfile='slfeng'//suffeng
       open(unit=71,file=engfile,status='new')
       if(cntdst.eq.1) then
          do iduv=1,ermax
             call repval(iduv,factor,pti,'intn')
             write(71,'(g15.7,i5,g25.15)') factor,pti,edens(iduv)
          enddo
       endif
       if(cntdst.eq.2) then
          do iduv=1,ermax
             do iduvp=1,ermax
                factor=ecorr(iduvp,iduv)
                if(factor.gt.tiny) write(71,'(g25.15)') factor
                if(factor.le.tiny) write(71,'(i1)') 0
             end do
          end do
       endif
       if(cntdst.eq.3) then
          do iduv=1,esmax
             call repval(iduv,factor,pti,'self')
             write(71,'(g15.7,g25.15)') factor,eself(iduv)
          end do
       endif
       endfile(71)
       close(71)
7199   continue
    end do
7999 continue                                                         ! MPI
    !
    return
  end subroutine engstore

  ! Calculate interaction energy between solute and solvent
  subroutine get_uv_energy(stnum, weighting, uvengy, has_error)
    use engmain, only: nummol,numatm,maxcnf,skpcnf,corrcal,slttype,wgtslf,&
         estype,sluvid,temp,volume,plmode,&
         maxins,ermax,numslv,esmax,uvspec,&
         edens,ecorr,eself,&
         slnuv,avslf,minuv,maxuv,numslt,sltlist,&
         engnorm,engsmpl,voffset,&
         boxshp, cltype, cell, &
         SYS_NONPERIODIC, SYS_PERIODIC, &
         EL_COULOMB, EL_PME, &
         CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX, &
         ES_NVT, ES_NPT, sitepos ! FIXME: remove
    use ptinsrt, only: instslt
    use realcal, only: realcal_proc, realcal_prepare, realcal_cleanup, realcal_bare, realcal_self
    use reciprocal, only: recpcal_init, &
         recpcal_prepare_solute, recpcal_prepare_solvent, recpcal_energy, recpcal_spline_greenfunc, &
         recpcal_self_energy
    use mpiproc                                                      ! MPI

    implicit none
    integer, intent(in) :: stnum
    real, intent(inout) :: uvengy(0:slvmax), weighting
    logical, intent(out) :: has_error

    integer :: i, k, q
    integer(8) :: current_solute_hash
    integer(8) :: solute_hash = 0
    real :: pairep, residual, factor
    real, save :: usreal
    logical, save :: initialized = .false.

    has_error = .false.
    ! determine / pick solute structure
    select case(slttype) 
    case(CAL_SOLN)
       tagslt=sltlist(cntdst)
       call check_mol_configuration(has_error)
       if(has_error) return
    case(CAL_REFS_RIGID, CAL_REFS_FLEX)
       tagslt=sltlist(1)
       if(.not. initialized) then
          call instslt(weighting,'init')
       endif
       call instslt(weighting,'proc')
       if((stnum.eq.maxcnf).and.(cntdst.eq.maxdst)) then
          call instslt(weighting,'last')
       endif
       initialized = .true.
    end select

    ! At this moment all coordinate in the system is determined
    call realcal_prepare

    uvengy(:) = 0
    ! Calculate system-wide values
    if(cltype == EL_PME) then
       call recpcal_prepare_solute(tagslt)
       call realcal_proc(tagslt, tagpt, slvmax, uvengy)
       call recpcal_self_energy(uvengy(0))
    endif

    ! solute-solute self energy
    pairep = 0.0
    residual = 0.0
    current_solute_hash = get_solute_hash() ! FIXME: if this tuns into a bottleneck, add conditionals
    if(current_solute_hash == solute_hash .or. &
         (slttype == CAL_REFS_RIGID .and. solute_hash /= 0)) then 
       ! For refs part, the structure of solute may change because of rotation & translation upon insertion, 
       ! though the self energy will not change.
       pairep = usreal ! reuse
    else
       call realcal_self(tagslt, pairep) ! calculate self-interaction
       usreal = pairep
    endif
    solute_hash = current_solute_hash
    call residual_ene(tagslt, tagslt, residual)
    uvengy(0) = uvengy(0) + pairep + residual

    ! solute-solvent pair
    do k=1,slvmax
       i=tagpt(k)
       if(i.eq.tagslt) cycle

       pairep = 0
       factor = 0
       if(cltype == EL_PME) then ! called only when PME, non-self interaction
          call residual_ene(tagslt, i, pairep)
          call recpcal_energy(tagslt, i, factor)
          pairep = pairep + factor
       else
          call realcal_bare(tagslt,i,pairep) ! Bare coulomb solute-solvent interaction
       endif
       uvengy(k) = uvengy(k) + pairep
    enddo

    call realcal_cleanup
  end subroutine get_uv_energy

  subroutine update_histogram(stnum, stat_weight, uvengy)
    use engmain, only: wgtslf, plmode, estype, slttype, corrcal, volume, temp, uvspec, &
         ermax, numslv, &
         slnuv, avslf,&
         minuv, maxuv, &
         edens, ecorr, eself, &
         engnorm, engsmpl, voffset, &
         io_flcuv, &
         CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX,&
         ES_NVT, ES_NPT
    use mpiproc
    implicit none
    integer, intent(in) :: stnum
    real, intent(in) :: uvengy(0:slvmax), stat_weight

    integer, allocatable :: insdst(:), engdst(:)
    real, allocatable :: svfl(:)

    integer :: i, k, q, iduv, iduvp, pti
    real :: factor, engnmfc, pairep

    logical, save :: voffset_initialized = .false.

    allocate(insdst(ermax), engdst(ermax))

    if(wgtslf.eq.0) engnmfc=1.0e0
    if(wgtslf.eq.1) then
       factor=uvengy(0)
       if(.not. voffset_initialized) then
          voffset=factor
          voffset_initialized = .true.
       endif
       factor=factor-voffset               ! shifted by offset
       select case(slttype)
       case (CAL_SOLN)
          engnmfc=exp(factor/temp)
       case (CAL_REFS_RIGID, CAL_REFS_FLEX)
          engnmfc=exp(-factor/temp)
       end select
    endif
    if(estype == ES_NPT) call volcorrect(engnmfc)
    if(slttype == CAL_REFS_RIGID .or. slttype == CAL_REFS_FLEX) engnmfc = engnmfc*stat_weight
    !
    engnorm=engnorm+engnmfc               ! normalization factor
    engsmpl=engsmpl+1.0e0                 ! number of sampling

    select case(estype)
    case (ES_NVT)
       avslf=avslf+1.0e0
    case (ES_NPT)
       avslf=avslf+volume
    end select


    ! self energy histogram
    call getiduv(0, uvengy(0), iduv)
    eself(iduv) = eself(iduv) + engnmfc
    minuv(0) = min(minuv(0), uvengy(0))
    maxuv(0) = max(maxuv(0), uvengy(0))

    insdst(:) = 0
    flceng(:, cntdst) = 0.0e0                     ! sum of solute-solvent energy
    flceng_stored(cntdst) = .true.
    ! interaction energy histogram
    do k = 1, slvmax
       i=tagpt(k)
       if(i.eq.tagslt) cycle
       pti=uvspec(i)
       if(pti.le.0) call halt_with_error('eng')

       pairep=uvengy(k)
       call getiduv(pti,pairep,iduv)

       insdst(iduv)=insdst(iduv)+1
       flceng(pti, cntdst) = flceng(pti, cntdst) + pairep    ! sum of solute-solvent energy

       minuv(pti) = min(minuv(pti), pairep)
       maxuv(pti) = max(maxuv(pti), pairep)
    end do

    if(slttype == CAL_SOLN) then
       slnuv(:) = slnuv(:) + flceng(:, cntdst) * engnmfc
    endif

    do iduv=1,ermax
       k=insdst(iduv)
       if(k.gt.0) edens(iduv)=edens(iduv)+engnmfc*real(k)
    enddo
    if(corrcal.eq.1) then
       do iduv=1,ermax
          k=insdst(iduv)
          if(k == 0) cycle

          do iduvp=1,ermax
             q=insdst(iduvp)
             if(q == 0) cycle

             ecorr(iduvp,iduv) = ecorr(iduvp,iduv) + engnmfc*real(k)*real(q)
          end do
       end do
    endif

    deallocate(insdst, engdst)
  end subroutine update_histogram

  !
  subroutine residual_ene(i, j, pairep)
    use engmain, only: screen, volume, numsite, mol_charge, cltype
    implicit none
    integer, intent(in) :: i, j
    real, intent(inout) :: pairep
    real :: rtp1, rtp2, epcl
    integer :: is, js, ismax, jsmax, ati, atj
    real, parameter :: pi = 3.141592653589793283462
    if(cltype == 0) return

    rtp1 = mol_charge(i)
    rtp2 = mol_charge(j)
    epcl=pi*rtp1*rtp2/screen/screen/volume
    if(i.eq.j) epcl=epcl/2.0e0 ! self-interaction
    pairep=pairep-epcl
  end subroutine residual_ene
  !
  subroutine volcorrect(engnmfc)
    use engmain, only:  nummol,maxsite,numatm,temp,numsite,sluvid,&
         cltype,screen,charge,specatm,volume
    implicit none
    integer i,ati,sid,stmax
    real pi,factor,engnmfc
    engnmfc=volume*engnmfc
    if(cltype.ne.0) then                                 ! Ewald and PME
       pi=real(4)*atan(real(1))
       factor=0.0e0
       do i=1,nummol
          if(sluvid(i).le.1) then                          ! physical particle
             stmax=numsite(i)
             do sid=1,stmax
                ati=specatm(sid,i)
                factor=factor+charge(ati)
             end do
          endif
       end do
       factor=pi*factor*factor/screen/screen/volume/2.0e0
       engnmfc=engnmfc*exp(factor/temp)
    endif
    return
  end subroutine volcorrect
  !
  subroutine update_cell_info()
    use engmain, only: cell, celllen
    implicit none
    integer :: i
    do i = 1, 3
       celllen(i) = sqrt(sum(cell(:, i) ** 2))
    enddo
  end subroutine update_cell_info

  subroutine get_inverted_cell
    use engmain, only:  cell,invcl,volume
    implicit none
    integer m,k
    volume=cell(1,1)*cell(2,2)*cell(3,3)&
         +cell(1,2)*cell(2,3)*cell(3,1)+cell(1,3)*cell(2,1)*cell(3,2)&
         -cell(1,3)*cell(2,2)*cell(3,1)&
         -cell(1,2)*cell(2,1)*cell(3,3)-cell(1,1)*cell(2,3)*cell(3,2)
    invcl(1,1)=cell(2,2)*cell(3,3)-cell(2,3)*cell(3,2)
    invcl(1,2)=cell(1,3)*cell(3,2)-cell(1,2)*cell(3,3)
    invcl(1,3)=cell(1,2)*cell(2,3)-cell(1,3)*cell(2,2)
    invcl(2,1)=cell(2,3)*cell(3,1)-cell(2,1)*cell(3,3)
    invcl(2,2)=cell(1,1)*cell(3,3)-cell(1,3)*cell(3,1)
    invcl(2,3)=cell(1,3)*cell(2,1)-cell(1,1)*cell(2,3)
    invcl(3,1)=cell(2,1)*cell(3,2)-cell(2,2)*cell(3,1)
    invcl(3,2)=cell(1,2)*cell(3,1)-cell(1,1)*cell(3,2)
    invcl(3,3)=cell(1,1)*cell(2,2)-cell(1,2)*cell(2,1)
    
    do m=1,3
       do k=1,3
          invcl(k,m)=invcl(k,m)/volume
       end do
    end do
    call update_cell_info
  end subroutine get_inverted_cell


  ! binsearch returns the smallest index (ret) which satisfies
  ! coord(ret) >= v
  ! Special cases are as follows:
  ! v < coord(1)  ==>  ret = 0 (out of range)
  ! v > coord(n)  ==>  ret = n (infinite)
  subroutine binsearch(coord, n, v, ret)
    implicit none 
    real, intent(in) :: coord(n)
    integer, intent(out) :: ret
    real, intent(in) :: v
    integer, intent(in) :: n
    integer :: rmin, rmax, rmid
    if(v < coord(1)) then
       ret = 0
       return
    endif
    if(v > coord(n)) then
       ret = n
       return
    endif

    rmin = 1
    rmax = n + 1
    do
       if(rmax - rmin <= 1) then
          exit
       endif
       rmid = (rmin + rmax - 1) / 2
       if(v > coord(rmid)) then
          rmin = rmid + 1
       else
          rmax = rmid + 1
       endif
    enddo
    ret = rmin - 1
  end subroutine binsearch

  ! returns the position of the bin corresponding to energy coordinate value
  subroutine getiduv(pti,engcoord,iduv)
    use engmain, only: ermax,numslv,uvmax,uvcrd,esmax,escrd,stdout
    use mpiproc, only: halt_with_error, warning
    implicit none
    integer, intent(in) :: pti
    real, intent(in) :: engcoord
    integer, intent(out) :: iduv
    integer :: k,idpick,idmax,picktest
    real :: egcrd
    real, parameter :: warn_threshold = 1e+8
    if(pti.eq.0) idmax=esmax               ! solute self-energy
    if(pti.gt.0) idmax=uvmax(pti)          ! solute-solvent interaction
    idpick=0
    if(pti.gt.1) then
       do k=1,pti-1
          idpick=idpick+uvmax(k)
       end do
    endif
    iduv=idpick
    if(pti == 0) then
       call binsearch(escrd, idmax, engcoord, picktest)
    else
       call binsearch(uvcrd(idpick+1), idmax, engcoord, picktest)
    endif
    iduv = picktest + idpick

    if(picktest == idmax .and. engcoord < warn_threshold) then
       ! Feature #52: put a warning if energy exceeds max binning region and pecore = 0
       ! Since it is hard to distinguish pecore = 0 bin at this moment,
       ! it is determined by looking engcoord: if too low, it is sprious
       write(stdout, '(A,g12.4,A,i3,A)'), '  energy of ', engcoord, ' for ', pti, '-th species'
       call warning('mbin')
    endif

    ! FIXME: clean up the following
    if(iduv.le.idpick) then
       iduv=idpick+1                                 ! smallest energy mesh
       write(stdout,199) engcoord,pti
199    format('  energy of ',g12.4,' for ',i3,'-th species')
       call halt_with_error('min')
    endif
    if(iduv.gt.(idpick+idmax)) iduv=idpick+idmax    ! largest energy mesh

    return
  end subroutine getiduv

  ! Check system consistency: either test particle or solvent must exist
  subroutine sanity_check_sluvid()
    use engmain, only: slttype, nummol, sluvid, CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX
    use mpiproc, only: halt_with_error
    implicit none

    ! sanity check
    if(any(sluvid(:) < 0) .or. any(sluvid(:) > 3)) call halt_with_error('bug')
    select case(slttype)
    case(CAL_SOLN)
       ! sluvid should be 0 (solvent) or 1 (solute)
       if(any(sluvid(:) >= 2)) call halt_with_error('par')
    case(CAL_REFS_RIGID, CAL_REFS_FLEX)
       ! sluvid should be 0 (solvent), 2, 3 (test particles)
       if(any(sluvid(:) == 1)) call halt_with_error('par')
    end select

    ! solute / test particle must exist
    if(all(sluvid(:) == 0)) call halt_with_error('par')
    ! solvent must exist
    if(all(sluvid(:) /= 0)) call halt_with_error('par')
  end subroutine sanity_check_sluvid

  ! Check whether molecule is within specified region (of sltcnd)
  subroutine check_mol_configuration(is_invalid)
    use ptinsrt, only: sltpstn, get_molecule_com
    implicit none
    logical, intent(out) :: is_invalid
    integer :: sltstat
    real :: com(3)
    is_invalid = .false.
    
    call get_molecule_com(tagslt, com)
    call check_mol_configuration_impl(com, is_invalid)
  end subroutine check_mol_configuration

  subroutine check_mol_configuration_impl(com, is_invalid)
    use engmain, only: inscnd, lwreg, upreg, boxshp, SYS_NONPERIODIC, invcl, celllen
    use ptinsrt, only: get_system_com
    use mpiproc, only: halt_with_error
    implicit none
    real, intent(in) :: com(3)
    logical, intent(out) :: is_invalid
    real :: system_com(3), dx(3)
    real :: distance

    is_invalid = .false.

    select case(inscnd)
    case(0)
       return
    case(1) ! sphere geometry
       call get_system_com(system_com)
       dx(:) = com(:) - system_com(:)
       distance = sqrt(dot_product(dx, dx))
    case(2) ! slab (only z-axis is constrained) configuration
       if(boxshp == SYS_NONPERIODIC) call halt_with_error('slb')
       call get_system_com(system_com)
       dx(:) = com(:) - system_com(:)
       distance = abs(dot_product(invcl(3,:), dx(:))) * celllen(3)
    end select
    
    if(distance > lwreg .and. distance < upreg) return
    is_invalid = .true.
  end subroutine check_mol_configuration_impl

  ! get the hashed function of solute coordinate
  integer(8) function get_solute_hash()
    use utility, only: hash
    use engmain, only: sitepos, mol_begin_index, numsite
    implicit none

    get_solute_hash = hash(sitepos(1:3, mol_begin_index(tagslt):(mol_begin_index(tagslt+1) - 1)), numsite(tagslt) * 3)
  end function get_solute_hash

  
  subroutine repval(iduv,factor,pti,caltype)
    use engmain, only: ermax,numslv,uvmax,uvsoft,uvcrd,esmax,escrd
    use mpiproc, only: halt_with_error
    implicit none
    integer iduv,idpt,pti,cnt,idpick,idmax,idsoft
    real factor
    character*4 caltype
    if(caltype.eq.'self') then
       if(iduv.lt.esmax) factor=(escrd(iduv)+escrd(iduv+1))/2.0e0
       if(iduv.eq.esmax) factor=escrd(esmax)
    endif
    if(caltype.eq.'intn') then
       idpick=0
       do cnt=1,numslv
          idpick=idpick+uvmax(cnt)
          if(iduv.le.idpick) exit
       enddo
       pti=cnt
       idpick=idpick-uvmax(pti)
       idsoft=uvsoft(pti)
       idmax=uvmax(pti)
       idpt=iduv-idpick
       if((idpt.lt.0).or.(idpt.gt.idmax)) call halt_with_error('ecd')
       if(idpt.le.idsoft) factor=(uvcrd(iduv)+uvcrd(iduv+1))/2.0e0
       if((idpt.gt.idsoft).and.(idpt.lt.idmax)) then
          factor=sqrt(uvcrd(iduv)*uvcrd(iduv+1))
       endif
       if(idpt.eq.idmax) factor=uvcrd(iduv)
    endif
    return
  end subroutine repval
end module engproc
