
module engproc
  implicit none
contains
  !
  !  procedure for constructing energy distribution functions
  !
  subroutine enginit
    !
    use engmain, only: numtype,nummol,engdiv,corrcal,slttype,&
         moltype,sluvid,&
         ermax,numslv,uvmax,uvsoft,esmax,uvspec,&
         uvcrd,edens,ecorr,escrd,eself,&
         aveuv,slnuv,avediv,minuv,maxuv,numslt,sltlist,&
         ene_param, ene_confname
    !
    ! FIXME
    implicit real(a-h,k-z)
    implicit integer(i,j)

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
          call eng_stop('typ')
       endif
    end do
    iduv=0
    if(numslt.le.0) iduv=9
    if((slttype.ge.2).and.(numslt.ne.1)) iduv=9
    if(iduv.ne.0) call eng_stop('num')
    allocate( sltlist(numslt) )
    do i=1,numslt
       sltlist(i)=tplst(i)
    end do
    if((slttype.ge.2).and.(sltlist(1).ne.nummol)) call eng_stop('ins')
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
          if(pti.eq.sltmltp) call eng_stop('typ')
          if(pti.gt.sltmltp) uvspec(i)=pti-1
       endif
       if(sluvid(i).ne.0) then            ! solute
          if(pti.ne.sltmltp) call eng_stop('typ')
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
#       include "param_eng"
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
                if(pecore.eq.1) call eng_stop('ecd')
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
          if(int(factor).lt.1) call eng_stop('ecd')
          pesoft=pesoft+nint(factor)
       end do

       pemax=pesoft+pecore
       if(pemax.gt.large) call eng_stop('siz')

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
          if(regn.eq.(rglmax+1)) then
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

    call engclear
    !
    return
  end subroutine enginit
  !
  !
  subroutine engclear
    use engmain, only: corrcal,slttype,ermax,numslv,esmax,&
         edens,ecorr,eself,slnuv,avslf,engnorm,engsmpl
    integer iduv,iduvp,pti
    do iduv=1,ermax
       edens(iduv)=0.0e0
    end do
    if(corrcal.eq.1) then
       do iduv=1,ermax
          do iduvp=1,ermax
             ecorr(iduvp,iduv)=0.0e0
          end do
       end do
    endif
    do iduv=1,esmax
       eself(iduv)=0.0e0
    end do
    if(slttype.eq.1) then
       do pti=1,numslv
          slnuv(pti)=0.0e0
       end do
    endif
    avslf=0.0e0
    engnorm=0.0e0
    engsmpl=0.0e0
    return
  end subroutine engclear
  !
  !
  subroutine engconst(stnum)
    !
    use engmain, only: nummol,maxcnf,skpcnf,corrcal,slttype,wgtslf,&
         estype,sluvid,temp,volume,plmode,&
         maxins,ermax,numslv,esmax,uvspec,&
         edens,ecorr,eself,&
         slnuv,avslf,minuv,maxuv,numslt,sltlist,&
         engnorm,engsmpl,voffset,&
         boxshp, cltype, cell
    use ptinsrt, only: instslt
    use realcal_blk, only: realcal_proc
    use reciprocal, only: recpcal
    use mpiproc                                                      ! MPI
    integer, parameter :: flcio=91                    ! IO unit for flcuv
    integer stnum,cntdst,maxdst,tagslt,slvmax,i,pti,iduv,iduvp,k,q
    integer ptinit,ptskip,dsinit,dsskip
    real engnmfc,pairep,wgtslcf,factor
    integer, dimension(:), allocatable :: insdst,engdst,tagpt,tplst
    real, dimension(:),    allocatable :: uvengy,flceng,svfl
    real, save :: prevcl(3, 3), usreal
    real, parameter :: tiny = 1.0e-20
    call mpi_info                                                    ! MPI
    !
    if((slttype.eq.1).and.(myrank.eq.0).and.(stnum.eq.skpcnf)) then
       open(unit=flcio,file='flcuv.tt',status='new')   ! open flcuv file
    endif
    !
    call sltcnd(q,0,'sys')
    if((slttype.eq.1).and.(q.ne.1)) q=9
    if((slttype.ge.2).and.(q.ne.2)) q=9
    if(q.eq.9) call eng_stop('par')
    if(slttype.eq.1) maxdst=numslt
    if(slttype.ge.2) maxdst=maxins
    !
    if(plmode.eq.0) then
       ptinit=myrank ; ptskip=nprocs
       dsinit=0 ; dsskip=1
    endif
    if(plmode.eq.1) then
       ptinit=0 ; ptskip=1
       dsinit=myrank ; dsskip=nprocs
    endif
    !
    allocate( tplst(nummol) )
    slvmax=0
    do i=1+ptinit,nummol,ptskip
       if((slttype.eq.1) .or. &
            ((slttype.ge.2).and.(sluvid(i).eq.0))) then
          slvmax=slvmax+1
          tplst(slvmax)=i ! which particle is treated by this node?
       end if
    end do
    allocate( insdst(ermax),engdst(ermax),tagpt(slvmax) )
    allocate( uvengy(0:slvmax),flceng(numslv),svfl(numslv) )
    do k=1,slvmax
       tagpt(k)=tplst(k) ! and copied from tplst
    end do
    deallocate( tplst )
    !
    call cellinfo
    !
    if(stnum.eq.skpcnf) then                ! Ewald and PME initialization
       call recpcal(0,0,factor,slvmax,tagpt,'alloct')
    endif

    ! check whether cell size changes
    q=1
    if(stnum.eq.skpcnf) q=0
    if(stnum.gt.skpcnf) then
       factor=0.0e0
       do i=1,3
          do k=1,3
             pairep=abs(prevcl(k,i)-cell(k,i))
             if(factor.lt.pairep) factor=pairep
          end do
       end do
       if(factor.gt.tiny) q=0
    endif
    do i=1,3
       do k=1,3
          prevcl(k,i)=cell(k,i)
       end do
    end do
    ! recpcal is called only when cell size differ
    if(q.eq.0) call recpcal(0,0,factor,slvmax,tagpt,'preeng')
    !
    do k=1,slvmax
       i=tagpt(k)
       call recpcal(i,i,factor,slvmax,tagpt,'slvenv')
    end do
    !
    do cntdst=1,maxdst
       if(slttype.eq.1) then                 ! solution system
          tagslt=sltlist(cntdst)
          call sltcnd(q,tagslt,'pos')
          if(q.eq.0) go to 99999
          if((q.ne.0).and.(q.ne.1)) call eng_stop('slt')
       endif
       if(slttype.ge.2) then                 ! solute insertion
          tagslt=sltlist(1)
          if((stnum.eq.skpcnf).and.(cntdst.eq.1)) then
             call instslt(wgtslcf,'init')
          endif
          call instslt(wgtslcf,'proc')
          if((stnum.eq.maxcnf).and.(cntdst.eq.maxdst)) then
             call instslt(wgtslcf,'last')
          endif
          if(slttype.eq.2) then               ! rigid solute
             if((stnum.eq.skpcnf).and.(cntdst.eq.1)) then   ! initialization
                call realcal(tagslt,tagslt,usreal)           ! self real part
             endif
          endif
          if(mod(cntdst-1,dsskip).ne.dsinit) go to 99999
       endif
       !
       call recpcal(tagslt,tagslt,factor,slvmax,tagpt,'sltsys')

       uvengy(:) = 0
       if(cltype /= 0) then ! called only when ewald-type real part
          if(boxshp == 0) stop "Ewald / PME is selected, but box is not periodic!"
          call realcal_proc(tagslt, tagpt, slvmax, uvengy)
       endif

       do k=0,slvmax
          if(k.eq.0) i=tagslt   ! solute self
          if(k.gt.0) then       ! solute-solvent pair
             i=tagpt(k)
             if(i.eq.tagslt) cycle
          endif

          if(cltype /= 0 .and. k /= 0) then ! called only when ewald-type, pair interaction
             pairep = 0
             call residual_ene(tagslt, i, pairep)
          else
             if((slttype.eq.2).and.(k.eq.0)) then ! rigid solute self real part
                pairep=usreal
             else
                call realcal(tagslt,i,pairep) ! usual case or self-interaction
             endif
          endif
          call recpcal(tagslt,i,factor,slvmax,tagpt,'energy')
          pairep = pairep + factor
          uvengy(k) = uvengy(k) + pairep
       enddo
       !
       if(wgtslf.eq.0) engnmfc=1.0e0
       if(wgtslf.eq.1) then
          factor=uvengy(0)
          if((stnum.eq.skpcnf).and.(cntdst.eq.(1+dsinit))) then
             if(dsinit.eq.0) voffset=factor
#ifndef noMPI
             if(plmode.eq.1) call mpi_bcast(voffset,1,&
                  mpi_double_precision,0,mpi_comm_world,ierror)    ! MPI
#endif
          endif
          factor=factor-voffset               ! shifted by offset
          if(slttype.eq.1) engnmfc=exp(factor/temp)
          if(slttype.ge.2) engnmfc=exp(-factor/temp)
       endif
       if(estype.eq.2) call volcorrect(engnmfc)
       if(slttype.ge.2) engnmfc=engnmfc*wgtslcf
       !
       engnorm=engnorm+engnmfc               ! normalization factor
       engsmpl=engsmpl+1.0e0                 ! number of sampling
       if(estype.le.1) avslf=avslf+1.0e0
       if(estype.eq.2) avslf=avslf+volume
       !
       do iduv=1,ermax
          insdst(iduv)=0
       enddo
       do pti=1,numslv
          flceng(pti)=0.0e0                   ! sum of solute-solvent energy
       enddo
       do k=0,slvmax
          if(k.eq.0) pti=0                    ! solute self
          if(k.gt.0) then                     ! solute-solvent pair
             i=tagpt(k)
             if(i.eq.tagslt) cycle
             pti=uvspec(i)
             if(pti.le.0) call eng_stop('eng')
          endif
          pairep=uvengy(k)
          call getiduv(pti,pairep,iduv)
          if(k.eq.0) eself(iduv)=eself(iduv)+engnmfc
          if(k.gt.0) then
             insdst(iduv)=insdst(iduv)+1
             flceng(pti)=flceng(pti)+pairep    ! sum of solute-solvent energy
          endif
          if(pairep.lt.minuv(pti)) minuv(pti)=pairep        ! minimum energy
          if(pairep.gt.maxuv(pti)) maxuv(pti)=pairep        ! maximum energy
       end do
       !
#ifndef noMPI
       if(plmode.eq.0) then
          call mpi_allreduce(insdst,engdst,ermax,&
               mpi_integer,mpi_sum,mpi_comm_world,ierror)       ! MPI
          call mpi_allreduce(flceng,svfl,numslv,mpi_double_precision,&
               mpi_sum,mpi_comm_world,ierror)                   ! MPI
          do iduv=1,ermax                                         ! MPI
             insdst(iduv)=engdst(iduv)                                  ! MPI
          end do
          do pti=1,numslv                                         ! MPI
             flceng(pti)=svfl(pti)                                      ! MPI
          end do
       endif
#endif
       if(slttype.eq.1) then
          do pti=1,numslv
             slnuv(pti)=slnuv(pti)+flceng(pti)*engnmfc
          end do
          if(myrank.eq.0) then
             if(maxdst.eq.1) then
                write(flcio, 911) stnum,(flceng(pti), pti=1,numslv)
             endif
             if(maxdst.gt.1) then
                write(flcio, 912) cntdst,stnum, (flceng(pti), pti=1,numslv)
             endif
          endif
911       format(i9,999999999f15.5)
912       format(2i9,999999999f15.5)
       endif
       !
       do iduv=1+ptinit,ermax,ptskip
          k=insdst(iduv)
          if(k.gt.0) edens(iduv)=edens(iduv)+engnmfc*real(k)
       enddo
       if(corrcal.eq.1) then
          do iduv=1+ptinit,ermax,ptskip
             k=insdst(iduv)
             if(k.gt.0) then
                do iduvp=1,ermax
                   q=insdst(iduvp)
                   if(q.gt.0) then
                      factor=engnmfc*real(k)*real(q)
                      ecorr(iduvp,iduv)=ecorr(iduvp,iduv)+factor
                   endif
                end do
             endif
          end do
       endif
       !
99999  continue
90000  continue
    end do
    !
    if((slttype.eq.1).and.(myrank.eq.0).and.(stnum.eq.maxcnf)) then
       endfile(flcio) ; close(flcio)                   ! close flcuv file
    endif
    !
    deallocate( insdst,engdst,tagpt,uvengy,flceng,svfl )
    !
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
         engnorm,engsmpl,voffset
    use mpiproc                                                      ! MPI
    integer stnum,i,pti,j,iduv,iduvp,k,q,cntdst
    character*10, parameter :: numbers='0123456789'
    character*9 engfile
    character*3 suffeng
    real factor
    real, parameter :: tiny=1.0e-30
    real, dimension(:), allocatable :: sve1,sve2
    call mpi_info                                                    ! MPI
    !
#ifndef noMPI
    if(plmode.eq.1) then                                             ! MPI
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
       do iduv=1,esmax                                           ! MPI
          sve1(iduv)=eself(iduv)                                       ! MPI
       end do
       call mpi_reduce(sve1,eself,esmax,&
            mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
       deallocate( sve1 )                                             ! MPI
    endif                                                            ! MPI
    allocate( sve1(0:numslv),sve2(0:numslv) )                        ! MPI
    do pti=0,numslv                                             ! MPI
       sve1(pti)=minuv(pti)                                           ! MPI
       sve2(pti)=maxuv(pti)                                           ! MPI
    end do
    call mpi_reduce(sve1,minuv,numslv+1,&                             ! MPI
         mpi_double_precision,mpi_min,0,mpi_comm_world,ierror)       ! MPI
    call mpi_reduce(sve2,maxuv,numslv+1,&                             ! MPI
         mpi_double_precision,mpi_max,0,mpi_comm_world,ierror)       ! MPI
    deallocate( sve1,sve2 )                                          ! MPI
    allocate( sve1(ermax) )                                          ! MPI
    do iduv=1,ermax                                             ! MPI
       sve1(iduv)=edens(iduv)                                         ! MPI
    end do
    call mpi_reduce(sve1,edens,ermax,&                                ! MPI
         mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)       ! MPI
    deallocate( sve1 )                                               ! MPI
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
    i=stnum/(maxcnf/engdiv)
    if(slttype.eq.1) then
       do pti=1,numslv
          aveuv(i,pti)=slnuv(pti)/engnorm
       end do
    endif
    avediv(i,1)=engnorm/engsmpl
    if(slttype.eq.1) avediv(i,2)=voffset-temp*log(avslf)
    if(slttype.ge.2) avediv(i,2)=voffset+temp*log(avslf)
    !
    if(stnum.eq.maxcnf) then
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
    j=i/10
    k=i-10*j
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
7011   continue
    end do
7999 continue                                                         ! MPI
    !
    return
  end subroutine engstore
  !
  !
  subroutine eng_stop(type)
    use engmain, only: io6
    use mpiproc                                                      ! MPI
    character*3 type
    if(type.eq.'typ') write(io6,991)
    if(type.eq.'num') write(io6,992)
    if(type.eq.'ins') write(io6,993)
    if(type.eq.'par') write(io6,994)
    if(type.eq.'slt') write(io6,995)
    if(type.eq.'crd') write(io6,996)
    if(type.eq.'eng') write(io6,997)
    if(type.eq.'siz') write(io6,998)
    if(type.eq.'min') write(io6,999)
    if(type.eq.'ecd') write(io6,981)
    if(type.eq.'fst') write(io6,982)
991 format(' The number of solute types is incorrectly set')
992 format(' The number of solute molecules is incorrectly set')
993 format(' The solute numbering is incorrect for insertion')
994 format(' The input parameter is incorrectly set')
995 format(' The input parameter is incorrect for solute')
996 format(' The coordinate system is incorrectly set')
997 format(' Inconsistency is present in the program')
998 format(' The number of energy-coordinate meshes is too large')
999 format(' The minimum of the energy coordinate is too large')
981 format(' The energy-coordinate system is inconsistent')
982 format(' The first particle needs to be the solute')
    call mpi_abend()                                                ! MPI
    stop
  end subroutine eng_stop
  !
  !
  subroutine realcal(i,j,pairep)
    !
    use engmain, only:  nummol,maxsite,numatm,boxshp,numsite,&
         elecut,lwljcut,upljcut,cmbrule,cltype,screen,&
         charge,ljene,ljlen,specatm,sitepos,&
         cell,invcl,volume
    integer i,j,is,js,ismax,jsmax,ati,atj,m,k
    real pi,reelcut,pairep,ljeps,ljsgm,chr2,rst,dis2,rtp1,rtp2
    real eplj,epcl,xst(3),clm(3),swth
    real, parameter :: infty=1.0e50      ! essentially equal to infinity
    !
    pi=real(4)*atan(real(1))
    if(i.eq.j) reelcut=infty             ! no cutoff for self-interaction
    if(i.ne.j) then
       if(cltype.eq.0) then                  ! bare Coulomb
          if(boxshp.eq.0) reelcut=infty
          if(boxshp.ne.0) reelcut=elecut
       endif
       if(cltype.ne.0) reelcut=elecut        ! Ewald and PME
    endif
    !
    pairep=0.0e0
    ismax=numsite(i)
    jsmax=numsite(j)
    !
    do is=1,ismax
       do js=1,jsmax
          if((i.eq.j).and.(is.gt.js)) go to 1599
          ati=specatm(is,i)
          atj=specatm(js,j)
          do m=1,3
             xst(m)=sitepos(m,ati)-sitepos(m,atj)
          end do
          if(boxshp.ne.0) then              ! when the system is periodic
             do k=1,3
                rst=0.0e0
                do m=1,3
                   rst=rst+invcl(k,m)*xst(m)
                end do
                clm(k)=real(nint(rst))
             end do
             do m=1,3
                rst=0.0e0
                do k=1,3
                   rst=rst+cell(m,k)*clm(k)
                end do
                xst(m)=xst(m)-rst             ! get the nearest distance between i,j
             end do
          endif
          dis2=xst(1)*xst(1)+xst(2)*xst(2)+xst(3)*xst(3)
          rst=sqrt(dis2)
          if((i.eq.j).or.(rst.gt.upljcut)) eplj=0.0e0
          if((i.ne.j).and.(rst.le.upljcut)) then
             ljeps=sqrt(ljene(ati)*ljene(atj))
             if(cmbrule.eq.0) ljsgm=(ljlen(ati)+ljlen(atj))/2.0e0
             if(cmbrule.eq.1) ljsgm=sqrt(ljlen(ati)*ljlen(atj))
             rtp1=ljsgm*ljsgm/dis2
             rtp2=rtp1*rtp1*rtp1
             eplj=4.0e0*ljeps*rtp2*(rtp2-1.0e0)
             if(rst.gt.lwljcut) then    ! CHARMM form of switching function
                rtp1=lwljcut*lwljcut
                rtp2=upljcut*upljcut
                swth=(2.0e0*dis2+rtp2-3.0e0*rtp1)*(dis2-rtp2)*(dis2-rtp2)&
                     /(rtp2-rtp1)/(rtp2-rtp1)/(rtp2-rtp1)
                eplj=swth*eplj
             endif
          endif
          if(rst.ge.reelcut) epcl=0.0e0
          if(rst.lt.reelcut) then
             chr2=charge(ati)*charge(atj)
             if((i.eq.j).and.(is.eq.js)) then
                if(cltype.eq.0) rtp1=real(0)                  ! bare Coulomb
                if(cltype.ne.0) rtp1=screen                   ! Ewald and PME
                epcl=-chr2*rtp1/sqrt(pi)
             endif
             if((i.ne.j).or.(is.ne.js)) then
                if(cltype.eq.0) rtp1=real(0)                  ! bare Coulomb
                if(cltype.ne.0) rtp1=screen*rst               ! Ewald and PME
                if(i.eq.j) epcl=-chr2*derf(rtp1)/rst
                if(i.ne.j) epcl=chr2*(1.0e0-derf(rtp1))/rst
             endif
          endif
          pairep=pairep+(eplj+epcl)
1599      continue
1502      continue
       end do
1501   continue
    end do
    !
    if(cltype.ne.0) then                                 ! Ewald and PME
       call residual_ene(i, j, pairep)
    endif
    !
    return
  end subroutine realcal
  !
  subroutine residual_ene(i, j, pairep)
    use engmain, only: screen, volume, specatm, numsite, charge, cltype
    implicit none
    integer, intent(in) :: i, j
    real, intent(inout) :: pairep
    real :: rtp1, rtp2, epcl
    integer :: is, js, ismax, jsmax, ati, atj
    real, parameter :: pi = 3.141592653589793283462
    if(cltype == 0) stop "Error: residual_ene: called when cltype == 0, cannot happen"
    ismax = numsite(i)
    jsmax = numsite(j)
    rtp1=0.0e0
    do is=1,ismax
       ati=specatm(is,i)
       rtp1=rtp1+charge(ati)
    enddo
    rtp2=0.0e0
    do js=1,jsmax
       atj=specatm(js,j)
       rtp2=rtp2+charge(atj)
    enddo
    epcl=pi*rtp1*rtp2/screen/screen/volume
    if(i.eq.j) epcl=epcl/2.0e0 ! self-interaction
    pairep=pairep-epcl
  end subroutine residual_ene
  !
  subroutine volcorrect(engnmfc)
    use engmain, only:  nummol,maxsite,numatm,temp,numsite,sluvid,&
         cltype,screen,charge,specatm,volume
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
  !
  subroutine cellinfo
    use engmain, only:  cell,invcl,volume
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
    return
  end subroutine cellinfo
  !
  !
  !
  subroutine binsearch(coord, n, value, ret)
    real, intent(in) :: coord(n)
    integer, intent(out) :: ret
    real, intent(in) :: value
    integer, intent(in) :: n
    integer :: rmin, rmax, rmid
    if(value < coord(1)) then
       ret = 0
       return
    endif
    if(value > coord(n)) then
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
       if(value > coord(rmid)) then
          rmin = rmid + 1
       else
          rmax = rmid + 1
       endif
    enddo
    ret = rmin - 1
  end subroutine binsearch
  !
  subroutine getiduv(pti,factor,iduv)
    use engmain, only: ermax,numslv,uvmax,uvcrd,esmax,escrd,io6
    integer pti,iduv,k,idpick,idmax,picktest
    real factor,egcrd
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
       call binsearch(escrd, idmax, factor, picktest)
    else
       call binsearch(uvcrd(idpick+1), idmax, factor, picktest)
    endif
    iduv = picktest + idpick

    if(iduv.le.idpick) then
       iduv=idpick+1                                 ! smallest energy mesh
       write(io6,199) factor,pti
199    format('  energy of ',g12.4,' for ',i3,'-th species')
       call eng_stop('min')
    endif
    if(iduv.gt.(idpick+idmax)) iduv=idpick+idmax    ! largest energy mesh

    return
  end subroutine getiduv

  subroutine sltcnd(systype,tagslt,type)
    use engmain, only: nummol,sluvid
    use ptinsrt, only: sltpstn
    integer systype,tagslt,i,uvi,cntuv(2)
    real xst(3)
    character*3 type
    if(type.eq.'sys') then
       systype=9
       do uvi=1,2
          cntuv(uvi)=0
       enddo
       do i=1,nummol
          uvi=sluvid(i)
          if(uvi.eq.3) uvi=2
          if(uvi.gt.0) cntuv(uvi)=cntuv(uvi)+1
       end do
       if((cntuv(1).ne.0).and.(cntuv(2).eq.0)) systype=1
       if((cntuv(1).eq.0).and.(cntuv(2).ne.0)) systype=2
    endif
    if(type.eq.'pos') call sltpstn(systype,xst,'solutn',tagslt)
    return
  end subroutine sltcnd

  subroutine repval(iduv,factor,pti,caltype)
    use engmain, only: ermax,numslv,uvmax,uvsoft,uvcrd,esmax,escrd
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
       if((idpt.lt.0).or.(idpt.gt.idmax)) call eng_stop('ecd')
       if(idpt.le.idsoft) factor=(uvcrd(iduv)+uvcrd(iduv+1))/2.0e0
       if((idpt.gt.idsoft).and.(idpt.lt.idmax)) then
          factor=sqrt(uvcrd(iduv)*uvcrd(iduv+1))
       endif
       if(idpt.eq.idmax) factor=uvcrd(iduv)
    endif
    return
  end subroutine repval
end module engproc
