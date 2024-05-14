!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! RelEPro: a program to calculate relativistic one-electron properties. Written by
!
! Wenli Zou,  Email: qcband@gmail.com
! Institute of Modern Physics, Northwest University, Xi'an, China
!
! Geometry, basis functions, and MOs are read from a MOLDEN or FCHK file.
! Only the standardized MOLDEN file by Molden2AIM or MultiWFN (since Version 3.8) is fully supported. In addition, the basis
!   functions for each atom should be sorted in ascending order according to angular quantum numbers. Thus the SP functions have
!   to be edited manually.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program RelEPro
  use Constant, only : maxlq, Zero, One
  implicit real(kind=8) (a-h,o-z)
  parameter (NOption=10, NJobTyp=20, iinp=5, iout=6, idat=41, igto=42)
  dimension     :: IOP(NOption), JobTyp(NJobTyp)
  character*200 :: ctmp, tag
  character*10  :: symm
  logical       :: ifham, ifpro, ifecp=.false., ifdbg=.false.
  dimension     :: dumb(1)
  allocatable   :: iza(:), za(:), xyz(:,:), rms(:), mapatm(:), mxla(:), nshlla(:,:,:), ialbs(:,:), ialbc(:,:), iptexp(:,:),  &
                   iptcon(:,:), gtoexp(:), gtocon(:), c2scoef(:)
  allocatable   :: pchg(:), pxyz(:,:), qdat(:,:)
  allocatable   :: den(:), pdens(:), pdenv(:), pdenw(:), smat(:), tmat(:), vmat(:), wmat(:), xmat(:), rmat(:), occ(:), ene(:),  &
                   ispn(:), symm(:), cmo(:), cmor(:), cmou(:), scr(:,:)

  call prt_title(iout)

  call fdate(ctmp)
  write(iout,"(/,1x,a)") trim(ctmp)

!---------------------------------------------------------------------------------------------------------------------------------
!  1. read input and open data file
!---------------------------------------------------------------------------------------------------------------------------------
  call open_dfile(iinp,iout,ifham,ifpro,ifdbg,irelc,NOption,IOP,NJobTyp,JobTyp,idat,ctmp)

!---------------------------------------------------------------------------------------------------------------------------------
! 2. Read data-1: natom & nmo
!---------------------------------------------------------------------------------------------------------------------------------
  if(IOP(1) == 1) then
    call rd_molden_1(iout,idat,ctmp,natom,nmo)
  else if(IOP(1) == 2) then
    call rd_fchk_1(iout,idat,ctmp,natom,nmo)
  end if
  ! reset IOP(5) = loca
  if(natom == 1 .and. IOP(5) /= 1) then
    IOP(5) = 1
    write(iout,"(/,' ### Warning: LOCA has been reset to 1.',/)")
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 3. Read data-2: maxtyp, ncbas, nsbas, npcar, npsph, nshell, ncoeff
!---------------------------------------------------------------------------------------------------------------------------------
  open(igto,file='gtoscr.dat')
  rewind(igto)
  if(IOP(1) == 1) then
    call rd_molden_2(iout,idat,igto,ctmp,tag, natom, maxtyp, ispher, ncbas, nsbas, npcar, npsph, nshell, ncoeff)
  else if(IOP(1) == 2) then
    call rd_fchk_2(iout,idat,igto,irelc,ctmp, natom, nmo, maxtyp, ispher, ncbas, nsbas, npcar, npsph, nshell, ncoeff)
  end if

  write(iout,"(' Number of atoms                      :',i8)") natom
  write(iout,"(' Max L                                :',i8)") maxtyp
  write(iout,"(' Number of primitive shells           :',i8)") nshell
  write(iout,"(' Number of contracted basis functions :',i8,' (car.)',i8,' (sph.)')") ncbas, nsbas
  write(iout,"(' Number of primitive basis functions  :',i8,' (car.)',i8,' (sph.)')") npcar, npsph
  if(ispher == 0) then
    write(iout,"(' MOs are printed in                   :','     Cartesian functions')")
    write(iout,*)
    write(iout,"(' ### Warning: Cartesian functions may lead to numerical errors.')")
  else
    write(iout,"(' MOs are printed in                   :','     Spherical functions')")
  end if

  ncss = max(npcar*npcar,800)
  ndm2 = irelc*irelc
  ndm3 = irelc*ndm2
  nscr = 4
  if(irelc == 2) nscr = nscr*8  ! see ZHEEV
  npatt=(maxtyp+1)*(maxtyp+2)*(maxtyp+3)*(maxtyp*3+2)
  npatt= npatt/12
  allocate(iza(natom), za(natom), xyz(3,natom),  &
           ! RMS radii; map of equivalent atoms; max_L of each atom; numbers of primitive and contracted shells of each L;
           rms(natom), mapatm(natom), mxla(natom), nshlla(1+maxlq,natom,2),  &
           ! starting position of each L in primitive sperical & Cartesian basis functions;
           ialbs(1+maxlq,natom), ialbc(1+maxlq,natom),  &
           ! starting positions of GTO exponents and contraction coefficients for each L;
           iptexp(1+maxlq,natom), iptcon(1+maxlq,natom),  &
           ! GTO exponents; contraction coefficients; coefficients of Cartesian to spherical functions
           gtoexp(nshell), gtocon(ncoeff), c2scoef(npatt), stat=ierr)
  if(ierr /= 0) then
    write(iout,"(' ### Insufficient Memory (1)!')")
    call estop(1)
  end if
  allocate(den(nsbas*nsbas*ndm3), occ(nmo), ene(nmo), ispn(nmo), symm(nmo), cmo(nsbas*nmo*ndm2), pdens(ncss*ndm3),  &
    smat(ncss), scr(ncss,nscr), stat=ierr)
  if(ierr /= 0) then
    write(iout,"(' ### Insufficient Memory (2)!')")
    call estop(1)
  end if
  if(ifham) then
    allocate(pdenv(ncss*ndm3), pdenw(ncss*ndm3), tmat(ncss), vmat(ncss), wmat(ncss*ndm2), xmat(ncss*ndm3), rmat(ncss*ndm3),  &
      stat=ierr)
    if(ierr /= 0) then
      write(iout,"(' ### Insufficient Memory (3)!')")
      call estop(1)
    end if
  end if
  if(IOP(6) == 1) then
    allocate(cmor(npsph*nmo*ndm2), cmou(npsph*nmo*ndm2), stat=ierr)
    if(ierr /= 0) then
      write(iout,"(' ### Insufficient Memory (4)!')")
      call estop(1)
    end if
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 4. Read data-3: xyz coordinates
!    In the case of ECP, za = iza - ncore
!---------------------------------------------------------------------------------------------------------------------------------
  if(IOP(1) == 1) then
    call rd_molden_3(iout,idat,ctmp,natom,iza,za,xyz)
  else if(IOP(1) == 2) then
    call rd_fchk_3(iout,idat,ctmp,natom,ifecp,iza,za,xyz)
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 5. Read external point charges
!---------------------------------------------------------------------------------------------------------------------------------
  npchg = 0
  call rd_pchar(iinp,iout,npchg,dumb,scr,ctmp)
  if(npchg > 0) then
    allocate(pchg(npchg), pxyz(3,npchg))
    call rd_pchar(iinp,iout,npchg,pchg,pxyz,ctmp)
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 6. Read basis functions from the backup in igto
!---------------------------------------------------------------------------------------------------------------------------------
  call rd_basfun(iout,igto,natom,mxla,nshlla,ialbs,ialbc,iptexp,iptcon,gtoexp,gtocon,scr)

!---------------------------------------------------------------------------------------------------------------------------------
! 7. Read data-4: occupation, MOs (in spherical functions), and density matrix (in spherical functions)
!---------------------------------------------------------------------------------------------------------------------------------
  nbas = ncbas
  if(ispher == 1) nbas = nsbas

  if(IOP(1) == 1) then

    call rd_molden_4(iout,idat,nbas,nmo,occtot,occ,ene,ispn,symm,ctmp,scr)
    ! MO: Cartesian to spherical
    if(ispher == 1) then
      call acopy(nsbas*nmo,scr,cmo)
    else
      call filc2s(1,(IOP(1)-2),maxtyp,c2scoef)
      call MOCar2Sph(natom,mxla,nshlla(1,1,2),ncbas,nsbas,nmo,scr,cmo,c2scoef)
    end if

  else if(IOP(1) == 2) then

    if(ispher == 0) then
      write(iout,"(' ### Error: FCHK in Cartesian functions is not allowed.')")
      call estop(1)
    end if

    call rd_fchk_4(iout,idat,irelc,(ndm2*nbas),nmo,occtot,occ,ene,ispn,symm,ctmp,cmo,scr)
    call rd_fchk_den(iout,idat,IOP(7),irelc,nbas,ctmp,scr)
    if(irelc == 1) then
      ! density matrix: l.t. --> sqr.
      call lt2sqr(0,nbas,scr,den)
    else
      ! density matrix: l.t. --> sqr.
      call lt2sqrC(.true.,2*nbas,scr,den,den)
      ! real and imaginary parts of MO coefficients are saved separately
      call acopy(4*nsbas*nmo,cmo,scr)
      call ReImSplit(nsbas*2,nmo,scr,cmo)
    end if

    !! MO: Cartesian to spherical
    !if(ispher == 0) then
    !  call filc2s(1,(IOP(1)-2),maxtyp,c2scoef)
    !  if(irelc == 1) then
    !    call acopy(ncbas*nmo,cmo,scr)
    !    call MOCar2Sph(natom,mxla,nshlla(1,1,2),ncbas,nsbas,nmo,scr,cmo,c2scoef)
    !  else if(irelc == 2) then
    !    call MOCar2Sph2C(natom,mxla,nshlla(1,1,2),ncbas,nsbas,nmo,cmo,scr,c2scoef)
    !    ! real and imaginary parts of MO coefficients are saved separately
    !    call ReImSplit(nsbas*2,nmo,scr,cmo)
    !  end if
    !end if

  end if

  ! density matrix
  if(IOP(7) == -1) then
    write(iout,"(/,' Computing density matrix in spherical functions ...')")
    if(irelc == 1) then
      call calc_den(nsbas,nmo,cmo,occ,den)
    else if(irelc == 2) then
      call calc_den2c(2*nsbas,nmo,cmo,nint(occtot),den)
    end if
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 8. S, T, V, and W in primitive spherical functions
!    The basis functions are ordered according to -L,...,+L for each atom.
!
!    NOTE. One-center V & W will be calculated for ALR & AU. Since the Hamiltonian is not finally constructed and used in this
!    program, this will not lead to problem.
!---------------------------------------------------------------------------------------------------------------------------------
  if(ifham) then
    write(iout,"(/,' Computing one-electron integrals ...')")
    if(ifecp .and. IOP(2) /= 1) then
      write(iout,"(/,' All-electron relativistic Hamiltonian is not allowed for ECP.')")
      call estop(1)
    end if
    call oneel_drv(iout,(IOP(1)-2),IOP(2),irelc,natom,iza,za,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,  &
      smat,tmat,vmat,wmat,xmat,rmat,scr,ncss,nscr)
  else
    write(iout,"(/,' Computing overlap integrals ...')")
    ! rms is also setup here
    call oneel_drv(iout,(IOP(1)-2),1,1,natom,iza,dumb,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,  &
      smat,dumb,dumb,dumb,pdens,scr(1,1),scr(1,2),ncss,nscr-1)
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 9. X2C.
!
!    NOTE. This is not a complete version of X2C program. Only the X and R matrices will be computed.
!---------------------------------------------------------------------------------------------------------------------------------
  if(ifham) then
    write(iout,"(/,' Constructing Hamiltonian ...')")
    call x2c_drv(iout,irelc,IOP(2),natom,nuqatm,mapatm,iza,mxla,nshlla,iptexp,gtoexp,ialbs,  &
      npsph,smat,tmat,vmat,wmat,xmat,rmat,scr)
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 10. contracted spherical-dimensional density matrix (den) and MOs (cmo) to primitive spherical-dimensional density matrix
!    (pdens) and MOs (cmou).
!---------------------------------------------------------------------------------------------------------------------------------
  write(iout,"(/,' Computing density matrix in the dimension of primitive spherical functions ...')")
  if(irelc == 1) then
    call drv_contra(1,natom,maxlq,nsbas,npsph,mxla,nshlla,iptcon,gtocon,den,pdens,scr(1,1),scr(1,2))
  else if(irelc == 2) then
    call drv_contra2(1,natom,maxlq,nsbas,npsph,mxla,nshlla,iptcon,gtocon,den,pdens,scr(1,1),scr(1,5))
  end if

  ! check density matrix
  call checkden(iout,natom,npsph,irelc,nchag,occtot,za,smat,pdens,scr)

  ! population
  if(ifdbg) call population(iout,natom,npsph,irelc,maxtyp,iza,za,mxla,ialbs,nshlla,smat,pdens,scr)

  ! LOCA = 3: delete atomic-overlap blocks in P
  if(ifpro .and. IOP(5) == 3) call abdens(irelc,natom,npsph,mxla,ialbs,nshlla,pdens,scr)

  if(ifpro .and. IOP(6) == 1) then
    write(iout,"(/,' Computing MOs in the dimension of primitive spherical functions ...')")
    call drv_mocontra(1,irelc,natom,maxlq,nsbas,npsph,nmo,mxla,nshlla,iptcon,gtocon,cmo,cmou,scr)
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 11. X2C transformed density matrix and MOs in the primitive spherical dimension
!---------------------------------------------------------------------------------------------------------------------------------
  if(ifham) then
    write(iout,"(/,' Performing density matrix transformation for properties ...')")

    if(IOP(2) == 1) then
      pdenv=pdens
      if(IOP(6) == 1) cmor=cmou
    else
      call X2CDenTran(irelc,npsph,xmat,rmat,pdens,pdenv,pdenw,scr)

      if(IOP(6) == 1) then
        write(iout,"(/,' Performing MO transformation ...')")
        call X2CMOTran(irelc,npsph,nmo,xmat,rmat,cmou,cmor,cmou)
      end if
    end if
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 12. property calculations
!---------------------------------------------------------------------------------------------------------------------------------
  if(ifpro) then
    write(iout,"(//,23x,'.',23('-'),'.',/,23x,'|  P R O P E R T I E S  |',/,23x,'`',23('-'),'`')")

    ! DIP
    if(JobTyp(1) == 1) then
      write(iout,"(//,1x,'electric dipole moment (DIP)',/,1x,28('='))")
      call dipcalc(iout,(IOP(1)-2),irelc,natom,IOP(2),za,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,pdenv,pdenw,  &
        npchg,pchg,pxyz, scr,ncss*nscr)
    end if

    ! APT
    if(JobTyp(2) == 1) then
      write(iout,"(//,1x,'atomic polar tensor (APT) without CPHF',/,1x,38('='))")
      write(iout,"(/,1x,'* NOTE. They are not translation invariant!')")
      call aptcalc(iout,(IOP(1)-2),irelc,natom,IOP(2),nchag,iza,za,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,  &
        smat,pdenv,pdenw,scr,ncss*nscr)
    end if

    ! DDD
    if(JobTyp(3) == 1) then
      write(iout,"(//,1x,'atomic DDD charges (DDD)',/,1x,24('='))")
      call dddcalc(iout,(IOP(1)-2),irelc,natom,IOP(2),nchag,iza,za,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,  &
        smat,pdenv,pdenw,scr,ncss*nscr)
    end if

    ! ED
    if(JobTyp(11) == 1) then
      write(iout,"(//,1x,'effective contact density (ED)',/,1x,30('='))")
      if(ifecp) then
        write(iout,"(/,' ED is incompatible with ECP.')")
        call estop(1)
      end if
      call filc2s(0,(IOP(1)-2),maxtyp,c2scoef)
      if(IOP(5) /= 1) then
        call edcalclo(iout,(IOP(1)-2),irelc,natom,IOP(2),IOP(3),IOP(5),iza,za,rms,npsph,mxla,nshlla,ialbs,iptexp, &
          gtoexp,c2scoef,pdenv,pdenw, scr,ncss*nscr)
      else
        call edcalc(iout,(IOP(1)-2),irelc,natom,IOP(2),IOP(3),IOP(6),iza,za,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,  &
          gtoexp,c2scoef,pdenv,pdenw, nmo,occ,ene,ispn,symm,cmor,cmou, scr,ncss*nscr)
      end if
    end if

    ! CD (scalar)
    if(JobTyp(12) == 1 .and. irelc == 1) then
      write(iout,"(//,1x,'contact density with PC effects (CD)',/,1x,36('='))")
      if(ifecp) then
        write(iout,"(/,' CD is incompatible with ECP.')")
        call estop(1)
      end if
      call filc2s(0,(IOP(1)-2),maxtyp,c2scoef)
      if(IOP(5) /= 1) then
        !???
        !call cdcalclo(iout,(IOP(1)-2),irelc,natom,IOP(2),IOP(3),IOP(5),iza,za,rms,npsph,mxla,nshlla,ialbs,iptexp,gtoexp,c2scoef,  &
        !  pdenv,pdens,scr,ncss*nscr)
      else
        if(IOP(2) == 1) then
          call cdcalc(iout,(IOP(1)-2),irelc,natom,IOP(2),IOP(3),iza,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,  &
            c2scoef,pdens,dumb,dumb, scr,ncss*nscr)
        else
          call cdcalc(iout,(IOP(1)-2),irelc,natom,IOP(2),IOP(3),iza,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,  &
            c2scoef,pdenv,pdenw,xmat, scr,ncss*nscr)
        end if
      end if
    end if

    ! EFG
    if(JobTyp(13) == 1) then
      write(iout,"(//,1x,'electric field gradient (EFG)',/,1x,29('='))")
      if(ifecp) then
        write(iout,"(/,' EFG is incompatible with ECP.')")
        call estop(1)
      end if
      ! read custom NQM values (in millibarn)
      allocate(qdat(2,natom))
      call rd_nqm(iinp,iout,natom,qdat,ctmp)

      call filc2s(0,(IOP(1)-2),maxtyp,c2scoef)
      if(IOP(5) /= 1) then
        call efgcalclo(iout,(IOP(1)-2),irelc,natom,IOP(2),IOP(3),IOP(5),iza,za,rms,xyz,npsph,mxla,nshlla,ialbs,iptexp, &
          gtoexp,c2scoef,pdenv,pdenw, npchg,pchg,pxyz,qdat, scr,ncss*nscr)
      else
        call efgcalc(iout,(IOP(1)-2),irelc,natom,IOP(2),IOP(3),IOP(6),iza,za,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,  &
          gtoexp,c2scoef,pdenv,pdenw, npchg,pchg,pxyz,qdat, nmo,occ,ene,ispn,symm,cmor,cmou, scr,ncss*nscr)
      end if
      deallocate(qdat)
    end if

  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 99. the last step
!---------------------------------------------------------------------------------------------------------------------------------
  99 continue
  close(idat)
  close(igto,status='delete')
  deallocate(iza, za, xyz, rms, mapatm, mxla, nshlla, ialbs, ialbc, iptexp, iptcon, gtoexp, gtocon, c2scoef)
  deallocate(den, occ, ene, ispn, symm, cmo, pdens, smat, scr)
  if(npchg > 0) deallocate(pchg, pxyz)
  if(ifham) deallocate(pdenv, pdenw, tmat, vmat, wmat, xmat, rmat)
  if(IOP(6) == 1) deallocate(cmor, cmou)
  call fdate(ctmp)
  write(iout,"(/,1x,a)") trim(ctmp)

end

!--- END
