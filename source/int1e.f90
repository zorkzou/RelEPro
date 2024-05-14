!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver of one electron integral calculation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine oneel_drv(iout,mpatt,iham,irelc,natom,iza,za,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,  &
  smat,tmat,vmat,wmat,amat,bmat,scr,ncss,nscr)
  use Constant, only : One, f_4cc
  implicit real(kind=8) (a-h,o-z)
  dimension :: iza(*), za(*), rms(*), xyz(*), mxla(*), nshlla(*), gtoexp(*), c2scoef(*), smat(npsph*npsph), tmat(npsph*npsph), &
               vmat(npsph*npsph), wmat(npsph*npsph,*), amat(*), bmat(*), scr(ncss,nscr)  ! actually nscr=4 is enough for W_so
  dimension :: dumb(1)

  npss = npsph*npsph
  npct = npcar*(npcar+1)/2

  call filc2s(0,mpatt,maxtyp,c2scoef)
  call setrms(natom,iza,rms)

  select case(iham)
    case(1)  ! non-rel. Hamiltonian. T and V are not needed at present.
      ! s
      call pdst_drv(iout,0,0,.false.,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,amat,dumb,dumb,scr,ncss*nscr)
      call lt2sqr(0,npcar,amat,bmat)
      call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,smat,c2scoef,scr(1,1),scr(1,2))

      ! v
      !call pnai_drv(iout,0,0,0,0,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,amat,scr,ncss*nscr)
      !call lt2sqr(0,npcar,amat,bmat)
      !call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,vmat,c2scoef,scr(1,1),scr(1,2))

    case(2,3,5)  ! rel. Hamiltonian with full V & W
      ! s
      call pdst_drv(iout,0,0,.false.,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,amat,dumb,dumb,scr,ncss*nscr)
      call lt2sqr(0,npcar,amat,bmat)
      call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,smat,c2scoef,scr(1,1),scr(1,2))

      ! t
      call pdst_drv(iout,1,0,.false.,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,amat,dumb,dumb,scr,ncss*nscr)
      call lt2sqr(0,npcar,amat,bmat)
      call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,tmat,c2scoef,scr(1,1),scr(1,2))

      ! v
      call pnai_drv(iout,0,0,0,0,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,amat,scr,ncss*nscr)
      call lt2sqr(0,npcar,amat,bmat)
      call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,vmat,c2scoef,scr(1,1),scr(1,2))

      ! w
      call pnai_drv(iout,irelc,0,0,0,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,amat,scr,ncss*nscr)
      do ic = 1, irelc*irelc
        i = npct*(ic-1) + 1
        j = (ic+1)/3
        call lt2sqr(j,npcar,amat(i),bmat)
        call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,wmat(1,ic),c2scoef,scr(1,1),scr(1,2))
        !write(*,"(' ic = ',i1)") ic
        !write(*,"(5d16.6)")wmat(:,ic)
        call ascale(npss,(One/f_4cc),wmat(1,ic),wmat(1,ic))
      end do

    case(4,6)  ! rel. Hamiltonian with one-center V & W
      ! s
      call pdst_drv(iout,0,0,.false.,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,amat,dumb,dumb,scr,ncss*nscr)
      call lt2sqr(0,npcar,amat,bmat)
      call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,smat,c2scoef,scr(1,1),scr(1,2))

      ! t
      call pdst_drv(iout,1,0,.false.,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,amat,dumb,dumb,scr,ncss*nscr)
      call lt2sqr(0,npcar,amat,bmat)
      call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,tmat,c2scoef,scr(1,1),scr(1,2))

      ! v
      call pnai_drv(iout,0,0,1,0,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,amat,scr,ncss*nscr)
      call lt2sqr(0,npcar,amat,bmat)
      call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,vmat,c2scoef,scr(1,1),scr(1,2))

      ! w
      call pnai_drv(iout,irelc,0,1,0,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,amat,scr,ncss*nscr)
      do ic = 1, irelc*irelc
        i = npct*(ic-1) + 1
        j = (ic+1)/3
        call lt2sqr(j,npcar,amat(i),bmat)
        call car2sph(0,natom,mxla,nshlla,npcar,npsph,bmat,wmat(1,ic),c2scoef,scr(1,1),scr(1,2))
        !write(*,"(' ic = ',i1)") ic
        !write(*,"(10d14.6)")wmat(:,ic)
        call ascale(npss,(One/f_4cc),wmat(1,ic),wmat(1,ic))
      end do

    case default
      write(iout,"(/,' Unknown Hamiltonian.')")
      call estop(1)

  end select

  return
end subroutine oneel_drv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine calculates dipole moment integrals in primitive Cartesian Gaussian functions.
! The basic method was described in
! H. Taketa, S. Huzinaga, and K. O-ohata, J. Phys. Soc. Jap. 21(11), 2313, 1966.
!
! Input:
!   LSD            : 0 for S (for test only), 1 for dV/dF, 2 for dW_sf/dF, 3 for dW_sf/dF, dW_x/dF, dW_y/dF, dW_z/dF
!   NDer           : 0/1, the order of derivative
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   natom          : number of atoms
!   xyzmol         : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   nbasis         : number of primitive Cartesian functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   NTT            : nbasis*(nbasis+1)/2
!
! Output:
!   DMat(NTT,NMat) : L.T. matrices in primitive Cartesian basis functions
!
! Scratch:
!   V(Mem)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine pdip_drv(iout,LSD,NDer,mpatt,natom,xyzmol,maxtyp,nshell,nbasis,mxla,nshlla,gtoexp,NTT,DMat,V,Mem)
  use Constant, only : maxlq, Zero
  implicit real(kind=8) (a-h,o-z)
  dimension :: xyzmol(*),mxla(natom),nshlla(0:maxlq,natom),gtoexp(*), DMat(NTT,*), V(Mem)
  allocatable :: LPat(:),K(:)
  logical :: IfSO

  ! check
  if(NDer < 0 .or. NDer > 1) then
    write(iout,"(/,' Error in sub. pdip_drv: NDer has been tested only for 0 and 1.')")
    call estop(1)
  end if
  if(LSD == 0) then
    NMat = (NDer+1)*(NDer+2)/2
    LSDD = 0
    IfSO = .false.
  else if(LSD == 1) then
    NMat = 3*(NDer+1)*(NDer+2)/2
    LSDD = 1
    IfSO = .false.
  else if(LSD == 2) then
    if(NDer == 1) then
      write(iout,"(/,' Error in sub. pdip_drv: dW/dF gradient has not been implemented yet.')")
      call estop(1)
    end if
    NMat = 3
    LSDD = 2
    IfSO = .false.
  else if(LSD == 3) then
    if(NDer == 1) then
      write(iout,"(/,' Error in sub. pdip_drv: dW/dF gradient has not been implemented yet.')")
      call estop(1)
    end if
    NMat = 12
    LSDD = 2
    IfSO = .true.
  else
    write(iout,"(/,' Error in sub. pdip_drv: LSD must be 0, 1, 2, or 3.')")
    call estop(1)
  end if

  ! ==================
  ! | initialization |
  ! ==================
  MaxI1 = MaxTyp + (LSDD/2)
  MaxI2 = MaxI1  + (LSDD/2)
  MaxSum=(MaxTyp+1)*(MaxTyp+2)*(MaxTyp+3)/6
  MaxXYZ=(MaxTyp+1)*(MaxTyp+2)/2
  MaxXYZ2=MaxXYZ*MaxXYZ

  !NSS=NBasis*NBasis
  LGex= 11
  if(LSDD < 2) then
    Ltmp= 4*3
    if(NDer==1) Ltmp= 8*3
  else
    Ltmp= 28*3
  end if
  Lprf=  6

  INrm = 1
  IE34 = INrm + MaxSum
  IRa  = IE34 + NSHELL
  IRb  = IRa  + 3
  Iab  = IRb  + 3
  Iexp = Iab  + 3
  Ifdt = Iexp + 1
  I2If = Ifdt + MaxI2+1
  IIoG = I2If + MaxI1+1
  IClv = IIoG + MaxI1+1
  Ipab = IClv + (MaxI2+1)*(MaxI2+1)
  Iprf = Ipab + (MaxI2+1)*(MaxI2+1)*6
  IGex = Iprf + Lprf
  Iddc = IGex + LGex
  if(NDer==0)then
    Itmp = Iddc
  else if(NDer==1)then
    Itmp = Iddc + 3
  else if(NDer==2)then
    Itmp = Iddc + 6
  end if
  ISlmn= Itmp + Ltmp
  ISshl= ISlmn+ MaxXYZ2*NMat
  LenV = ISshl+ MaxSum*MaxSum*NMat - 1
  if(LenV > Mem) then
    write(iout,"(/,' Error in sub. pdip_drv: LenV > Mem.')")
    call estop(1)
  end if

  KOff = 1
  KMap = KOff + NSHELL
  KLqn = KMap + NSHELL
  LenK = KLqn + NSHELL

  allocate(LPat(MaxSum*3),K(LenK))
  V(1:LenV) = Zero
  K = 0

  ! N!
  call FacCal(MaxI2,V(Ifdt))

  ! Generate (x,y,z) patterns in MOLDEN or Gaussian mode
  call GnPatt(MaxTyp,mpatt,LPat)

  ! common normalization factors
  call NormFac(MaxTyp,LPat,V(INrm))

  ! primitive exponents EXX and common factors EXXP=EXX^0.75
  call ExpP34(NSHELL,gtoexp,V(IE34))

  ! setup indices of primitive shells
  call SetShl(natom,nshell,mxla,nshlla,K(KOff),K(KMap),K(KLqn))

  DMat(:,1:NMat) = Zero

  call DDIPmain(MaxTyp,MaxXYZ,MaxSum,NSHELL,NTT,NDer,NMat,LSDD,IfSO,xyzmol,K(KOff),K(KMap),K(KLqn),LPat,V(Iprf),V(INrm),     &
    V(Ifdt),gtoexp,V(IE34),V(Iexp),V(IRa),V(IRb),V(Iab),V(Ipab),V(I2If),V(IIoG),V(IClv),V(IGex),V(Iddc),V(Itmp),V(ISlmn),    &
    V(ISshl),DMat)

  deallocate(LPat,K)

  Return
End subroutine pdip_drv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main subroutine to calculate the S, dV/dF and dW/dF matrices in primitive Cartesian Gaussian functions.
!
! Input:
!   MaxTyp    The highest angular quantum number present in the basis.
!   MaxXYZ    (MaxTyp+1)*(MaxTyp+2)/2
!   MaxSum    (MaxTyp+1)*(MaxTyp+2)*(MaxTyp+3)/6
!   NSHELL    Number of primitive shells
!   NTT       NTT=NBAS*(NBAS+1)/2. NBAS is the number of primitive Cartesian Gaussian functions.
!   NDer      0/1/2, the order of derivative. NDer=2 has not been tested
!   NMat      the number of (derivative) matrices
!   LSDD      0 for S, 1 for dV/dF, 2 for dW/dF
!   IfSO      whether compute dW_so/dF or not
!   Coord     Atomic Cartesian coordinates
!   NOff      offset of each primitive shell
!   MapA      map relationship of each shell
!   LQnm      L-quantum number of each shell
!   LPat      (x,y,z) patterns of Cartesian functions
!   FNorm     common normalization factors (alpha independent part)
!   facdat    it saves (N)!
!   EXX       primitive exponents
!   EXX34     EXX^(3/4)
!
! Output:
!   SMat      NDer-th derivative matrices
!
! Scratch:
!   pref      prefactors of T integrals
!   fexp      it saves (pi/gamma)^3/2 * exp[(-alpha*beta/gamma)*Rab^2]
!   Ra,Rb     Positions of the centers (atoms) a and b
!   AB        Ra - Rb
!   pab       powers of PA and PB: 0,...,Lmu or Lnu + 2*abs(LSDD)
!   f2If      it saves (2i-1)!!
!   fIoG      it saves (2i-1)!!/(2gamma)^i
!   fClv      it saves c(l,v)=l!/[v!(l-v)!]
!   Gexp      alpha, beta, and gamma related values
!   ddc       1st and second derivatives of fexp.
!   tmp       tmp. values
!   Slmn      integrals for a given (L,M,N)
!   Sshl      saves integrals for S, P, D, F, ..., shells.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DDIPmain(MaxTyp,MaxXYZ,MaxSum,NSHELL,NTT,NDer,NMat,LSDD,IfSO,Coord,NOff,MapA,LQnm,LPat,pref,FNorm,facdat,EXX,EXX34, &
  fexp,RA,RB,AB,pab,f2If,fIoG,fClv,Gexp,ddc,tmp,Slmn,Sshl,SMat)
  use Constant, only : pi, Zero, One, Two
  implicit real(kind=8) (a-h,o-z)
  dimension Coord(3,*),NOff(*),MapA(*),LQnm(*),LPat(*),pref(*),FNorm(*),facdat(*),EXX(*),EXX34(*),RA(3),RB(3),AB(3),pab(*),   &
    f2If(*),fIoG(*),fClv(*),Gexp(*),ddc(*),tmp(*),Slmn(NMat,MaxXYZ,MaxXYZ),Sshl(NMat,MaxSum,MaxSum),SMat(NTT,NMat)
  logical IfSO

  LDD = LSDD/2
  MaxI1 = MaxTyp + LDD
  MaxI2 = MaxI1  + LDD

  ! (2n-1)!!, n=0,1,...,MaxI1
  do i=1,MaxI1+1
    f2If(i)=Tlm1ff(i)
  end do

  ! c(l,v)
  Call CombiFac(MaxI2,fClv,facdat)

  loop_ISHELL: DO ISHELL=1,NSHELL
    ! IA center to which ISHELL is attached
    IA=MapA(ISHELL)
    Ra(1) = Coord(1,IA)
    Ra(2) = Coord(2,IA)
    Ra(3) = Coord(3,IA)
    ITYPE = LQnm(ISHELL)
    IST = NOff(ISHELL)
    IStart=(ITYPE+0)*(ITYPE+1)*(ITYPE+2)/6+1
    IEnd=(ITYPE+1)*(ITYPE+2)*(ITYPE+3)/6

    Gexp(1)=EXX(ISHELL)        ! (1) Alpha
    Gexp(4)=EXX34(ISHELL)      ! (4) Alpha^3/4

    if(LSDD > 1) pref(5) = Gexp(1) + Gexp(1)   ! 2 * alpha

    loop_JSHELL: DO JSHELL=1,ISHELL
      ! JA center to which JSHELL is attached
      JA=MapA(JSHELL)
      ! If IA=JA, the derivatives are zero, i.e., dS/dIA + dS/dJA = 0 and so on.
      ! In this subroutine only the first part is calculated.
      if(IA == JA .and. NDer > 0)cycle

      Rb(1) = Coord(1,JA)
      Rb(2) = Coord(2,JA)
      Rb(3) = Coord(3,JA)
      JTYPE = LQnm(JSHELL)
      JST = NOff(JSHELL)
      JStart=(JTYPE+0)*(JTYPE+1)*(JTYPE+2)/6+1
      JEnd=(JTYPE+1)*(JTYPE+2)*(JTYPE+3)/6

      IMJ=ABS(ISHELL-JSHELL)

      ! Rab
      AB(1)=Ra(1)-Rb(1)
      AB(2)=Ra(2)-Rb(2)
      AB(3)=Ra(3)-Rb(3)

      Gexp(2)=EXX(JSHELL)        ! (2) Beta
      Gexp(5)=EXX34(JSHELL)      ! (5) Beta^3/4

      ! Gamma and other aa,bb related values
      Gexp(3)=Gexp(1)+Gexp(2)                                   ! (3) Gamma

      if(NDer > 0)then
        Gexp(6)=Gexp(1)/Gexp(3)                                 ! (6) Alpha/Gamma
        Gexp(7)=Gexp(2)/Gexp(3)                                 ! (7) Beta/Gamma
        Gexp(8)=Two*Gexp(1)*Gexp(2)/Gexp(3)                     ! (8) 2*Alpha*Beta/Gamma
      end if
      if(NDer > 1)then
        Gexp(9)=Gexp(6)*Gexp(6)                                 ! (9) (6)*(6)
        Gexp(10)=Gexp(7)*Gexp(7)                                ! (10) (7)*(7)
        Gexp(11)=Two*Gexp(6)*Gexp(7)                            ! (11) 2*(6)*(7)
      end if

      ! (2i-1)!!/(2gamma)^i
      do i=0,MaxI1
        fIoG(i+1)=f2If(i+1)*(Two*Gexp(3))**(-i)
      end do

      ! powers of PA and PB
      call PowerPAB(.false.,MaxI2,ITYPE+LDD,JTYPE+LDD,Gexp(1),Gexp(2),Gexp(3),Ra,Rb,tmp,pab)

      ! (pi/gamma)^3/2 * exp[(-alpha*beta/gamma)*Rab^2]
      tmp(1) = sqrt(pi/Gexp(3))
      tmp(2) =-Gexp(1)*Gexp(2)/Gexp(3)
      tmp(3) = AB(1) * AB(1) + AB(2) * AB(2) + AB(3) * AB(3)
      fexp = tmp(1)*tmp(1)*tmp(1) * exp(tmp(2)*tmp(3))
      ! Derivatives of fexp
      if(NDer > 0)then
        ddc(1)=-AB(1)*Gexp(8)
        ddc(2)=-AB(2)*Gexp(8)
        ddc(3)=-AB(3)*Gexp(8)
      end if
      if(NDer > 1)then
        ddc(4)=ddc(1)*ddc(1)-Gexp(8)
        ddc(5)=ddc(2)*ddc(2)-Gexp(8)
        ddc(6)=ddc(3)*ddc(3)-Gexp(8)
      end if

      ! prefactors
      if(LSDD > 0)then
        pref(1) = (Gexp(1)*Ra(1) + Gexp(2)*Rb(1)) / Gexp(3)  ! PC (= PO)
        pref(2) = (Gexp(1)*Ra(2) + Gexp(2)*Rb(2)) / Gexp(3)
        pref(3) = (Gexp(1)*Ra(3) + Gexp(2)*Rb(3)) / Gexp(3)
        pref(4) = 0.5d0 / Gexp(3)      ! 1/(2p)
      end if
      if(LSDD > 1) pref(6) = Gexp(2) + Gexp(2)   ! 2 * beta

      ! clean S
      Slmn=Zero

      ! Integrals calculation
      call DDIPmunu(MaxTyp,MaxXYZ,MaxI2,LSDD,IfSO,NDer,NMat,LPat,pref,Gexp,ITYPE,JTYPE,fClv,pab,fIoG,ddc,tmp,Slmn)

      ! Slmn is scaled by fexp, and saved to Sshl.
      do I=IStart,IEnd
        do J=JStart,JEnd
          do k=1,NMat
            Sshl(K,J,I) = fexp * Slmn(K,J-JStart+1,I-IStart+1)
          end do
        end do
      end do

      ! calculate SMat matrices
      Do I = IStart,IEnd
        MU=IST+I-IStart+1
        MUS=(MU*(MU-1))/2
        ! tmp(1) and tmp(2): normalization factors
        tmp(1)=FNorm(I)*Gexp(4)
        If(I > 1)tmp(1)=tmp(1)*Sqrt(Gexp(1)**ITYPE)
        JND=JEnd
        If(IMJ == 0) JND=I
        Do J =JStart,JND
          NU=JST+J-JStart+1
          MUNU=MUS+NU
          tmp(2)=FNorm(J)*Gexp(5)
          If(J > 1)tmp(2)=tmp(2)*Sqrt(Gexp(2)**JTYPE)
          tmp(3)=tmp(1)*tmp(2)
          do K=1,NMat
            SMat(MUNU,K)=SMat(MUNU,K)+Sshl(K,J,I)*tmp(3)
          end do
        end do
      end do

    end do loop_JSHELL
  end do loop_ISHELL

  Return
End subroutine DDIPmain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the integral matrices of S, dV/dF, or dW/dF between two primitive Cartesian Gaussian functions.
!
! Input:
!   LSDD:     0 for S, 1 for dV/dF, 2 for dW/dF
!   IfSO      whether compute dW_so/dF or not
!   MaxTyp:   Max-L. see commonb.inc.
!   MaxXYZ:   (MaxTyp+1)*(MaxTyp+2)/2
!   LQmu,LQnu:quantum numbers of Cartesian Gaussian functions <mu|
!             and |nu>, e.g. s=0, p=1, d=2, ...
!             NOTE: sp has been split into s and p
!   NDer:     the order of derivative
!   NMat:     the number of matrices
!   LPat:     (lx,ly,lz) patterns
!
! Output:
!   S:        integrals for given (L,M,N)
!
! Scratch:
!   pref:     prefactors
!   Gexp:     alpha, beta, and gamma related values
!   pab:      powers of PA and PB: 0,...,Lmu or Lnu + 2*abs(LDD)
!   fIoG:     it saves (2i-1)!!/(2gamma)^i
!   fClv:     it saves c(l,v)=l!/[v!(l-v)!]
!   ddc:      1st and second derivatives of fexp.
!   tmp:      tmp. values
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DDIPmunu(MaxTyp,MaxXYZ,MaxI2,LSDD,IfSO,NDer,NMat,LPat,pref,Gexp, LQmu,LQnu,  fClv,pab,fIoG,ddc,tmp,S)
  implicit real(kind=8) (a-h,o-z)
  dimension LPat(*), pref(*),pab(MaxI2+1,2,3),fIoG(*),Gexp(*),ddc(*),fClv(MaxI2+1,*),tmp(3,*),S(NMat,MaxXYZ,MaxXYZ)
  logical IfSO

! do loop in cartesian shells of mu and nu
  NKK=(LQmu+1)*(LQmu+2)/2
  NLL=(LQnu+1)*(LQnu+2)/2

  do KK=1,NKK
    call LPatt(LQmu, KK, LPat, La, Ma, Na)

    do LL=1,NLL
      call LPatt(LQnu, LL, LPat, Lb, Mb, Nb)

      if(LSDD == 0)then
!       for S
!       1) <0|0>    : <La|Lb>, <Ma|Mb>, <Na|Nb>
        tmp(1,1) = S0Int(La,  Lb,  fClv(1,La+1),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,2))
        tmp(2,1) = S0Int(Ma,  Mb,  fClv(1,Ma+1),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,2))
        tmp(3,1) = S0Int(Na,  Nb,  fClv(1,Na+1),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,2))
        if(NDer == 0) then
          S(1,LL,KK) = S(1,LL,KK) + tmp(1,1) * tmp(2,1) * tmp(3,1)
          cycle
        end if
!       2) <0|0>'   :
        tmp(1,2) = S1Int(La,  Lb,  fClv(1,La+1),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,Gexp,tmp(1,3))
        tmp(2,2) = S1Int(Ma,  Mb,  fClv(1,Ma+1),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,Gexp,tmp(1,3))
        tmp(3,2) = S1Int(Na,  Nb,  fClv(1,Na+1),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,Gexp,tmp(1,3))
        S(1,LL,KK) = S(1,LL,KK) - (ddc(1) * tmp(1,1) + tmp(1,2)) * tmp(2,1) * tmp(3,1)
        S(2,LL,KK) = S(2,LL,KK) - (ddc(2) * tmp(2,1) + tmp(2,2)) * tmp(1,1) * tmp(3,1)
        S(3,LL,KK) = S(3,LL,KK) - (ddc(3) * tmp(3,1) + tmp(3,2)) * tmp(1,1) * tmp(2,1)
      else if(LSDD == 1)then
!       for D = dV/dF
!       1) <0|0>    : <La|Lb>, <Ma|Mb>, <Na|Nb>
        tmp(1,1) = S0Int(La,  Lb,  fClv(1,La+1),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,4))
        tmp(2,1) = S0Int(Ma,  Mb,  fClv(1,Ma+1),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,4))
        tmp(3,1) = S0Int(Na,  Nb,  fClv(1,Na+1),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,4))
!       2) <-1|0>   : <La-1|Lb>, <Ma-1|Mb>, <Na-1|Nb>
        tmp(1,2) = S0Int(La-1,Lb,  fClv(1,La  ),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,4))
        tmp(2,2) = S0Int(Ma-1,Mb,  fClv(1,Ma  ),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,4))
        tmp(3,2) = S0Int(Na-1,Nb,  fClv(1,Na  ),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,4))
!       3) <0|-1>   : <La|Lb-1>, <Ma|Mb-1>, <Na|Nb-1>
        tmp(1,3) = S0Int(La,  Lb-1,fClv(1,La+1),fClv(1,Lb  ), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,4))
        tmp(2,3) = S0Int(Ma,  Mb-1,fClv(1,Ma+1),fClv(1,Mb  ), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,4))
        tmp(3,3) = S0Int(Na,  Nb-1,fClv(1,Na+1),fClv(1,Nb  ), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,4))
!       4) {0|0} = <0|0> + <-1|0> + <0|-1>
        tmp(1,4) = pref(1) * tmp(1,1) + La * pref(4) * tmp(1,2) + Lb * pref(4) * tmp(1,3)
        tmp(2,4) = pref(2) * tmp(2,1) + Ma * pref(4) * tmp(2,2) + Mb * pref(4) * tmp(2,3)
        tmp(3,4) = pref(3) * tmp(3,1) + Na * pref(4) * tmp(3,2) + Nb * pref(4) * tmp(3,3)
        if(NDer == 0) then
          S(1,LL,KK) = S(1,LL,KK) + tmp(1,4) * tmp(2,1) * tmp(3,1)
          S(2,LL,KK) = S(2,LL,KK) + tmp(1,1) * tmp(2,4) * tmp(3,1)
          S(3,LL,KK) = S(3,LL,KK) + tmp(1,1) * tmp(2,1) * tmp(3,4)
          cycle
        end if
!       5) <0|0>'   :
        tmp(1,5) = S1Int(La,  Lb,  fClv(1,La+1),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,Gexp,tmp(1,8))
        tmp(2,5) = S1Int(Ma,  Mb,  fClv(1,Ma+1),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,Gexp,tmp(1,8))
        tmp(3,5) = S1Int(Na,  Nb,  fClv(1,Na+1),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,Gexp,tmp(1,8))
!       6) <-1|0>'  :
        tmp(1,6) = S1Int(La-1,Lb,  fClv(1,La  ),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,Gexp,tmp(1,8))
        tmp(2,6) = S1Int(Ma-1,Mb,  fClv(1,Ma  ),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,Gexp,tmp(1,8))
        tmp(3,6) = S1Int(Na-1,Nb,  fClv(1,Na  ),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,Gexp,tmp(1,8))
!       7) <0|-1>'  :
        tmp(1,7) = S1Int(La,  Lb-1,fClv(1,La+1),fClv(1,Lb  ), pab(1,1,1),pab(1,2,1),fIoG,Gexp,tmp(1,8))
        tmp(2,7) = S1Int(Ma,  Mb-1,fClv(1,Ma+1),fClv(1,Mb  ), pab(1,1,2),pab(1,2,2),fIoG,Gexp,tmp(1,8))
        tmp(3,7) = S1Int(Na,  Nb-1,fClv(1,Na+1),fClv(1,Nb  ), pab(1,1,3),pab(1,2,3),fIoG,Gexp,tmp(1,8))
!       8) {0|0}'   :
        tmp(1,8) = Gexp(6) * tmp(1,1) + pref(1) * tmp(1,5) + La * pref(4) * tmp(1,6) + Lb * pref(4) * tmp(1,7)
        tmp(2,8) = Gexp(6) * tmp(2,1) + pref(2) * tmp(2,5) + Ma * pref(4) * tmp(2,6) + Mb * pref(4) * tmp(2,7)
        tmp(3,8) = Gexp(6) * tmp(3,1) + pref(3) * tmp(3,5) + Na * pref(4) * tmp(3,6) + Nb * pref(4) * tmp(3,7)
!       dDx/dx, dDy/dx, dDz/dx
        S(1,LL,KK) = S(1,LL,KK) + (ddc(1) * tmp(1,4) + tmp(1,8)) * tmp(2,1) * tmp(3,1)
        S(2,LL,KK) = S(2,LL,KK) + (ddc(1) * tmp(1,1) + tmp(1,5)) * tmp(2,4) * tmp(3,1)
        S(3,LL,KK) = S(3,LL,KK) + (ddc(1) * tmp(1,1) + tmp(1,5)) * tmp(2,1) * tmp(3,4)
!       dDx/dy, dDy/dy, dDz/dy
        S(4,LL,KK) = S(4,LL,KK) + (ddc(2) * tmp(2,1) + tmp(2,5)) * tmp(1,4) * tmp(3,1)
        S(5,LL,KK) = S(5,LL,KK) + (ddc(2) * tmp(2,4) + tmp(2,8)) * tmp(1,1) * tmp(3,1)
        S(6,LL,KK) = S(6,LL,KK) + (ddc(2) * tmp(2,1) + tmp(2,5)) * tmp(1,1) * tmp(3,4)
!       dDx/dz, dDy/dz, dDz/dz
        S(7,LL,KK) = S(7,LL,KK) + (ddc(3) * tmp(3,1) + tmp(3,5)) * tmp(1,4) * tmp(2,1)
        S(8,LL,KK) = S(8,LL,KK) + (ddc(3) * tmp(3,1) + tmp(3,5)) * tmp(1,1) * tmp(2,4)
        S(9,LL,KK) = S(9,LL,KK) + (ddc(3) * tmp(3,4) + tmp(3,8)) * tmp(1,1) * tmp(2,1)
      else if(LSDD == 2)then
!       for dW/dF
!       1) <0|+1>   : <La|Lb+1>, <Ma|Mb+1>, <Na|Nb+1>
        tmp(1, 1) = S0Int(La,  Lb+1,fClv(1,La+1),fClv(1,Lb+2), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 1) = S0Int(Ma,  Mb+1,fClv(1,Ma+1),fClv(1,Mb+2), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 1) = S0Int(Na,  Nb+1,fClv(1,Na+1),fClv(1,Nb+2), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       2) <0|0>    : <La|Lb>, <Ma|Mb>, <Na|Nb>
        tmp(1, 2) = S0Int(La,  Lb,  fClv(1,La+1),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 2) = S0Int(Ma,  Mb,  fClv(1,Ma+1),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 2) = S0Int(Na,  Nb,  fClv(1,Na+1),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       3) <0|-1>   : <La|Lb-1>, <Ma|Mb-1>, <Na|Nb-1>
        tmp(1, 3) = S0Int(La,  Lb-1,fClv(1,La+1),fClv(1,Lb  ), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 3) = S0Int(Ma,  Mb-1,fClv(1,Ma+1),fClv(1,Mb  ), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 3) = S0Int(Na,  Nb-1,fClv(1,Na+1),fClv(1,Nb  ), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       4) <0|-2>   :
!       5) <-1|+1>  : <La-1|Lb+1>, <Ma-1|Mb+1>, <Na-1|Nb+1>
        tmp(1, 5) = S0Int(La-1,Lb+1,fClv(1,La  ),fClv(1,Lb+2), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 5) = S0Int(Ma-1,Mb+1,fClv(1,Ma  ),fClv(1,Mb+2), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 5) = S0Int(Na-1,Nb+1,fClv(1,Na  ),fClv(1,Nb+2), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       6) <-1|0>   : <La-1|Lb>, <Ma-1|Mb>, <Na-1|Nb>
        tmp(1, 6) = S0Int(La-1,Lb,  fClv(1,La  ),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 6) = S0Int(Ma-1,Mb,  fClv(1,Ma  ),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 6) = S0Int(Na-1,Nb,  fClv(1,Na  ),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       7) <-1|-1>  : <La-1|Lb-1>, <Ma-1|Mb-1>, <Na-1|Nb-1>
        tmp(1, 7) = S0Int(La-1,Lb-1,fClv(1,La  ),fClv(1,Lb  ), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 7) = S0Int(Ma-1,Mb-1,fClv(1,Ma  ),fClv(1,Mb  ), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 7) = S0Int(Na-1,Nb-1,fClv(1,Na  ),fClv(1,Nb  ), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       8) <-1|-2>  : <La-1|Lb-2>, <Ma-1|Mb-2>, <Na-1|Nb-2>
        tmp(1, 8) = S0Int(La-1,Lb-2,fClv(1,La  ),fClv(1,Lb-1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 8) = S0Int(Ma-1,Mb-2,fClv(1,Ma  ),fClv(1,Mb-1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 8) = S0Int(Na-1,Nb-2,fClv(1,Na  ),fClv(1,Nb-1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       9) <+1|+1>  : <La+1|Lb+1>, <Ma+1|Mb+1>, <Na+1|Nb+1>
        tmp(1, 9) = S0Int(La+1,Lb+1,fClv(1,La+2),fClv(1,Lb+2), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 9) = S0Int(Ma+1,Mb+1,fClv(1,Ma+2),fClv(1,Mb+2), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 9) = S0Int(Na+1,Nb+1,fClv(1,Na+2),fClv(1,Nb+2), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       10) <+1|0>  : <La+1|Lb>, <Ma+1|Mb>, <Na+1|Nb>
        tmp(1,10) = S0Int(La+1,Lb  ,fClv(1,La+2),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2,10) = S0Int(Ma+1,Mb  ,fClv(1,Ma+2),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3,10) = S0Int(Na+1,Nb  ,fClv(1,Na+2),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       11) <+1|-1> : <La+1|Lb-1>, <Ma+1|Mb-1>, <Na+1|Nb-1>
        tmp(1,11) = S0Int(La+1,Lb-1,fClv(1,La+2),fClv(1,Lb  ), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2,11) = S0Int(Ma+1,Mb-1,fClv(1,Ma+2),fClv(1,Mb  ), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3,11) = S0Int(Na+1,Nb-1,fClv(1,Na+2),fClv(1,Nb  ), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       12) <+1|-2> : <La+1|Lb-2>, <Ma+1|Mb-2>, <Na+1|Nb-2>
        tmp(1,12) = S0Int(La+1,Lb-2,fClv(1,La+2),fClv(1,Lb-1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2,12) = S0Int(Ma+1,Mb-2,fClv(1,Ma+2),fClv(1,Mb-1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3,12) = S0Int(Na+1,Nb-2,fClv(1,Na+2),fClv(1,Nb-1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       13) <-2|+1> : <La-2|Lb+1>, <Ma-2|Mb+1>, <Na-2|Nb+1>
        tmp(1,13) = S0Int(La-2,Lb+1,fClv(1,La-1),fClv(1,Lb+2), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2,13) = S0Int(Ma-2,Mb+1,fClv(1,Ma-1),fClv(1,Mb+2), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3,13) = S0Int(Na-2,Nb+1,fClv(1,Na-1),fClv(1,Nb+2), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       14) <-2|0>  :
!       15) <-2|-1> : <La-2|Lb-1>, <Ma-2|Mb-1>, <Na-2|Nb-1>
        tmp(1,15) = S0Int(La-2,Lb-1,fClv(1,La-1),fClv(1,Lb  ), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2,15) = S0Int(Ma-2,Mb-1,fClv(1,Ma-1),fClv(1,Mb  ), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3,15) = S0Int(Na-2,Nb-1,fClv(1,Na-1),fClv(1,Nb  ), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       16) {0|+1}
!       17) {0|0} = <0|0> + <-1|0> + <0|-1>
        tmp(1,17) = pref(1) * tmp(1, 2) +  La    * pref(4) * tmp(1, 6) +  Lb    * pref(4) * tmp(1, 3)
        tmp(2,17) = pref(2) * tmp(2, 2) +  Ma    * pref(4) * tmp(2, 6) +  Mb    * pref(4) * tmp(2, 3)
        tmp(3,17) = pref(3) * tmp(3, 2) +  Na    * pref(4) * tmp(3, 6) +  Nb    * pref(4) * tmp(3, 3)
!       18) {0|-1}
!       19) {+1|+1} = <+1|+1> + <0|+1> + <+1|0>
        tmp(1,19) = pref(1) * tmp(1, 9) + (La+1) * pref(4) * tmp(1, 1) + (Lb+1) * pref(4) * tmp(1,10)
        tmp(2,19) = pref(2) * tmp(2, 9) + (Ma+1) * pref(4) * tmp(2, 1) + (Mb+1) * pref(4) * tmp(2,10)
        tmp(3,19) = pref(3) * tmp(3, 9) + (Na+1) * pref(4) * tmp(3, 1) + (Nb+1) * pref(4) * tmp(3,10)
!       20) {+1|0}
!       21) {+1|-1} = <+1|-1> + <0|-1> + <+1|-2>
        tmp(1,21) = pref(1) * tmp(1,11) + (La+1) * pref(4) * tmp(1, 3) + (Lb-1) * pref(4) * tmp(1,12)
        tmp(2,21) = pref(2) * tmp(2,11) + (Ma+1) * pref(4) * tmp(2, 3) + (Mb-1) * pref(4) * tmp(2,12)
        tmp(3,21) = pref(3) * tmp(3,11) + (Na+1) * pref(4) * tmp(3, 3) + (Nb-1) * pref(4) * tmp(3,12)
!       22) {-1|+1} = <-1|+1> + <-2|+1> + <-1|0>
        tmp(1,22) = pref(1) * tmp(1, 5) + (La-1) * pref(4) * tmp(1,13) + (Lb+1) * pref(4) * tmp(1, 6)
        tmp(2,22) = pref(2) * tmp(2, 5) + (Ma-1) * pref(4) * tmp(2,13) + (Mb+1) * pref(4) * tmp(2, 6)
        tmp(3,22) = pref(3) * tmp(3, 5) + (Na-1) * pref(4) * tmp(3,13) + (Nb+1) * pref(4) * tmp(3, 6)
!       23) {-1|0}
!       24) {-1|-1} = <-1|-1> + <-2|-1> + <-1|-2>
        tmp(1,24) = pref(1) * tmp(1, 7) + (La-1) * pref(4) * tmp(1,15) + (Lb-1) * pref(4) * tmp(1, 8)
        tmp(2,24) = pref(2) * tmp(2, 7) + (Ma-1) * pref(4) * tmp(2,15) + (Mb-1) * pref(4) * tmp(2, 8)
        tmp(3,24) = pref(3) * tmp(3, 7) + (Na-1) * pref(4) * tmp(3,15) + (Nb-1) * pref(4) * tmp(3, 8)
!       tmp(:,25)
        tmp(1,25) = pref(5) * pref(6) * tmp(1,19)   - pref(5) * Lb     * tmp(1,21)    &
                  - La      * pref(6) * tmp(1,22)   + La      * Lb     * tmp(1,24)
        tmp(2,25) = pref(5) * pref(6) * tmp(2,19)   - pref(5) * Mb     * tmp(2,21)    &
                  - Ma      * pref(6) * tmp(2,22)   + Ma      * Mb     * tmp(2,24)
        tmp(3,25) = pref(5) * pref(6) * tmp(3,19)   - pref(5) * Nb     * tmp(3,21)    &
                  - Na      * pref(6) * tmp(3,22)   + Na      * Nb     * tmp(3,24)
!       tmp(:,26)
        tmp(1,26) = pref(5) * pref(6) * tmp(1, 9)   - pref(5) * Lb     * tmp(1,11)    &
                  - La      * pref(6) * tmp(1, 5)   + La      * Lb     * tmp(1, 7)
        tmp(2,26) = pref(5) * pref(6) * tmp(2, 9)   - pref(5) * Mb     * tmp(2,11)    &
                  - Ma      * pref(6) * tmp(2, 5)   + Ma      * Mb     * tmp(2, 7)
        tmp(3,26) = pref(5) * pref(6) * tmp(3, 9)   - pref(5) * Nb     * tmp(3,11)    &
                  - Na      * pref(6) * tmp(3, 5)   + Na      * Nb     * tmp(3, 7)
!       dW0/dF
        S(1,LL,KK) = S(1,LL,KK) + tmp(1,25) * tmp(2, 2) * tmp(3, 2) + tmp(1,17) * (tmp(2,26) * tmp(3,2) + tmp(2,2) * tmp(3,26))
        S(2,LL,KK) = S(2,LL,KK) + tmp(1, 2) * tmp(2,25) * tmp(3, 2) + tmp(2,17) * (tmp(1,26) * tmp(3,2) + tmp(1,2) * tmp(3,26))
        S(3,LL,KK) = S(3,LL,KK) + tmp(1, 2) * tmp(2, 2) * tmp(3,25) + tmp(3,17) * (tmp(1,26) * tmp(2,2) + tmp(1,2) * tmp(2,26))

        if(.NOT. IfSO) cycle

!       For dWso/dF only.
!       4) <0|-2>   : <La|Lb-2>, <Ma|Mb-2>, <Na|Nb-2>
        tmp(1, 4) = S0Int(La,  Lb-2,fClv(1,La+1),fClv(1,Lb-1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2, 4) = S0Int(Ma,  Mb-2,fClv(1,Ma+1),fClv(1,Mb-1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3, 4) = S0Int(Na,  Nb-2,fClv(1,Na+1),fClv(1,Nb-1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       14) <-2|0>  : <La-2|Lb>, <Ma-2|Mb>, <Na-2|Nb>
        tmp(1,14) = S0Int(La-2,Lb  ,fClv(1,La-1),fClv(1,Lb+1), pab(1,1,1),pab(1,2,1),fIoG,tmp(1,28))
        tmp(2,14) = S0Int(Ma-2,Mb  ,fClv(1,Ma-1),fClv(1,Mb+1), pab(1,1,2),pab(1,2,2),fIoG,tmp(1,28))
        tmp(3,14) = S0Int(Na-2,Nb  ,fClv(1,Na-1),fClv(1,Nb+1), pab(1,1,3),pab(1,2,3),fIoG,tmp(1,28))
!       16) {0|+1} = <0|+1> + <-1|+1> + <0|0>
        tmp(1,16) = pref(1) * tmp(1, 1) +  La    * pref(4) * tmp(1, 5) + (Lb+1) * pref(4) * tmp(1, 2)
        tmp(2,16) = pref(2) * tmp(2, 1) +  Ma    * pref(4) * tmp(2, 5) + (Mb+1) * pref(4) * tmp(2, 2)
        tmp(3,16) = pref(3) * tmp(3, 1) +  Na    * pref(4) * tmp(3, 5) + (Nb+1) * pref(4) * tmp(3, 2)
!       18) {0|-1} = <0|-1> + <-1|-1> + <0|-2>
        tmp(1,18) = pref(1) * tmp(1, 3) +  La    * pref(4) * tmp(1, 7) + (Lb-1) * pref(4) * tmp(1, 4)
        tmp(2,18) = pref(2) * tmp(2, 3) +  Ma    * pref(4) * tmp(2, 7) + (Mb-1) * pref(4) * tmp(2, 4)
        tmp(3,18) = pref(3) * tmp(3, 3) +  Na    * pref(4) * tmp(3, 7) + (Nb-1) * pref(4) * tmp(3, 4)
!       20) {+1|0} = <+1|0> + <0|0> + <+1|-1>
        tmp(1,20) = pref(1) * tmp(1,10) + (La+1) * pref(4) * tmp(1, 2) +  Lb    * pref(4) * tmp(1,11)
        tmp(2,20) = pref(2) * tmp(2,10) + (Ma+1) * pref(4) * tmp(2, 2) +  Mb    * pref(4) * tmp(2,11)
        tmp(3,20) = pref(3) * tmp(3,10) + (Na+1) * pref(4) * tmp(3, 2) +  Nb    * pref(4) * tmp(3,11)
!       23) {-1|0} = <-1|0> + <-2|0> + <-1|-1>
        tmp(1,23) = pref(1) * tmp(1, 6) + (La-1) * pref(4) * tmp(1,14) +  Lb    * pref(4) * tmp(1, 7)
        tmp(2,23) = pref(2) * tmp(2, 6) + (Ma-1) * pref(4) * tmp(2,14) +  Mb    * pref(4) * tmp(2, 7)
        tmp(3,23) = pref(3) * tmp(3, 6) + (Na-1) * pref(4) * tmp(3,14) +  Nb    * pref(4) * tmp(3, 7)
!       tmp(:,25)
        tmp(1,25) = pref(5) * tmp(1,20) - La * tmp(1,23)
        tmp(2,25) = pref(5) * tmp(2,20) - Ma * tmp(2,23)
        tmp(3,25) = pref(5) * tmp(3,20) - Na * tmp(3,23)
!       tmp(:,26)
        tmp(1,26) = pref(6) * tmp(1,16) - Lb * tmp(1,18)
        tmp(2,26) = pref(6) * tmp(2,16) - Mb * tmp(2,18)
        tmp(3,26) = pref(6) * tmp(3,16) - Nb * tmp(3,18)
!       tmp(:,27)
        tmp(1,27) = pref(5) * tmp(1,10) - La * tmp(1, 6)
        tmp(2,27) = pref(5) * tmp(2,10) - Ma * tmp(2, 6)
        tmp(3,27) = pref(5) * tmp(3,10) - Na * tmp(3, 6)
!       tmp(:,28)
        tmp(1,28) = pref(6) * tmp(1, 1) - Lb * tmp(1, 3)
        tmp(2,28) = pref(6) * tmp(2, 1) - Mb * tmp(2, 3)
        tmp(3,28) = pref(6) * tmp(3, 1) - Nb * tmp(3, 3)
!       dWx/dF
        S( 4,LL,KK) = S( 4,LL,KK) + tmp(1,17) * (tmp(2,27) * tmp(3,28) - tmp(2,28) * tmp(3,27))
        S( 5,LL,KK) = S( 5,LL,KK) + tmp(1, 2) * (tmp(2,25) * tmp(3,28) - tmp(2,26) * tmp(3,27))
        S( 6,LL,KK) = S( 6,LL,KK) + tmp(1, 2) * (tmp(2,27) * tmp(3,26) - tmp(2,28) * tmp(3,25))
!       dWy/dF
        S( 7,LL,KK) = S( 7,LL,KK) + tmp(2, 2) * (tmp(1,26) * tmp(3,27) - tmp(1,25) * tmp(3,28))
        S( 8,LL,KK) = S( 8,LL,KK) + tmp(2,17) * (tmp(1,28) * tmp(3,27) - tmp(1,27) * tmp(3,28))
        S( 9,LL,KK) = S( 9,LL,KK) + tmp(2, 2) * (tmp(1,28) * tmp(3,25) - tmp(1,27) * tmp(3,26))
!       dWz/dF
        S(10,LL,KK) = S(10,LL,KK) + tmp(3, 2) * (tmp(1,25) * tmp(2,28) - tmp(1,26) * tmp(2,27))
        S(11,LL,KK) = S(11,LL,KK) + tmp(3, 2) * (tmp(1,27) * tmp(2,26) - tmp(1,28) * tmp(2,25))
        S(12,LL,KK) = S(12,LL,KK) + tmp(3,17) * (tmp(1,27) * tmp(2,28) - tmp(1,28) * tmp(2,27))
      end if
    end do
  end do

  return
end subroutine DDIPmunu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine calculates the nuclear attraction integrals V and pVp in primitive basis functions.
! The basic method was described in
! H. Taketa, S. Huzinaga, and K. O-ohata, J. Phys. Soc. Jap. 21(11), 2313, 1966.
!
! Input:
!   LpVp      =  0   for V
!             =  1   for pVp_sf
!             =  2   for pVp_sf, pVp_x, pVp_y, and pVp_z
!   ISO0      =      compute only pVp_sf (1), pVp_x (2), pVp_y (3), or pVp_z (4)
!                    for LpVp = 2, and compute all four compounents by default.
!   NAI1c            One-center V or pVp if non-zero
!   IRMSIA    =  0   V or pVp
!             >  0   d(V)/dR0 or d(pVp)/dR0 for atom IRMSIA
!             <  0   d(V)/dQ or d(pVp)/dQ for atom |IRMSIA|
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   natom          : number of atoms
!   za             : atomic nuclear charges
!   xyzmol         : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   nbasis         : number of primitive Cartesian functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   R0             : Nuclear RMS charge radii
!
! Output:
!   DMat(NTT*nm)     NAI matrices where nm = 1 (LpVp /= 2) or 4 (LpVp = 2).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine pnai_drv(iout,LpVp,ISO0,NAI1c,IRMSIA,mpatt,natom,za,xyzmol,maxtyp,nshell,nbasis,mxla,nshlla,gtoexp,R0,DMat,V,Mem)
  use Constant, only : maxlq, Zero
  implicit real(kind=8) (a-h,o-z)
  dimension :: za(*), xyzmol(*), mxla(natom), nshlla(0:maxlq,natom), gtoexp(*), R0(*), DMat(*), V(Mem)
  allocatable :: K(:), LPat(:)

! check
  if(IRMSIA == 0)then
    NDerRQ = 0
    NEla = 1
    NElb = 1
  else if(IRMSIA > 0 .and. IRMSIA <= natom)then
    NDerRQ = 1
    NEla = 1
    NElb = 1
  else if(IRMSIA < 0 .and. -IRMSIA <= natom)then
    NDerRQ = 2
    NEla = 6
    NElb = 4
  else
    write(iout,"(/,' Error in sub. pnai_drv: IRMSIA is out of range.')")
    call estop(1)
  end if
  if(IRMSIA /= 0 .and. NAI1c /= 0)then  ! for X2C-ARX and X2C-AU only
    write(iout,"(/,' Error in sub. pnai_drv: One-center dNAI does not work.')")
    call estop(1)
  end if

  if(LpVp == 0) then
    nm1   = 1
    nm2   = 1
    LpVp0 = 0
    MaxLI = MaxTyp
  else if(LpVp == 1) then
    nm1   = 1
    nm2   = 1
    LpVp0 = 1
    MaxLI = MaxTyp + 1
  else if(LpVp == 2) then
    nm1   = 1
    nm2   = 4
    if(ISO0 > 0 .and. ISO0 < 5) then
      nm1   = ISO0
      nm2   = ISO0
    end if
    LpVp0 = 1
    MaxLI = MaxTyp + 1
  else
    write(iout,"(/,' Error in sub. pnai_drv: Unknown LpVp = ',i2)") LpVp
    call estop(1)
  end if

! ==================
! | initialization |
! ==================
  MaxLD=2*MaxLI
  MaxSum=(MaxTyp+1)*(MaxTyp+2)*(MaxTyp+3)/6
  MaxSzA=(MaxLD+1)*(MaxLI+1)*(MaxLI+1)
  MaxXYZ=(MaxTyp+1)*(MaxTyp+2)/2
  MaxXYZ2=MaxXYZ*MaxXYZ
  NPatt =(MaxTyp+1)*(MaxTyp+2)*(MaxTyp+3)*(MaxTyp*3+2)
  NPatt = NPatt/12
  Lgg = 18
  Ltmp= 47
  Lprf= 12
  Lfau= 56+NDerRQ+MaxLD
  NTT=nbasis*(nbasis+1)/2

  INrm = 1
  I00  = INrm + MaxSum
  Iprf = I00  + 1
  Iamu = Iprf + Lprf
  Ibnu = Iamu + 1
  IEPW = Ibnu + 1
  IRa  = IEPW + NSHELL
  IRb  = IRa  + 3
  IRc  = IRb  + 3
  Idv  = IRc  + 3
  Idvc = Idv  + MaxXYZ2*NEla
  Ipc  = Idvc + MaxSum*MaxSum*NEla
  Ipcp = Ipc  + 3
  Iab  = Ipcp + (MaxLD + 1)*3
  Ipab = Iab  + 3
  Ieps = Ipab + (MaxLD + 1)*6
  Iiru = Ieps + MaxLD + 1
  Ijsv = Iiru + MaxSzA*NElb
  Iktw = Ijsv + MaxSzA*NElb
  Igk  = Iktw + MaxSzA*NElb
  Igg  = Igk  + MaxLD + 1 + NDerRQ
  Ifdt = Igg  + Lgg
  Ifau = Ifdt + (MaxLD+1)
  Itmp = Ifau + Lfau
  MemV = Itmp + Ltmp - 1
  if(MemV > Mem) then
    write(iout,"(/,' Error in sub. pnai_drv: MemV > Mem.')")
    call estop(1)
  end if

  KOff = 1
  KMap = KOff + NSHELL
  KLqn = KMap + NSHELL
  MemK = KLqn + NSHELL

  allocate(K(MemK),LPat(MaxSum*3))
  V(1:MemV) = Zero
  K = 0

  ! eps
  V(I00)=EpsFun(30)

  ! N!
  call FacCal(MaxLD,V(Ifdt))

  ! Generate (x,y,z) patterns in MOLDEN or Gaussian mode
  call GnPatt(MaxTyp,mpatt,LPat)

  ! common normalization factors
  call NormFac(MaxTyp,LPat,V(INrm))

  ! primitive exponents EXX and common factors EXXP=sqrt(EXX^LQ)
  call ExpPsq(natom,mxla,nshlla,gtoexp,V(IEPW))

  ! setup indices of primitive shells
  call SetShl(natom,nshell,mxla,nshlla,K(KOff),K(KMap),K(KLqn))

  DMat(1:NTT*NEla*(nm2-nm1+1)) = Zero
  ioff = 1
  do i=nm1,nm2
    ISO = i-1
    if(NAI1c == 0) then
      call NAImain(MaxXYZ,MaxSum,LpVp0,ISO,NTT,natom,NSHELL,IRMSIA,NEla,R0,za,xyzmol,V(IRa),V(IRb),V(IRc),V(INrm),     &
        V(I00),V(Iprf),V(Idv),V(Idvc),V(Iamu),V(Ibnu),gtoexp,V(IEPW),V(Igg),V(Ieps),V(Iab),V(Ipab),V(Ipc),V(Ipcp),     &
        V(Ifdt),V(Iiru),V(Ijsv),V(Iktw),V(Igk),V(Ifau),V(Itmp),K(KOff),K(KMap),K(KLqn),LPat,DMat(ioff))
    else
      call NAI1Cmain(MaxXYZ,MaxSum,LpVp0,ISO,natom,NSHELL,R0,za,V(IRa),V(IRb),V(IRc),V(INrm),V(I00),V(Iprf),V(Idv),    &
        V(Idvc),V(Iamu),V(Ibnu),gtoexp,V(IEPW),V(Igg),V(Ieps),V(Iab),V(Ipab),V(Ipc),V(Ipcp),V(Ifdt),V(Iiru),V(Ijsv),   &
        V(Iktw),V(Igk),V(Ifau),V(Itmp),K(KOff),K(KMap),K(KLqn),LPat,DMat(ioff))
    end if
    ioff = ioff + NTT*NEla
  end do

  deallocate(K,LPat)

  return
end subroutine pnai_drv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! One-center version of subroutine NAImain. Note that NAI<La|Lb> = 0 if La /= Lb.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NAI1Cmain(MaxXYZ,MaxSum,LpVp,ISO,NAtoms,NSHELL,R0,ZA,Ra,Rb,Rc,FNorm,e0,prefw,dV,dVc,aa,bb,EXX,EXXP,gg,Epsp,AB,pab,  &
  PC,PCp,facdat,Airu,Ajsv,Aktw,Gkp,faux,tmp,NOff,MapA,LQnm,LPat,DMat)
  use Constant, only : Zero, One
  implicit real(kind=8) (a-h,o-z)
  parameter(Lgg = 18, Ltmp= 47, Lprf= 12)
  dimension R0(*),ZA(*),Ra(3),Rb(3),Rc(3),FNorm(*),prefw(Lprf),dV(MaxXYZ,MaxXYZ),dVc(MaxSum,MaxSum),EXX(*),EXXP(*),gg(Lgg),   &
    Epsp(*),AB(3),pab(*),PC(3),PCp(*),facdat(*),Airu(*),Ajsv(*),Aktw(*),Gkp(*),faux(*),tmp(Ltmp),NOff(*),MapA(*),LQnm(*),    &
    LPat(*),DMat(*)

  NDerRQ = 0
  NEle = 1

  Ra = Zero
  Rb = Zero
  Rc = Zero
  AB = Zero

!-LOOP OVER ISHELL
  DO ISHELL=1,NSHELL
    Iatm   = MapA(ISHELL)

    aa     = EXX(ISHELL)
    ITYPE  = LQnm(ISHELL)
    IStart = (ITYPE+0)*(ITYPE+1)*(ITYPE+2)/6+1
    IEnd   = (ITYPE+1)*(ITYPE+2)*(ITYPE+3)/6
    IST    = NOff(ISHELL)

!---LOOP OVER JSHELL
    DO JSHELL=1,ISHELL
      if(Iatm /= MapA(JSHELL)) cycle

      bb     = EXX(JSHELL)
      JTYPE  = LQnm(JSHELL)
      JStart = (JTYPE+0)*(JTYPE+1)*(JTYPE+2)/6+1
      JEnd   = (JTYPE+1)*(JTYPE+2)*(JTYPE+3)/6
      IMJ    = ABS(ISHELL-JSHELL)
      JST    = NOff(JSHELL)

      if(ITYPE == JTYPE) then
        LQSum=ITYPE+JTYPE+LpVp*2

        ! Gamma and other aa,bb related values
        call exponval(aa,bb,gg)
        ! powers of 1/(4*Gamma)
        call EpsPower(LQSum,gg(1),Epsp)
        ! powers of PA and PB (1-center)
        call PowerPAB(.true.,LQSum,ITYPE+LpVp,JTYPE+LpVp,aa,bb,gg,Ra,Rb,tmp,pab)
        ! Integrals calculation
        dV = Zero
        call NAImunu(LQSum,MaxXYZ,NDerRQ,LpVp,ISO,e0,NEle,dV,gg,Epsp,facdat,prefw,aa,Ra,ITYPE, bb,Rb,JTYPE, LPat,R0(Iatm),  &
          Rc,AB,pab,PC,PCp,Airu,Ajsv,Aktw,Gkp,faux,tmp)

        ! dV is scaled by (-Zc)*(2pi/Gamma)*exp(...), Zc > 0, and saved to dVc.
        tmp(2) = -ZA(Iatm)
        Call ZExABG(aa,bb,gg(1),tmp(2),AB,tmp(1),tmp(3))
        do I=IStart,IEnd
          do J=JStart,JEnd
            dVc(J,I) = tmp(1) * dV(J-JStart+1,I-IStart+1)
          end do
        end do

        ! calculate DMat
        Do I = IStart,IEnd
          MU=IST+I-IStart+1
          MUS=(MU*(MU-1))/2
          tmp(1)=FNorm(I)*gg(15)
          If(I > 1)tmp(1)=tmp(1)*EXXP(ISHELL)
          JND=JEnd
          If(IMJ == 0) JND=I
          Do J =JStart,JND
            NU=JST+J-JStart+1
            MUNU=MUS+NU
            tmp(2)=FNorm(J)*gg(16)
            If(J > 1)tmp(2)=tmp(2)*EXXP(JSHELL)
            tmp(3)=tmp(1)*tmp(2)
            DMat(MUNU)=DMat(MUNU)+dVc(J,I)*tmp(3)
          end do
        end do
      else
        ! DMat(La,Lb) = 0 if La /= Lb
        Do I = IStart,IEnd
          MU=IST+I-IStart+1
          MUS=(MU*(MU-1))/2
          JND=JEnd
          If(IMJ == 0) JND=I
          Do J =JStart,JND
            NU=JST+J-JStart+1
            MUNU=MUS+NU
            DMat(MUNU) = Zero
          end do
        end do
      end if

!---JSHELL
    end do

!-ISHELL
  end do

  return
end subroutine NAI1Cmain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Main subroutine to calculate the integral NAI<LQmu(a)|V(c)|LQnu(b)> between two primitive Cartesian Gaussian functions.
!
! Input:
!   LpVp:   = 0 for V
!           = 1 for pVp
!   ISO:    = 0 for V or pVp_sf
!           = 1 for pVp_x
!           = 2 for pVp_y
!           = 3 for pVp_z
!   NAtoms:   Number of atoms
!   NSHELL:   Number of primitive shells
!   NTT:      NTT=NBAS*(NBAS+1)/2
!   R0:       radii of nucleus
!   ZA:       Atomic nuclear charges
!   Coord:    Atomic Cartesian coordinates
!   FNorm:    common normalization factors (alpha independent part)
!   EXX:      primitive exponents
!   EXXP:     sqrt(EXX^LQ)
!   MaxXYZ:   (MaxTyp+1)*(MaxTyp+2)/2
!   MaxSum:   (MaxTyp+1)*(MaxTyp+2)*(MaxTyp+3)/6
!   facdat:   it saves (N)!
!   IRMS:     Index of the nucleus
!             IRMS > 0: dV/dR0 or dW/dR0 will be calculated
!             IRMS < 0: dV/dG or dW/dG will be calculated
!   e0:       a small value which is close to zero.
!   NOff:     offset of each primitive shell
!   MapA:     map relationship of each shell
!   LQnm:     L-quantum number of each shell
!   LPat:     (x,y,z) patterns of Cartesian functions
!   NEle:     the first dimension of dV and dVc, which can be 1 (for V and dV/dR0)
!             or 6 (for dV/dQ)
!
! Output:
!   DMat:     NAI matrix
!
! Scratch:
!   dV:       NAI<mu(a)|V(c)|nu(b)> for given (L,M,N)
!   dVc:      saves dV for S, P, D, F, ... shells
!   Ra,Rb,Rc: Positions of the centers (atoms) a, b, and c
!   aa,bb:    exponents of primitive Gaussian functions <mu| and |nu>
!   gg:       Gamma = aa + bb (gg(1)), and other aa,bb related values
!   Epsp:     powers of 1/(4*Gamma)
!             max index = 0,...,LQSum for the whole molecule system
!   AB:       Distance Rab
!   pab:      powers of PA and PB: 0,...,Lmu or Lnu
!   PC:       Distance PC
!   PCp:      powers of PC
!             max index = 0,...,LQSum for the whole molecule system
!   Airu,Ajsv,Aktw:
!             Coefficients matrix A
!   Gkp:      G(kappa,T)=F(kappa,T)/(1+Gamma*r0^2)^(kappa+1/2)
!             max index = 0,...,LQSum (or LQSum+1) for the whole molecule system
!   prefw:    prefactors of the terms in pVp.
!   faux:     an auxiliary array
!   tmp:      tmp. values
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NAImain(MaxXYZ,MaxSum,LpVp,ISO,NTT,NAtoms,NSHELL,IRMS,NEle,R0,ZA,Coord,Ra,Rb,Rc,FNorm,e0,prefw,dV,dVc,aa,bb,  &
  EXX,EXXP,gg,Epsp,AB,pab,PC,PCp,facdat,Airu,Ajsv,Aktw,Gkp,faux,tmp,NOff,MapA,LQnm,LPat,DMat)
  use Constant, only : Zero, One
  implicit real(kind=8) (a-h,o-z)
  parameter(Lgg = 18, Ltmp= 47, Lprf= 12)
  dimension R0(*),ZA(*),Coord(3,*),Ra(3),Rb(3),Rc(3),FNorm(*),prefw(Lprf),dV(NEle,MaxXYZ,MaxXYZ),dVc(NEle,MaxSum,MaxSum),        &
    EXX(*),EXXP(*),gg(Lgg),Epsp(*),AB(3),pab(*),PC(3),PCp(*),facdat(*),Airu(*),Ajsv(*),Aktw(*),Gkp(*),faux(*),tmp(Ltmp),         &
    NOff(*),MapA(*),LQnm(*),LPat(*),DMat(NTT,NEle)

  if(IRMS == 0)then
    IAtom1 = 1
    IAtom2 = NAtoms
    NDerRQ = 0
  else if(IRMS > 0)then
    IAtom1 = IRMS
    IAtom2 = IRMS
    NDerRQ = 1
  else if(IRMS < 0)then
    IAtom1 = -IRMS
    IAtom2 = -IRMS
    NDerRQ = 2
  end if

!-LOOP OVER ISHELL
  DO ISHELL=1,NSHELL
    aa     = EXX(ISHELL)
    Ra(1)  = Coord(1,MapA(ISHELL))
    Ra(2)  = Coord(2,MapA(ISHELL))
    Ra(3)  = Coord(3,MapA(ISHELL))
    ITYPE  = LQnm(ISHELL)
    IStart = (ITYPE+0)*(ITYPE+1)*(ITYPE+2)/6+1
    IEnd   = (ITYPE+1)*(ITYPE+2)*(ITYPE+3)/6
    IST    = NOff(ISHELL)

!---LOOP OVER JSHELL
    DO JSHELL=1,ISHELL
      bb     = EXX(JSHELL)
      Rb(1)  = Coord(1,MapA(JSHELL))
      Rb(2)  = Coord(2,MapA(JSHELL))
      Rb(3)  = Coord(3,MapA(JSHELL))
      JTYPE  = LQnm(JSHELL)
      JStart = (JTYPE+0)*(JTYPE+1)*(JTYPE+2)/6+1
      JEnd   = (JTYPE+1)*(JTYPE+2)*(JTYPE+3)/6
      IMJ    = ABS(ISHELL-JSHELL)
      JST    = NOff(JSHELL)

      LQSum=ITYPE+JTYPE+LpVp*2
      ! Rab
      AB(1)=Ra(1)-Rb(1)
      AB(2)=Ra(2)-Rb(2)
      AB(3)=Ra(3)-Rb(3)

      ! Gamma and other aa,bb related values
      call exponval(aa,bb,gg)
      ! powers of 1/(4*Gamma)
      call EpsPower(LQSum,gg(1),Epsp)
      ! powers of PA and PB
      call PowerPAB(.false.,LQSum,ITYPE+LpVp,JTYPE+LpVp,aa,bb,gg,Ra,Rb,tmp,pab)

!-----LOOP OVER ATOMS
      DO IAtom = IAtom1,IAtom2

        Rc(1)=Coord(1,IAtom)
        Rc(2)=Coord(2,IAtom)
        Rc(3)=Coord(3,IAtom)

        ! Integrals calculation
        dV = Zero
        call NAImunu(LQSum,MaxXYZ,NDerRQ,LpVp,ISO,e0,NEle,dV,gg,Epsp,facdat,prefw,aa,Ra,ITYPE, bb,Rb,JTYPE, LPat,  &
          R0(IAtom),Rc,AB,pab,PC,PCp,Airu,Ajsv,Aktw,Gkp,faux,tmp)

        ! dV is scaled by (-Zc)*(2pi/Gamma)*exp(...), Zc > 0, and saved to dVc.
        ! NOTE. EFG is not multiplied by Zc!
        if(NDerRQ == 2) then
          tmp(2) = -One
        else
          tmp(2) = -ZA(IAtom)
        end if
        Call ZExABG(aa,bb,gg(1),tmp(2),AB,tmp(1),tmp(3))
        do I=IStart,IEnd
          do J=JStart,JEnd
            dVc(:,J,I) = tmp(1) * dV(:,J-JStart+1,I-IStart+1)
          end do
        end do

        ! calculate DMat
        Do I = IStart,IEnd
          MU=IST+I-IStart+1
          MUS=(MU*(MU-1))/2
          tmp(1)=FNorm(I)*gg(15)
          If(I > 1)tmp(1)=tmp(1)*EXXP(ISHELL)
          JND=JEnd
          If(IMJ == 0) JND=I
          Do J =JStart,JND
            NU=JST+J-JStart+1
            MUNU=MUS+NU
            tmp(2)=FNorm(J)*gg(16)
            If(J > 1)tmp(2)=tmp(2)*EXXP(JSHELL)
            tmp(3)=tmp(1)*tmp(2)
            Do K = 1, NEle
              DMat(MUNU,K)=DMat(MUNU,K)+dVc(K,J,I)*tmp(3)
            end do
          end do
        end do

!-----IAtom
      end do

!---JSHELL
    end do

!-ISHELL
  end do

  return
end subroutine NAImain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the integral NAI<LQmu(a)|V(c)|LQnu(b)> between two primitive Cartesian Gaussian functions.
!
! Input:
!   LQmu,LQnu:quantum numbers of Cartesian Gaussian functions <mu| and |nu>, e.g. s=0,
!             p=1, d=2, ... NOTE: sp has been split into s and p
!   LQSum:    LQmu+LQnu; For pVp integrals, it is LQmu+LQnu+2.
!   LpVp:   = 0 for V
!           = 1 for pVp
!   ISO:    = 0 for V or pVp_sf
!           = 1 for pVp_x
!           = 2 for pVp_y
!           = 3 for pVp_z
!   MaxXYZ:   (MaxTyp+1)*(MaxTyp+2)/2
!   aa,bb:    exponents of primitive Gaussian functions <mu| and |nu>
!   gg:       Gamma = aa + bb (=gg(1)), and other aa,bb related values
!   Ra,Rb:    Cartesian coordinates of center <mu| and |nu>
!   AB:       Distance Rab
!   pab:      powers of PA and PB: 0,...,Lmu or Lnu
!   r0:       radius of nucleus C
!   Rc:       Cartesian coordinates of nucleus C
!   LPat:     (x,y,z) patterns of Cartesian functions
!   Epsp:     powers of 1/(4*Gamma)
!             max index = 0,...,LQSum for the whole molecule system
!   e0:       a small value which is close to zero.
!   facdat:   it saves (N)!
!   NDerRQ:   0 calculate V or pVp
!             1 calculate d(V)/dR0 or d(pVp)/dR0
!             2 calculate d(V)/dQ or d(pVp)/dQ
!   NEle:     the first dimension of dV, which can be 1 (for V and dV/dR0) or 6 (for dV/dQ)
!
! Output:
!   dV:       matrix NAI<mu(a)|V(c)|nu(b)> for given (L,M,N)
!
! Scratch:
!   pref:     prefactors of the terms in pVp.
!   PC:       Distance PC
!   PCp:      powers of PC
!             max index = 0,...,LQSum for the whole molecule system
!   Airu,Ajsv,Aktw:
!             Coefficients matrix A
!   Gkp:      G(kappa,T)=F(kappa,T)/(1+Gamma*r0^2)^(kappa+1/2)
!             max index = 0,...,LQSum (or LQSum+1) for the whole molecule system
!   faux:     an auxiliary array
!   tmp:      tmp. values
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine NAImunu(LQSum,MaxXYZ,NDerRQ,LpVp,ISO,e0,NEle,dV,gg,Epsp,facdat,pref,  aa,Ra,LQmu, bb,Rb,LQnu, LPat, r0,Rc,AB,pab,   &
  PC,PCp,Airu,Ajsv,Aktw,Gkp,faux,tmp)
  use Constant, only : One, Three
  implicit real(kind=8) (a-h,o-z)
  dimension dV(NEle,MaxXYZ,MaxXYZ),gg(*),Ra(3),Rb(3),Rc(3),AB(3),pab(*),PC(3),PCp(LQSum+1,3),Epsp(*),Airu(*),Ajsv(*),Aktw(*),  &
    facdat(*),Gkp(*),pref(6,2),faux(*),LPat(*),tmp(*)

! PC
  do i=1,3
!   Cartesian coordinates of point P
    tmp(1)=(aa*Ra(i)+bb*Rb(i))/gg(1)
    PC(i)=tmp(1)-Rc(i)
  end do

! Powers of PC
  do ix=1,3
    PCp(1,ix)=one
    do im=2,LQSum+1
      PCp(im,ix)=PCp(im-1,ix)*PC(ix)
    end do
  end do

! Calculates G-integrals
  tmp(1)=PC(1)*PC(1)+PC(2)*PC(2)+PC(3)*PC(3)
! Gamma*PC^2/(1+Gamma*r0^2)
  tmp(2)=tmp(1)*gg(1)/(One+gg(1)*r0*r0)
! G(kappa+3,T) and G(kappa+2,T) may be needed
  LQSum1=LQSum+NDerRQ+1
  call Fintgl(tmp(2),gg(1),r0,Gkp,LQSum1,faux,tmp(3))

! do loop in cartesian shells of mu and nu
  NKK=(LQmu+1)*(LQmu+2)/2
  NLL=(LQnu+1)*(LQnu+2)/2

  do KK=1,NKK
    call LPatt(LQmu, KK, LPat, Lmu, Mmu, Nmu)
    if(LpVp /= 0)then
!     prefactors of the terms in pVp
      tmp(1)=aa+aa
      pref(1,1)=tmp(1)
      if(Lmu > 0)pref(2,1)=-dble(Lmu)
      pref(3,1)=tmp(1)
      if(Mmu > 0)pref(4,1)=-dble(Mmu)
      pref(5,1)=tmp(1)
      if(Nmu > 0)pref(6,1)=-dble(Nmu)
    end if

    do LL=1,NLL
      call LPatt(LQnu, LL, LPat, Lnu, Mnu, Nnu)
      if(LpVp /= 0)then
!       prefactors of the terms in pVp
        tmp(1)=bb+bb
        pref(1,2)=tmp(1)
        if(Lnu > 0)pref(2,2)=-dble(Lnu)
        pref(3,2)=tmp(1)
        if(Mnu > 0)pref(4,2)=-dble(Mnu)
        pref(5,2)=tmp(1)
        if(Nnu > 0)pref(6,2)=-dble(Nnu)
      end if

      if(LpVp == 0)then
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       for V
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NAI<Lmu,Mmu,Nmu(a)|V(c)|Lnu,Mnu,Nnu(b)>
        call NAIlmn(LQSum,NDerRQ,NEle,e0,One,dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,Nmu,  bb,Lnu,Mnu,Nnu,  &
          pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
      else
        if(ISO == 0) then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         for pVp_sf
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         1) NAI<Lmu+1,Mmu,Nmu(a)|V(c)|Lnu+1,Mnu,Nnu(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(1,1)*pref(1,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu+1),Mmu,Nmu,  bb,(Lnu+1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         2) NAI<Lmu,Mmu+1,Nmu(a)|V(c)|Lnu,Mnu+1,Nnu(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(3,1)*pref(3,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu+1),Nmu,  bb,Lnu,(Mnu+1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         3) NAI<Lmu,Mmu,Nmu+1(a)|V(c)|Lnu,Mnu,Nnu+1(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(5,1)*pref(5,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu+1),  bb,Lnu,Mnu,(Nnu+1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         4) NAI<Lmu-1,Mmu,Nmu(a)|V(c)|Lnu+1,Mnu,Nnu(b)>
          if(Lmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(2,1)*pref(1,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu-1),Mmu,Nmu,  bb,(Lnu+1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         5) NAI<Lmu+1,Mmu,Nmu(a)|V(c)|Lnu-1,Mnu,Nnu(b)>
          if(Lnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(1,1)*pref(2,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu+1),Mmu,Nmu,  bb,(Lnu-1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         6) NAI<Lmu,Mmu-1,Nmu(a)|V(c)|Lnu,Mnu+1,Nnu(b)>
          if(Mmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(4,1)*pref(3,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu-1),Nmu,  bb,Lnu,(Mnu+1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         7) NAI<Lmu,Mmu+1,Nmu(a)|V(c)|Lnu,Mnu-1,Nnu(b)>
          if(Mnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(3,1)*pref(4,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu+1),Nmu,  bb,Lnu,(Mnu-1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         8) NAI<Lmu,Mmu,Nmu-1(a)|V(c)|Lnu,Mnu,Nnu+1(b)>
          if(Nmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(6,1)*pref(5,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu-1),  bb,Lnu,Mnu,(Nnu+1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         9) NAI<Lmu,Mmu,Nmu+1(a)|V(c)|Lnu,Mnu,Nnu-1(b)>
          if(Nnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(5,1)*pref(6,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu+1),  bb,Lnu,Mnu,(Nnu-1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         10) NAI<Lmu-1,Mmu,Nmu(a)|V(c)|Lnu-1,Mnu,Nnu(b)>
          if(Lmu*Lnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(2,1)*pref(2,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu-1),Mmu,Nmu,  bb,(Lnu-1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         11) NAI<Lmu,Mmu-1,Nmu(a)|V(c)|Lnu,Mnu-1,Nnu(b)>
          if(Mmu*Mnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(4,1)*pref(4,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu-1),Nmu,  bb,Lnu,(Mnu-1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         12) NAI<Lmu,Mmu,Nmu-1(a)|V(c)|Lnu,Mnu,Nnu-1(b)>
          if(Nmu*Nnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,pref(6,1)*pref(6,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu-1),  bb,Lnu,Mnu,(Nnu-1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
        else if(ISO == 1)then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         for pVp_x
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         1) NAI<Lmu,Mmu+1,Nmu(a)|V(c)|Lnu,Mnu,Nnu+1(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(3,1)*pref(5,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu+1),Nmu,  bb,Lnu,Mnu,(Nnu+1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         2) NAI<Lmu,Mmu+1,Nmu(a)|V(c)|Lnu,Mnu,Nnu-1(b)>
          if(Nnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(3,1)*pref(6,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu+1),Nmu,  bb,Lnu,Mnu,(Nnu-1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         3) NAI<Lmu,Mmu-1,Nmu(a)|V(c)|Lnu,Mnu,Nnu+1(b)>
          if(Mmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(4,1)*pref(5,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu-1),Nmu,  bb,Lnu,Mnu,(Nnu+1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         4) NAI<Lmu,Mmu-1,Nmu(a)|V(c)|Lnu,Mnu,Nnu-1(b)>
          if(Mmu*Nnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(4,1)*pref(6,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu-1),Nmu,  bb,Lnu,Mnu,(Nnu-1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         5) NAI<Lmu,Mmu,Nmu+1(a)|V(c)|Lnu,Mnu+1,Nnu(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(5,1)*pref(3,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu+1),  bb,Lnu,(Mnu+1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         6) NAI<Lmu,Mmu,Nmu+1(a)|V(c)|Lnu,Mnu-1,Nnu(b)>
          if(Mnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(5,1)*pref(4,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu+1),  bb,Lnu,(Mnu-1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         7) NAI<Lmu,Mmu,Nmu-1(a)|V(c)|Lnu,Mnu+1,Nnu(b)>
          if(Nmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(6,1)*pref(3,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu-1),  bb,Lnu,(Mnu+1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         8) NAI<Lmu,Mmu,Nmu-1(a)|V(c)|Lnu,Mnu-1,Nnu(b)>
          if(Nmu*Mnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(6,1)*pref(4,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu-1),  bb,Lnu,(Mnu-1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
        else if(ISO == 2)then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         for pVp_y
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         1) NAI<Lmu,Mmu,Nmu+1(a)|V(c)|Lnu+1,Mnu,Nnu(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(5,1)*pref(1,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu+1),  bb,(Lnu+1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         2) NAI<Lmu,Mmu,Nmu+1(a)|V(c)|Lnu-1,Mnu,Nnu(b)>
          if(Lnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(5,1)*pref(2,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu+1),  bb,(Lnu-1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         3) NAI<Lmu,Mmu,Nmu-1(a)|V(c)|Lnu+1,Mnu,Nnu(b)>
          if(Nmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(6,1)*pref(1,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu-1),  bb,(Lnu+1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         4) NAI<Lmu,Mmu,Nmu-1(a)|V(c)|Lnu-1,Mnu,Nnu(b)>
          if(Nmu*Lnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(6,1)*pref(2,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,Mmu,(Nmu-1),  bb,(Lnu-1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         5) NAI<Lmu+1,Mmu,Nmu(a)|V(c)|Lnu,Mnu,Nnu+1(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(1,1)*pref(5,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu+1),Mmu,Nmu,  bb,Lnu,Mnu,(Nnu+1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         6) NAI<Lmu+1,Mmu,Nmu(a)|V(c)|Lnu,Mnu,Nnu-1(b)>
          if(Nnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(1,1)*pref(6,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu+1),Mmu,Nmu,  bb,Lnu,Mnu,(Nnu-1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         7) NAI<Lmu-1,Mmu,Nmu(a)|V(c)|Lnu,Mnu,Nnu+1(b)>
          if(Lmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(2,1)*pref(5,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu-1),Mmu,Nmu,  bb,Lnu,Mnu,(Nnu+1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         8) NAI<Lmu-1,Mmu,Nmu(a)|V(c)|Lnu,Mnu,Nnu-1(b)>
          if(Lmu*Nnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(2,1)*pref(6,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu-1),Mmu,Nmu,  bb,Lnu,Mnu,(Nnu-1),  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
        else if(ISO == 3)then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         for pVp_z
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         1) NAI<Lmu+1,Mmu,Nmu(a)|V(c)|Lnu,Mnu+1,Nnu(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(1,1)*pref(3,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu+1),Mmu,Nmu,  bb,Lnu,(Mnu+1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         2) NAI<Lmu+1,Mmu,Nmu(a)|V(c)|Lnu,Mnu-1,Nnu(b)>
          if(Mnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(1,1)*pref(4,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu+1),Mmu,Nmu,  bb,Lnu,(Mnu-1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         3) NAI<Lmu-1,Mmu,Nmu(a)|V(c)|Lnu,Mnu+1,Nnu(b)>
          if(Lmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(2,1)*pref(3,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu-1),Mmu,Nmu,  bb,Lnu,(Mnu+1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         4) NAI<Lmu-1,Mmu,Nmu(a)|V(c)|Lnu,Mnu-1,Nnu(b)>
          if(Lmu*Mnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0, pref(2,1)*pref(4,2),dV(1,LL,KK),facdat,gg,  aa,(Lmu-1),Mmu,Nmu,  bb,Lnu,(Mnu-1),Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         5) NAI<Lmu,Mmu+1,Nmu(a)|V(c)|Lnu+1,Mnu,Nnu(b)>
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(3,1)*pref(1,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu+1),Nmu,  bb,(Lnu+1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         6) NAI<Lmu,Mmu+1,Nmu(a)|V(c)|Lnu-1,Mnu,Nnu(b)>
          if(Lnu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(3,1)*pref(2,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu+1),Nmu,  bb,(Lnu-1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         7) NAI<Lmu,Mmu-1,Nmu(a)|V(c)|Lnu+1,Mnu,Nnu(b)>
          if(Mmu /= 0)      &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(4,1)*pref(1,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu-1),Nmu,  bb,(Lnu+1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
!         8) NAI<Lmu,Mmu-1,Nmu(a)|V(c)|Lnu-1,Mnu,Nnu(b)>
          if(Mmu*Lnu /= 0)  &
          call NAIlmn(LQSum,NDerRQ,NEle,e0,-pref(4,1)*pref(2,2),dV(1,LL,KK),facdat,gg,  aa,Lmu,(Mmu-1),Nmu,  bb,(Lnu-1),Mnu,Nnu,  &
            pab, r0,PC,PCp,Epsp,Airu,Ajsv,Aktw,Gkp,tmp)
         end if
      end if
      if(NDerRQ == 2) then
        tmp(1)=(dV(1,LL,KK)+dV(3,LL,KK)+dV(6,LL,KK))/Three
        dV(1,LL,KK)=dV(1,LL,KK)-tmp(1)
        dV(3,LL,KK)=dV(3,LL,KK)-tmp(1)
        dV(6,LL,KK)=dV(6,LL,KK)-tmp(1)
      end if
    end do  ! LL
  end do  ! KK

  return
end Subroutine NAImunu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the integral NAI<LQmu(a)|V(c)|LQnu(b)> between two primitive Cartesian Gaussian functions
!
! Input:
!   LQSum:    LQmu+LQnu; For pVp integrals, it is LQmu+LQnu+2.
!   (Lmu,Mmu,Nmu),(Lnu,Mnu,Nnu):
!             quantum numbers of Cartesian Gaussian functions <mu| and |nu>
!   aa,bb:    exponents of primitive Gaussian functions <mu| and |nu>
!   gg:       Gamma = aa + bb (gg(1)), and other aa,bb related values
!   r0:       radius of nucleus C
!   PC:       Distance PC
!   PCp:      powers of PC: 0,...,LQSum
!   pab:      powers of PA and PB: 0,...,Lmu or Lnu
!   Gkp:      G(kappa,T)=F(kappa,T)/(1+Gamma*r0^2)^(kappa+1/2); index 1 corresponds to 0
!   Eps:      powers of 1/(4*Gamma); index 1 corresponds to 0
!   Airu,Ajsv,Aktw:
!             Coefficients matrix A and its 1st & 2nd derivatives.
!             (1,*) ... Coefficients matrix A
!             (2,*) ... 4(I2R2U + 1/2) x A
!             (3,*) ... I2R2U x A'
!             (4,*) ... I2R2U x (I2R2U - 1) x A"
!   e0:       a small value which is close to zero.
!   pref:     prefactors of the terms in pVp.
!   facdat:   it saves (N)!
!   NDerRQ:   0 calculate V or pVp
!             1 calculate d(V)/dR0 or d(pVp)/dR0
!             2 calculate d(V)/dQ or d(pVp)/dQ
!
! Output:
!   dV:       matrix element NAI<mu(a)|V(c)|nu(b)> for given (L,M,N)
!
! Scratch:
!   x:        tmp. values
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine NAIlmn(LQSum,NDerRQ,NEle,e0,pref,dV,facdat,gg,  aa,Lmu,Mmu,Nmu,  bb,Lnu,Mnu,Nnu, pab, r0,PC,PCp,Eps,Airu,Ajsv,Aktw,  &
  Gkp,x)
  use Constant, only : Zero, One
  implicit real(kind=8) (a-h,o-z)
  dimension dV(NEle),gg(*),pab(LQSum+1,2,3),PC(3),PCp(LQSum+1,3),Eps(*),facdat(*),Gkp(*),x(*),  &
    Airu(LQSum+1,LQSum/2+1,LQSum/2+1,*),Ajsv(LQSum+1,LQSum/2+1,LQSum/2+1,*),Aktw(LQSum+1,LQSum/2+1,LQSum/2+1,*)
  integer r,s,t,u,v,w

  If(NDerRQ == 2) then
    Airu(:,:,:,1:4) = Zero
    Ajsv(:,:,:,1:4) = Zero
    Aktw(:,:,:,1:4) = Zero
  else
    Airu(:,:,:,1) = Zero
    Ajsv(:,:,:,1) = Zero
    Aktw(:,:,:,1) = Zero
  end if

! Airu;
! x(:)    1: (-1)^i;  2: f_i;  3: x(1)*x(2)*i!;  4: (-1)^u;  5: 1/r!  6: x(4)*(1/u!)*(1/(i-2r-2u)!)
  x(1)=-One
  LI = Lmu+Lnu
  do i = 0, LI
    x(1)=-x(1)
    Call Fkd0(x(2),e0,facdat,pab(1,1,1),Lmu,Lnu,i,LQSum,x(3))
    x(3)=x(1)*x(2)*facdat(i+1)

    LR = i/2
    do r = 0, LR
      x(4)=-One
      x(5)=x(3)/facdat(r+1)
      I2R=i-r-r

      LU = LR - r
      do u = 0, LU
        x(4)=-x(4)
        I2R2U=I2R-u-u
        x(6)=x(4)*x(5)*Eps(r+u+1)/(facdat(u+1)*facdat(I2R2U+1))
        Airu(i+1,r+1,u+1,1)=x(6)*PCp(I2R2U+1,1)
        If(NDerRQ == 2) then
          Airu(i+1,r+1,u+1,2)=Airu(i+1,r+1,u+1,1)*dble(I2R2U*4+2)
          if(I2R2U > 0) then
            Airu(i+1,r+1,u+1,3) = x(6)*PCp(I2R2U,1)*dble(I2R2U)
            if(I2R2U > 1) Airu(i+1,r+1,u+1,4) = x(6)*PCp(I2R2U-1,1)*dble(I2R2U*(I2R2U-1))
          end if
        end if
      end do
    end do
  end do

! Ajsv;
  x(1)=-One
  LJ = Mmu+Mnu
  do j = 0, LJ
    x(1)=-x(1)
    Call Fkd0(x(2),e0,facdat,pab(1,1,2),Mmu,Mnu,j,LQSum,x(3))
    x(3)=x(1)*x(2)*facdat(j+1)

    LS = j/2
    do s = 0, LS
      x(4)=-One
      x(5)=x(3)/facdat(s+1)
      J2S=j-s-s

      LV = LS - s
      do v = 0, LV
        x(4)=-x(4)
        J2S2V=J2S-v-v
        x(6)=x(4)*x(5)*Eps(s+v+1)/(facdat(v+1)*facdat(J2S2V+1))
        Ajsv(j+1,s+1,v+1,1)=x(6)*PCp(J2S2V+1,2)
        If(NDerRQ == 2) then
          Ajsv(j+1,s+1,v+1,2)=Ajsv(j+1,s+1,v+1,1)*dble(J2S2V*4+2)
          if(J2S2V > 0) then
            Ajsv(j+1,s+1,v+1,3) = x(6)*PCp(J2S2V,2)*dble(J2S2V)
            if(J2S2V > 1) Ajsv(j+1,s+1,v+1,4) = x(6)*PCp(J2S2V-1,2)*dble(J2S2V*(J2S2V-1))
          end if
        end if
      end do
    end do
  end do

! Aktw;
  x(1)=-One
  LK = Nmu+Nnu
  do k = 0, LK
    x(1)=-x(1)
    Call Fkd0(x(2),e0,facdat,pab(1,1,3),Nmu,Nnu,k,LQSum,x(3))
    x(3)=x(1)*x(2)*facdat(k+1)

    LT = k/2
    do t = 0, LT
      x(4)=-One
      x(5)=x(3)/facdat(t+1)
      K2T=k-t-t

      LW = LT - t
      do w = 0, LW
        x(4)=-x(4)
        K2T2W=K2T-w-w
        x(6)=x(4)*x(5)*Eps(t+w+1)/(facdat(w+1)*facdat(K2T2W+1))
        Aktw(k+1,t+1,w+1,1)=x(6)*PCp(K2T2W+1,3)
        If(NDerRQ == 2) then
          Aktw(k+1,t+1,w+1,2)=Aktw(k+1,t+1,w+1,1)*dble(K2T2W*4+2)
          if(K2T2W > 0) then
            Aktw(k+1,t+1,w+1,3) = x(6)*PCp(K2T2W,3)*dble(K2T2W)
            if(K2T2W > 1) Aktw(k+1,t+1,w+1,4) = x(6)*PCp(K2T2W-1,3)*dble(K2T2W*(K2T2W-1))
          end if
        end if
      end do
    end do
  end do

  If(NDerRQ == 2) then
!   2Gamma
    x(01) = gg(1) + gg(1)
!   4Gamma^2
    x(02) = x(01) * x(01)
    x(28:40) = Zero

!   integrals
    do i = 0, LI

     LR = i/2
     do r = 0, LR

      LU = LR - r
      do u = 0, LU
        x(03) = pref * Airu(i+1,r+1,u+1,1)   ! Ax
        x(04) = pref * Airu(i+1,r+1,u+1,2)   ! Ax*
        x(05) = pref * Airu(i+1,r+1,u+1,3)   ! Ax'
        x(06) = pref * Airu(i+1,r+1,u+1,4)   ! Ax"
        kx = i-r-r-u

        do j = 0, LJ

         LS = j/2
         do s = 0, LS

          LV = LS - s
          do v = 0, LV
            x(07) = x(03) * Ajsv(j+1,s+1,v+1,1)    ! Ax  Ay
            x(08) = x(03) * Ajsv(j+1,s+1,v+1,2)    ! Ax  Ay*
            x(09) = x(03) * Ajsv(j+1,s+1,v+1,3)    ! Ax  Ay'
            x(10) = x(03) * Ajsv(j+1,s+1,v+1,4)    ! Ax  Ay"
            x(11) = x(04) * Ajsv(j+1,s+1,v+1,1)    ! Ax* Ay
            x(12) = x(05) * Ajsv(j+1,s+1,v+1,1)    ! Ax' Ay
            x(13) = x(05) * Ajsv(j+1,s+1,v+1,3)    ! Ax' Ay'
            x(14) = x(06) * Ajsv(j+1,s+1,v+1,1)    ! Ax" Ay
            ky = j-s-s-v

            do k = 0, LK

             LT = k/2
             do t = 0, LT

              LW = LT - t
              do w = 0, LW
                x(15) = x(07) * Aktw(k+1,t+1,w+1,1)    ! Ax  Ay  Az
                x(16) = x(07) * Aktw(k+1,t+1,w+1,2)    ! Ax  Ay  Az*
                x(17) = x(07) * Aktw(k+1,t+1,w+1,3)    ! Ax  Ay  Az'
                x(18) = x(07) * Aktw(k+1,t+1,w+1,4)    ! Ax  Ay  Az"
                x(19) = x(08) * Aktw(k+1,t+1,w+1,1)    ! Ax  Ay* Az
                x(20) = x(09) * Aktw(k+1,t+1,w+1,1)    ! Ax  Ay' Az
                x(21) = x(09) * Aktw(k+1,t+1,w+1,3)    ! Ax  Ay' Az'
                x(22) = x(10) * Aktw(k+1,t+1,w+1,1)    ! Ax  Ay" Az
                x(23) = x(11) * Aktw(k+1,t+1,w+1,1)    ! Ax* Ay  Az
                x(24) = x(12) * Aktw(k+1,t+1,w+1,1)    ! Ax' Ay  Az
                x(25) = x(12) * Aktw(k+1,t+1,w+1,3)    ! Ax' Ay  Az'
                x(26) = x(13) * Aktw(k+1,t+1,w+1,1)    ! Ax' Ay' Az
                x(27) = x(14) * Aktw(k+1,t+1,w+1,1)    ! Ax" Ay  Az
                kz = k-t-t-w

                kp = kx+ky+kz
!               combine the terms together to get the integrals
                x(28) = x(28) + x(27) * Gkp(kp+1)
                x(29) = x(29) + x(26) * Gkp(kp+1)
                x(30) = x(30) + x(22) * Gkp(kp+1)
                x(31) = x(31) + x(25) * Gkp(kp+1)
                x(32) = x(32) + x(21) * Gkp(kp+1)
                x(33) = x(33) + x(18) * Gkp(kp+1)
                x(34) = x(34) + x(23) * Gkp(kp+2)
                x(35) = x(35) + x(20) * Gkp(kp+2)
                x(36) = x(36) + x(24) * Gkp(kp+2)
                x(37) = x(37) + x(19) * Gkp(kp+2)
                x(38) = x(38) + x(17) * Gkp(kp+2)
                x(39) = x(39) + x(16) * Gkp(kp+2)
                x(40) = x(40) + x(15) * Gkp(kp+3)

              end do
             end do
            end do

          end do
         end do
        end do

      end do
     end do
    end do

    ! XX
    dV(1) = dV(1) + x(28) - gg(1)*x(34) + x(02)*PC(1)*PC(1)*x(40)
    ! XY
    dV(2) = dV(2) + x(29) - x(01)*PC(1)*x(35) - x(01)*PC(2)*x(36) + x(02)*PC(1)*PC(2)*x(40)
    ! YY
    dV(3) = dV(3) + x(30) - gg(1)*x(37) + x(02)*PC(2)*PC(2)*x(40)
    ! XZ
    dV(4) = dV(4) + x(31) - x(01)*PC(1)*x(38) - x(01)*PC(3)*x(36) + x(02)*PC(1)*PC(3)*x(40)
    ! YZ
    dV(5) = dV(5) + x(32) - x(01)*PC(2)*x(38) - x(01)*PC(3)*x(35) + x(02)*PC(2)*PC(3)*x(40)
    ! ZZ
    dV(6) = dV(6) + x(33) - gg(1)*x(39) + x(02)*PC(3)*PC(3)*x(40)

  else If(NDerRQ == 1) then
!   x(01): gg*R0/(1+gg*R0*R0)
    x(01) = gg(1)*R0
    x(01) = x(01)/(One+x(01)*R0)
!   x(02): 2gg * Rpc^2
    x(02) = PC(1)*PC(1) + PC(2)*PC(2) + PC(3)*PC(3)
    x(02) = gg(1) * x(02)
    x(02) = x(02) + x(02)
    x(06:07) = Zero

!   integrals
    do i = 0, LI

     LR = i/2
     do r = 0, LR

      LU = LR - r
      do u = 0, LU
        x(03) = pref * Airu(i+1,r+1,u+1,1)   ! Ax
        kx = i-r-r-u

        do j = 0, LJ

         LS = j/2
         do s = 0, LS

          LV = LS - s
          do v = 0, LV
            x(04) = x(03) * Ajsv(j+1,s+1,v+1,1)    ! Ax  Ay
            ky = j-s-s-v

            do k = 0, LK

             LT = k/2
             do t = 0, LT

              LW = LT - t
              do w = 0, LW
                x(05) = x(04) * Aktw(k+1,t+1,w+1,1)    ! Ax  Ay  Az
                kz = k-t-t-w

                kp = kx+ky+kz
!               combine the terms together to get the integrals
                x(06) = x(06) + x(05) * Gkp(kp+2)
                x(07) = x(07) + dble(kp+kp+1) * x(05) * Gkp(kp+1)

              end do
             end do
            end do

          end do
         end do
        end do

      end do
     end do
    end do

    dV(1) = dV(1) + x(01) * (x(02) * x(06) - x(07))

  else If(NDerRQ == 0) then
    x(06) = Zero

!   integrals
    do i = 0, LI

     LR = i/2
     do r = 0, LR

      LU = LR - r
      do u = 0, LU
        x(03) = pref * Airu(i+1,r+1,u+1,1)   ! Ax
        kx = i-r-r-u

        do j = 0, LJ

         LS = j/2
         do s = 0, LS

          LV = LS - s
          do v = 0, LV
            x(04) = x(03) * Ajsv(j+1,s+1,v+1,1)    ! Ax  Ay
            ky = j-s-s-v

            do k = 0, LK

             LT = k/2
             do t = 0, LT

              LW = LT - t
              do w = 0, LW
                x(05) = x(04) * Aktw(k+1,t+1,w+1,1)    ! Ax  Ay  Az
                kz = k-t-t-w

                kp = kx+ky+kz
!               combine the terms together to get the integrals
                x(06) = x(06) + x(05) * Gkp(kp+1)

              end do
             end do
            end do

          end do
         end do
        end do

      end do
     end do
    end do

    dV(1) = dV(1) + x(06)

  end if

  Return
End Subroutine NAIlmn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine calculates the 0th, 1st, and 2nd derivatives of S and T matrices in primitive Cartesian Gaussian functions.
! The basic method was described in
! H. Taketa, S. Huzinaga, and K. O-ohata, J. Phys. Soc. Jap. 21(11), 2313, 1966.
!
! Input:
!   LST            : 0 for S and 1 for T
!   NDer           : 0/1/2, the order of derivative
!   Ifgrd          : (.true./.false.) whether calculate gradients for NDer = 1. Not used in this program.
!   Den            : square density matrix to calculate egrad (for NDer = 1 only if Ifgrd = .true.). Not used in this program.
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   natom          : number of atoms
!   xyzmol         : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   nbasis         : number of primitive Cartesian functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!
! Output:
!   DMat(NTT,NMat) : NMat l.t. matrices in primitive Cartesian basis functions, where NMat = (NDer+1)*(NDer+2)/2
!   egrad          : gradients (for NDer = 1 only if Ifgrd = .true.). Not used in this program.
!
! Scratch:
!   V(Mem)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine pdst_drv(iout,LST,NDer,Ifgrd,mpatt,natom,xyzmol,maxtyp,nshell,nbasis,mxla,nshlla,gtoexp,DMat,Den,egrad,V,Mem)
  use Constant, only : maxlq, Zero
  implicit real(kind=8) (a-h,o-z)
  dimension :: xyzmol(*),mxla(natom),nshlla(0:maxlq,natom),gtoexp(*), DMat(*), Den(*), egrad(*), V(Mem)
  logical :: Ifgrd
  allocatable :: LPat(:),K(:)

  ! check
  if(LST /= 0 .and. LST /= 1) then
    write(iout,"(/,' Error in sub. pdst_drv: LST must be 0 or 1.')")
    call estop(1)
  end if

  if(NDer /= 0) then
    write(iout,"(/,' Error in sub. pdst_drv: NDer is out of range.')")
    call estop(1)
  end if

  ! ==================
  ! | initialization |
  ! ==================
  NMat=(NDer+1)*(NDer+2)/2
  MaxI1=MaxTyp+abs(LST)
  MaxI2=MaxI1+abs(LST)
  MaxSum=(MaxTyp+1)*(MaxTyp+2)*(MaxTyp+3)/6
  MaxXYZ=(MaxTyp+1)*(MaxTyp+2)/2
  MaxXYZ2=MaxXYZ*MaxXYZ

  NTT=NBasis*(NBasis+1)/2
  !NSS=NBasis*NBasis
  LGex= 11
  Ltmp= 12

  INrm = 1
  IE34 = INrm + MaxSum
  IRa  = IE34 + NSHELL
  IRb  = IRa  + 3
  Iab  = IRb  + 3
  Iexp = Iab  + 3
  Ifdt = Iexp + 1
  I2If = Ifdt + MaxI2+1
  IIoG = I2If + MaxI1+1
  IClv = IIoG + MaxI1+1
  Ipab = IClv + (MaxI2+1)*(MaxI2+1)
  Iprf = Ipab + (MaxI2+1)*(MaxI2+1)*6
  IGex = Iprf + 5
  Iddc = IGex + LGex
  if(NDer==0)then
    Itmp = Iddc
  else if(NDer==1)then
    Itmp = Iddc + 3
  else if(NDer==2)then
    Itmp = Iddc + 6
  end if
  IST  = Itmp + Ltmp
  ISTc = IST  + MaxXYZ2*NMat
  LenV = ISTc + MaxSum*MaxSum*NMat - 1
  if(LenV > Mem) then
    write(iout,"(/,' Error in sub. pdst_drv: LenV > Mem.')")
    call estop(1)
  end if

  KOff = 1
  KMap = KOff + NSHELL
  KLqn = KMap + NSHELL
  LenK = KLqn + NSHELL

  allocate(LPat(MaxSum*3),K(LenK))
  V(1:LenV) = Zero
  K = 0

  ! N!
  call FacCal(MaxI2,V(Ifdt))

  ! Generate (x,y,z) patterns in MOLDEN or Gaussian mode
  call GnPatt(MaxTyp,mpatt,LPat)

  ! common normalization factors
  call NormFac(MaxTyp,LPat,V(INrm))

  ! primitive exponents EXX and common factors EXXP=EXX^0.75
  call ExpP34(NSHELL,gtoexp,V(IE34))

  ! setup indices of primitive shells
  call SetShl(natom,nshell,mxla,nshlla,K(KOff),K(KMap),K(KLqn))

  DMat(1:NTT*NMat) = Zero

  call DDSTmain(MaxTyp,MaxXYZ,MaxSum,NSHELL,NTT,NDer,NMat,LST,xyzmol,K(KOff),K(KMap),K(KLqn),LPat,  &
    V(Iprf),V(INrm),V(Ifdt),gtoexp,V(IE34),V(Iexp),V(IRa),V(IRb),V(Iab),V(Ipab),V(I2If),V(IIoG),V(IClv),  &
    V(IGex),V(Iddc),V(Itmp),V(IST),V(ISTc),DMat)

  deallocate(LPat,K)

  Return
End subroutine pdst_drv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main subroutine to calculate the 0th, 1st, and 2nd derivatives of S and T matrices in primitive Cartesian Gaussian functions.
!
! Input:
!   MaxTyp    The highest angular quantum number present in the basis.
!   MaxXYZ:   (MaxTyp+1)*(MaxTyp+2)/2
!   MaxSum:   (MaxTyp+1)*(MaxTyp+2)*(MaxTyp+3)/6
!   NSHELL:   Number of primitive shells
!   NTT:      NTT=NBAS*(NBAS+1)/2. NBAS is the number of primitive Cartesian Gaussian functions.
!   NDer:     0/1/2, the order of derivative
!   NMat:     =(NDer+1)*(NDer+2)/2; the number of (derivative) matrices
!   LST:      0 for S and 1 for T
!   Coord:    Atomic Cartesian coordinates
!   NOff:     offset of each primitive shell
!   MapA:     map relationship of each shell
!   LQnm:     L-quantum number of each shell
!   LPat:     (x,y,z) patterns of Cartesian functions
!   FNorm:    common normalization factors (alpha independent part)
!   facdat:   it saves (N)!
!   EXX:      primitive exponents
!   EXX34:    EXX^(3/4)
!
! Output:
!   SMat:     NMat NDer-th derivative matrices
!
! Scratch:
!   pref:     prefactors of T integrals
!   fexp:     it saves (pi/gamma)^3/2 * exp[(-alpha*beta/gamma)*Rab^2]
!   Ra,Rb:    Positions of the centers (atoms) a and b
!   AB:       Ra - Rb
!   pab:      powers of PA and PB: 0,...,Lmu or Lnu + 2*abs(LST)
!   f2If:     it saves (2i-1)!!
!   fIoG:     it saves (2i-1)!!/(2gamma)^i
!   fClv:     it saves c(l,v)=l!/[v!(l-v)!]
!   Gexp:     alpha, beta, and gamma related values
!   ddc:      1st and second derivatives of fexp.
!   tmp:      tmp. values
!   S:        integrals for a given (L,M,N)
!   Sc :      saves integrals for S, P, D, F, ..., shells.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DDSTmain(MaxTyp,MaxXYZ,MaxSum,NSHELL,NTT,NDer,NMat,LST,Coord,NOff,MapA,LQnm,LPat,pref,FNorm,facdat,EXX,EXX34,fexp, &
  RA,RB,AB,pab,f2If,fIoG,fClv,Gexp,ddc,tmp,S,Sc,SMat)
  use Constant, only : pi, Zero, Two
  implicit real(kind=8) (a-h,o-z)
  dimension Coord(3,*),NOff(*),MapA(*),LQnm(*),LPat(*),pref(*),FNorm(*),facdat(*),EXX(*),EXX34(*),RA(3),RB(3),AB(3),pab(*),   &
    f2If(*),fIoG(*),fClv(*),Gexp(*),ddc(*),tmp(*),S(NMat,MaxXYZ,MaxXYZ),Sc(NMat,MaxSum,MaxSum),SMat(NTT,*)

  MaxI1=MaxTyp+abs(LST)
  MaxI2=MaxI1+abs(LST)

  ! (2n-1)!!, n=0,1,...,MaxI1
  do i=1,MaxI1+1
    f2If(i)=Tlm1ff(i)
  end do

  ! c(l,v)
  Call CombiFac(MaxI2,fClv,facdat)

  ! LOOP OVER ISHELL.
  DO ISHELL=1,NSHELL
    ! IA center to which ISHELL is attached
    IA=MapA(ISHELL)
    Ra(1) = Coord(1,IA)
    Ra(2) = Coord(2,IA)
    Ra(3) = Coord(3,IA)
    ITYPE = LQnm(ISHELL)
    IST = NOff(ISHELL)
    IStart=(ITYPE+0)*(ITYPE+1)*(ITYPE+2)/6+1
    IEnd=(ITYPE+1)*(ITYPE+2)*(ITYPE+3)/6

    Gexp(1)=EXX(ISHELL)        ! (1) Alpha
    Gexp(4)=EXX34(ISHELL)      ! (4) Alpha^3/4

    ! LOOP OVER JSHELL.
    DO JSHELL=1,ISHELL
      ! JA center to which JSHELL is attached
      JA=MapA(JSHELL)
      ! If IA=JA, the derivatives are zero, i.e., dS/dIA + dS/dJA = 0 and so on.
      ! In this subroutine only the first part is calculated.
      if(IA == JA .and. NDer > 0)cycle

      Rb(1) = Coord(1,JA)
      Rb(2) = Coord(2,JA)
      Rb(3) = Coord(3,JA)
      JTYPE = LQnm(JSHELL)
      JST = NOff(JSHELL)
      JStart=(JTYPE+0)*(JTYPE+1)*(JTYPE+2)/6+1
      JEnd=(JTYPE+1)*(JTYPE+2)*(JTYPE+3)/6
      IMJ=ABS(ISHELL-JSHELL)

      LQSum=ITYPE+JTYPE+abs(LST)*2
      ! Rab
      AB(1)=Ra(1)-Rb(1)
      AB(2)=Ra(2)-Rb(2)
      AB(3)=Ra(3)-Rb(3)

      Gexp(2)=EXX(JSHELL)        ! (2) Beta
      Gexp(5)=EXX34(JSHELL)      ! (5) Beta^3/4

      ! Gamma and other aa,bb related values
      Gexp(3)=Gexp(1)+Gexp(2)                                   ! (3) Gamma

      if(NDer > 0)then
        Gexp(6)=Gexp(1)/Gexp(3)                                 ! (6) Alpha/Gamma
        Gexp(7)=Gexp(2)/Gexp(3)                                 ! (7) Beta/Gamma
        Gexp(8)=Two*Gexp(1)*Gexp(2)/Gexp(3)                     ! (8) 2*Alpha*Beta/Gamma
      end if
      if(NDer > 1)then
        Gexp(9)=Gexp(6)*Gexp(6)                                 ! (9) (6)*(6)
        Gexp(10)=Gexp(7)*Gexp(7)                                ! (10) (7)*(7)
        Gexp(11)=Two*Gexp(6)*Gexp(7)                            ! (11) 2*(6)*(7)
      end if

      ! (2i-1)!!/(2gamma)^i
      do i=0,MaxI1
        fIoG(i+1)=f2If(i+1)*(Two*Gexp(3))**(-i)
      end do

      ! powers of PA and PB
      call PowerPAB(.false.,MaxI2,ITYPE,JTYPE+2*LST,Gexp(1),Gexp(2),Gexp(3),Ra,Rb,tmp,pab)

      ! (pi/gamma)^3/2 * exp[(-alpha*beta/gamma)*Rab^2]
      tmp(1) = sqrt(pi/Gexp(3))
      tmp(2) =-Gexp(1)*Gexp(2)/Gexp(3)
      tmp(3) = AB(1) * AB(1) + AB(2) * AB(2) + AB(3) * AB(3)
      fexp = tmp(1)*tmp(1)*tmp(1) * exp(tmp(2)*tmp(3))
      ! Derivatives of fexp
      if(NDer > 0)then
        ddc(1)=-AB(1)*Gexp(8)
        ddc(2)=-AB(2)*Gexp(8)
        ddc(3)=-AB(3)*Gexp(8)
      end if
      if(NDer > 1)then
        ddc(4)=ddc(1)*ddc(1)-Gexp(8)
        ddc(5)=ddc(2)*ddc(2)-Gexp(8)
        ddc(6)=ddc(3)*ddc(3)-Gexp(8)
      end if

      ! clean S
      S=Zero

      ! prefactors for T
      if(LST == 1)then
        pref(1)=Gexp(2)*dble(3+2*JTYPE)
        pref(2)=-Two*Gexp(2)*Gexp(2)
      end if

      ! Integrals calculation
      call DSTmunu(MaxTyp,MaxXYZ,LST,NDer,NMat,LPat,pref,Gexp,ITYPE,JTYPE,AB,fClv,pab,fIoG,ddc,tmp,S)

      ! S is scaled by fexp, and saved to Sc.
      do I=IStart,IEnd
        do J=JStart,JEnd
          do k=1,NMat
            Sc(K,J,I) = fexp * S(K,J-JStart+1,I-IStart+1)
          end do
        end do
      end do

      ! calculate SMat matrices
      Do I = IStart,IEnd
        MU=IST+I-IStart+1
        MUS=(MU*(MU-1))/2
        tmp(1)=FNorm(I)*Gexp(4)
        If(I > 1)tmp(1)=tmp(1)*Sqrt(Gexp(1)**ITYPE)
        JND=JEnd
        If(IMJ == 0) JND=I
        Do J =JStart,JND
          NU=JST+J-JStart+1
          MUNU=MUS+NU
          tmp(2)=FNorm(J)*Gexp(5)
          If(J > 1)tmp(2)=tmp(2)*Sqrt(Gexp(2)**JTYPE)
          tmp(3)=tmp(1)*tmp(2)
          do K=1,NMat
            SMat(MUNU,K)=SMat(MUNU,K)+Sc(K,J,I)*tmp(3)
          end do
        end do
      end do

    ! JSHELL
    end do

  ! ISHELL
  end do

  Return
End subroutine DDSTmain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the 0th, 1st, and 2nd derivative matrices of S or T between two primitive Cartesian Gaussian functions.
!
! Input:
!   LST:      0 for S and 1 for T
!   MaxTyp:   Max-L. see commonb.inc.
!   MaxXYZ:   (MaxTyp+1)*(MaxTyp+2)/2
!   LQmu,LQnu:quantum numbers of Cartesian Gaussian functions <mu|
!             and |nu>, e.g. s=0, p=1, d=2, ...
!             NOTE: sp has been split into s and p
!   NDer:     the order of derivative
!   NMat:     the number of matrices
!   LPat:     (lx,ly,lz) patterns
!   AB:       Ra - Rb
!
! Output:
!   S:        integrals for given (L,M,N)
!
! Scratch:
!   pref:     prefactors of T integrals
!   Gexp:     alpha, beta, and gamma related values
!   pab:      powers of PA and PB: 0,...,Lmu or Lnu + 2*abs(LST)
!   fIoG:     it saves (2i-1)!!/(2gamma)^i
!   fClv:     it saves c(l,v)=l!/[v!(l-v)!]
!   ddc:      1st and second derivatives of fexp.
!   tmp:      tmp. values
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DSTmunu(MaxTyp,MaxXYZ,LST,NDer,NMat,LPat,pref,Gexp, LQmu,LQnu,AB,  fClv,pab,fIoG,ddc,tmp,S)
  use Constant, only : One
  implicit real(kind=8) (a-h,o-z)
  dimension LPat(*), pref(*),AB(3),pab(*),fIoG(*),Gexp(*),ddc(*),fClv(MaxTyp+2*LST+1,*),tmp(*),S(NMat,MaxXYZ,MaxXYZ)

! do loop in cartesian shells of mu and nu
  NKK=(LQmu+1)*(LQmu+2)/2
  NLL=(LQnu+1)*(LQnu+2)/2
  do KK=1,NKK
    call LPatt(LQmu, KK, LPat, Lmu, Mmu, Nmu)

    do LL=1,NLL
      call LPatt(LQnu, LL, LPat, Lnu, Mnu, Nnu)

      if(LST == 0)then
!       for S
!       <Lmu,Mmu,Nmu(a)|Lnu,Mnu,Nnu(b)>
        call DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu, Mnu, Nnu, One,    Gexp,AB,fClv,pab,fIoG,ddc,tmp,S(1,LL,KK))
      else
!       for T
!       1) <Lmu,Mmu,Nmu(a)|Lnu,Mnu,Nnu(b)>
        call DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu, Mnu, Nnu, pref(1),Gexp,AB,fClv,pab,fIoG,ddc,tmp,S(1,LL,KK))
!       2) <Lmu,Mmu,Nmu(a)|Lnu+2,Mnu,Nnu(b)>
        call DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu+2,Mnu,Nnu, pref(2),Gexp,AB,fClv,pab,fIoG,ddc,tmp,S(1,LL,KK))
!       3) <Lmu,Mmu,Nmu(a)|Lnu,Mnu+2,Nnu(b)>
        call DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu,Mnu+2,Nnu, pref(2),Gexp,AB,fClv,pab,fIoG,ddc,tmp,S(1,LL,KK))
!       4) <Lmu,Mmu,Nmu(a)|Lnu,Mnu,Nnu+2(b)>
        call DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu,Mnu,Nnu+2, pref(2),Gexp,AB,fClv,pab,fIoG,ddc,tmp,S(1,LL,KK))
        if(LQnu > 1)then
!         prefactors of the terms in T
          pref(3)=dble(-Lnu*(Lnu-1)/2)
          pref(4)=dble(-Mnu*(Mnu-1)/2)
          pref(5)=dble(-Nnu*(Nnu-1)/2)
        else
          cycle
        end if
!       5) <Lmu,Mmu,Nmu(a)|Lnu-2,Mnu,Nnu(b)>
        if(Lnu > 1) call DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu-2,Mnu,Nnu, pref(3),Gexp,AB,fClv,pab,fIoG,ddc,tmp,S(1,LL,KK))
!       6) <Lmu,Mmu,Nmu(a)|Lnu,Mnu-2,Nnu(b)>
        if(Mnu > 1) call DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu,Mnu-2,Nnu, pref(4),Gexp,AB,fClv,pab,fIoG,ddc,tmp,S(1,LL,KK))
!       7) <Lmu,Mmu,Nmu(a)|Lnu,Mnu,Nnu-2(b)>
        if(Nnu > 1) call DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu,Mnu,Nnu-2, pref(5),Gexp,AB,fClv,pab,fIoG,ddc,tmp,S(1,LL,KK))
      end if
    end do
  end do

  return
end subroutine DSTmunu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the 0th, 1st, and 2nd derivative matrices of S or T between two primitive Cartesian Gaussian functions.
!
! Input:
!   MaxTyp:   Max-L. see commonb.inc.
!   LST:      0 for S and 1 for T
!   NDer:     the order of derivative
!   (Lmu,Mmu,Nmu),(Lnu,Mnu,Nnu):
!             quantum numbers of Cartesian Gaussian functions
!   AB:       Ra - Rb
!   pref:     prefactor of T integral
!   Gexp:     alpha, beta, and gamma related values
!   pab:      powers of PA and PB: 0,...,Lmu or Lnu + 2*abs(LST)
!   fIoG:     it saves (2i-1)!!/(2gamma)^i
!   fClv:     it saves c(l,v)=l!/[v!(l-v)!]
!   ddc:      1st and second derivatives of fexp.
!
! Output:
!   S:        integrals
!
! Scratch:
!   tmp:      tmp. values
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine DSTlmn(MaxTyp,LST,NDer, Lmu,Mmu,Nmu, Lnu,Mnu,Nnu, pref,Gexp,AB,fClv,pab,fIoG,ddc,tmp,S)
  implicit real(kind=8) (a-h,o-z)
  dimension Gexp(*),AB(3),ddc(*),fClv(MaxTyp+2*LST+1,*),pab(MaxTyp+2*LST+1,2,3),fIoG(*),tmp(*),S(*)

! 0th derivative integrals Sx, Sy, and Sz
  tmp(1) = S0Int(Lmu,Lnu,fClv(1,Lmu+1),fClv(1,Lnu+1), pab(1,1,1),pab(1,2,1),fIoG,tmp(4))
  tmp(2) = S0Int(Mmu,Mnu,fClv(1,Mmu+1),fClv(1,Mnu+1), pab(1,1,2),pab(1,2,2),fIoG,tmp(4))
  tmp(3) = S0Int(Nmu,Nnu,fClv(1,Nmu+1),fClv(1,Nnu+1), pab(1,1,3),pab(1,2,3),fIoG,tmp(4))

  if(NDer == 0)then
    S(1)=S(1)+pref*tmp(1)*tmp(2)*tmp(3)
    return
  else
    write(*,"(/,' Error in sub. DSTlmn: derivatives are not implemented in this program.')")
    call estop(1)
  end if

  Return
End Subroutine DSTlmn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate powers of PA and PB: from 0 to MaxLa or MaxLb
! if1c = .true. : one-center calculation; only MaxSize and pab are needed.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine PowerPAB(if1c,MaxSize,MaxLa,MaxLb,aa,bb,gg,Ra,Rb,tmp,pab)
  use Constant, only : Zero, One
  implicit real(kind=8) (a-h,o-z)
  logical if1c
  dimension Ra(3),Rb(3),tmp(*),pab(MaxSize+1,2,3)

  if(if1c) pab = Zero
  do ix=1,3                               ! x,y,z
    pab(1,1,ix) = One
    pab(1,2,ix) = One
    if(if1c) cycle

    tmp(1)=(aa*Ra(ix)+bb*Rb(ix))/gg       ! P
    tmp(2)=tmp(1)-Ra(ix)                  ! PA
    tmp(3)=tmp(1)-Rb(ix)                  ! PB
    do im=2,MaxLa+1
      pab(im,1,ix)=pab(im-1,1,ix)*tmp(2)
    end do
    do im=2,MaxLb+1
      pab(im,2,ix)=pab(im-1,2,ix)*tmp(3)
    end do
  end do

Return
End Subroutine PowerPAB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the 0th derivative of an S integral.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S0Int(La,Lb,clu,clv,pa,pb,fIoG,tmp)
  use Constant, only : Zero
  implicit real(kind=8) (a-h,o-z)
  dimension clu(*),clv(*),pa(*),pb(*),fIoG(*),tmp(3)
  integer u,v

  tmp(1)=Zero
  if(La > -1 .and. Lb > -1) then
    do i=0,(La+Lb)/2
      tmp(2)=Zero
      do u=0,La
        v=2*i-u
        if(v < 0) cycle
        if(v > Lb)cycle
        tmp(3) = clu(u+1)*clv(v+1)*pa(La-u+1)*pb(Lb-v+1)
        tmp(2) = tmp(2)+tmp(3)
      end do
      ! F_2i = tmp(2)
      tmp(1)=tmp(1)+tmp(2)*fIoG(i+1)
    end do
  end if
  S0Int = tmp(1)

  Return
End function S0Int

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the 1st derivative of an S integral.
!
!  NOTE.                        |-------|--|---------|
!                               |       |//|         |
!  In the block I (JGauss <=    |   0   | I|    0    |
!  IGauss), the taken values    |       |//|         |
!  must be multiplied by -1.    |-------|--|---------|
!                               |///I///|IA|////II///|
!                               |-------|--|---------|
!                               |       |//|         |
!                               |       |//|         |
!                               |   0   |II|    0    |
!                               |       |//|         |
!                               |-------|--|---------|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S1Int(La,Lb,clu,clv,pa,pb,fIoG,Gexp,tmp)
  use Constant, only : Zero
  implicit real*8(a-h,o-z)
  dimension clu(*),clv(*),pa(*),pb(*),fIoG(*),Gexp(*),tmp(*)
  integer u,v

  tmp(1)=Zero
  if(La > -1 .and. Lb > -1) then
    do i=0,(La+Lb)/2
      tmp(2)=Zero
      do u=0,La
        v=2*i-u
        if(v < 0) cycle
        if(v > Lb)cycle
        if(La-u > 0)then
          tmp(3) = clu(u+1)*clv(v+1)*pa(La-u)*pb(Lb-v+1)
          tmp(3) =-tmp(3)*dble(La-u)*Gexp(7)
          tmp(2) = tmp(2)+tmp(3)
        end if
        if(Lb-v > 0)then
          tmp(3) = clu(u+1)*clv(v+1)*pa(La-u+1)*pb(Lb-v)
          tmp(3) = tmp(3)*dble(Lb-v)*Gexp(6)
          tmp(2) = tmp(2)+tmp(3)
        end if
      end do
      ! F_2i' = tmp(2)
      tmp(1)=tmp(1)+tmp(2)*fIoG(i+1)
    end do
  end if
  S1Int = tmp(1)

  Return
End function S1Int

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generate (x,y,z) patterns for L <= MaxL
!
! Input
!   MaxL      the maximum angular momentum
!   IMod      0 for Gaussian mode (Z first and X last), > 0 for BDF mode (X first and Z last), and < 0 for MOLDEN mode
!
! Output
!   LPat      (x,y,z) patterns
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GnPatt(MaxL,IMod,LPat)
  implicit real(kind=8) (a-h,o-z)
  dimension LPat(3,*)

  Ixyz = 0
  if(IMod == 0) then
    ! Gaussian: same as MOLDEN for s,p,d,f-functions
    if(MaxL >= 0) then
      LPat(:, 1) = (/0,0,0/)
    end if
    if(MaxL >= 1) then
      LPat(:, 2) = (/1,0,0/)
      LPat(:, 3) = (/0,1,0/)
      LPat(:, 4) = (/0,0,1/)
    end if
    if(MaxL >= 2) then
      LPat(:, 5) = (/2,0,0/)
      LPat(:, 6) = (/0,2,0/)
      LPat(:, 7) = (/0,0,2/)
      LPat(:, 8) = (/1,1,0/)
      LPat(:, 9) = (/1,0,1/)
      LPat(:,10) = (/0,1,1/)
    end if
    if(MaxL >= 3) then
      LPat(:,11) = (/3,0,0/)
      LPat(:,12) = (/0,3,0/)
      LPat(:,13) = (/0,0,3/)
      LPat(:,14) = (/1,2,0/)
      LPat(:,15) = (/2,1,0/)
      LPat(:,16) = (/2,0,1/)
      LPat(:,17) = (/1,0,2/)
      LPat(:,18) = (/0,1,2/)
      LPat(:,19) = (/0,2,1/)
      LPat(:,20) = (/1,1,1/)
    end if
    Ixyz=20
    ! g, h, ...
    do L=4,MaxL
      do Lx=0,L
        Lyz=L-Lx
        do Lz=Lyz,0,-1
          Ly=Lyz-Lz
          Ixyz=Ixyz+1
          LPat(:,Ixyz) = (/Lx,Ly,Lz/)
        end do
      end do
    end do
  else if(IMod > 0) then
    ! BDF
    do L=0,MaxL
      do Lx=L,0,-1
        Lyz=L-Lx
        do Ly=Lyz,0,-1
          Lz=Lyz-Ly
          Ixyz=Ixyz+1
          LPat(:,Ixyz) = (/Lx,Ly,Lz/)
        end do
      end do
    end do
  else
    ! MOLDEN: same as Gaussian for s,p,d,f-functions
    if(MaxL >= 0) then
      LPat(:, 1) = (/0,0,0/)
    end if
    if(MaxL >= 1) then
      LPat(:, 2) = (/1,0,0/)
      LPat(:, 3) = (/0,1,0/)
      LPat(:, 4) = (/0,0,1/)
    end if
    if(MaxL >= 2) then
      LPat(:, 5) = (/2,0,0/)
      LPat(:, 6) = (/0,2,0/)
      LPat(:, 7) = (/0,0,2/)
      LPat(:, 8) = (/1,1,0/)
      LPat(:, 9) = (/1,0,1/)
      LPat(:,10) = (/0,1,1/)
    end if
    if(MaxL >= 3) then
      LPat(:,11) = (/3,0,0/)
      LPat(:,12) = (/0,3,0/)
      LPat(:,13) = (/0,0,3/)
      LPat(:,14) = (/1,2,0/)
      LPat(:,15) = (/2,1,0/)
      LPat(:,16) = (/2,0,1/)
      LPat(:,17) = (/1,0,2/)
      LPat(:,18) = (/0,1,2/)
      LPat(:,19) = (/0,2,1/)
      LPat(:,20) = (/1,1,1/)
    end if
    if(MaxL >= 4) then
      LPat(:,21) = (/4,0,0/)
      LPat(:,22) = (/0,4,0/)
      LPat(:,23) = (/0,0,4/)
      LPat(:,24) = (/3,1,0/)
      LPat(:,25) = (/3,0,1/)
      LPat(:,26) = (/1,3,0/)
      LPat(:,27) = (/0,3,1/)
      LPat(:,28) = (/1,0,3/)
      LPat(:,29) = (/0,1,3/)
      LPat(:,30) = (/2,2,0/)
      LPat(:,31) = (/2,0,2/)
      LPat(:,32) = (/0,2,2/)
      LPat(:,33) = (/2,1,1/)
      LPat(:,34) = (/1,2,1/)
      LPat(:,35) = (/1,1,2/)
    end if
    if(MaxL >= 5) then    ! same as Gaussian for H-functions (not officially supported by MOLDEN)
      Ixyz = 35
      do Lx=0,5
        Lyz=5-Lx
        do Lz=Lyz,0,-1
          Ly=Lyz-Lz
          Ixyz=Ixyz+1
          LPat(:,Ixyz) = (/Lx,Ly,Lz/)
        end do
      end do
    end if
    if(MaxL > 5) then
      write(*,"(/,' Error in sub. GnPatt: L > 5 has not been supported by MOLDEN now.')")
      call estop(1)
    end if
  end if

  Return
end subroutine GnPatt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the common normalization factors (alpha independent part)
! (2/Pi)^(3/4) * 2^L / [(2l-1)!! * (2m-1)!! * (2n-1)!!]^(1/2)
! where L = l + m + n
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NormFac(MaxL,LPat,Fn)
  use Constant, only : One, Two
  implicit real(kind=8) (a-h,o-z)
  dimension LPat(*), Fn(*)

 f1=sqrt(sqrt(Two/acos(-One)))
 f1=f1*f1*f1
 i=0
 do l=0,MaxL
   f2=Two**l
   do n=1,(l+1)*(l+2)/2
     i=i+1
     call LPatt(l, n, LPat, Lx, Ly, Lz)
     Fn(i)=Tlm1ff(Lx+1)*Tlm1ff(Ly+1)*Tlm1ff(Lz+1)
     Fn(i)=f1*f2/sqrt(Fn(i))
   end do
 end do

return
end subroutine NormFac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! input: Lquant, Lxyz, LPat
! output: (Lx, Ly, Lz), ie, (l, m, n)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LPatt(Lquant, Lxyz, LPat, Lx, Ly, Lz)
  implicit real(kind=8) (a-h,o-z)
  dimension LPat(3,*)
  dimension NOff(0:7)
  ! = (L'+1)*(L'+2)*(L'+3)/6 where L'=Lquant-1
  data NOff(:)/0,1,4,10,20,35,56,84/
  save NOff

  if(Lquant<0 .or. Lquant>7)then
    write(*,"(/,' LPatt is available only for Lq = 0 ~ 7!')")
    call estop(1)
  end if

  Ixyz = Lxyz + NOff(Lquant)

  Lx = LPat(1,Ixyz)
  Ly = LPat(2,Ixyz)
  Lz = LPat(3,Ixyz)

  Return
end subroutine LPatt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! EXX and EXXP=EXX^(3/4)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine ExpP34(N,EXX,EXXP)
  implicit real(kind=8) (a-h,o-z)
  dimension :: EXX(N),EXXP(N)

  do i = 1,N
    EXXP(i) = Sqrt(EXX(i))
    EXXP(i) = Sqrt(EXXP(i)*EXXP(i)*EXXP(i))
  end do

  Return
End Subroutine ExpP34

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1. set offset for each primitive shell, i.e. the total number of Cartesian functions before this shell.
! 2. set map relationship of each shell.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SetShl(natom,nshell,mxla,nshlla,NOff,MapA,LQnm)
  use Constant, only : maxlq
  implicit real(kind=8) (a-h,o-z)
  dimension mxla(natom),nshlla(0:maxlq,natom),NOff(nshell),MapA(nshell),LQnm(nshell)

  NOff(1)=0
  ishell=0
  do iatom=1,natom
    do li=0,mxla(iatom)
      ncar=(li+1)*(li+2)/2
      do ip=1,nshlla(li,iatom)
        ishell=ishell+1
        if(ishell < nshell) NOff(ishell+1)=NOff(ishell)+ncar
        MapA(ishell) = iatom
        LQnm(ishell) = li
      end do
    enddo
  enddo

return
end subroutine SetShl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! gtoex and gtoexpp=sqrt(gtoex^LQ)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine ExpPsq(natom,mxla,nshlla,gtoex,gtoexpp)
  use Constant, only : maxlq, One
  implicit real(kind=8) (a-h,o-z)
  dimension :: mxla(natom),nshlla(0:maxlq,natom),gtoex(*),gtoexpp(*)

  ishell = 0
  do iatom=1,natom
    do li=0,mxla(iatom)
      do ip = 1, nshlla(li,iatom)
        ishell = ishell + 1
        if (li == 0) then
          gtoexpp(ishell) = One
        else
          gtoexpp(ishell) = Sqrt(gtoex(ishell)**li)
        end if
      end do
    enddo
  enddo

  Return
End Subroutine ExpPsq

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate a,b related values:
! (1) =g,       (2) =2*a,     (3) =2*b,     (4) =2*g,     (5) =a/g,
! (6) =b/g,     (7) =(2)*(2), (8) =(3)*(3), (9) =(4)*(4), (10)=(2)*(5),
! (11)=(3)*(6), (12)=(5)*(5), (13)=(6)*(6), (14)=(5)*(6), (15)=a^(3/4),
! (16)=b^(3/4), (17)=(2)*(4), (18)=sqrt(a*b)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine exponval(a,b,g)
  implicit real(kind=8) (a-h,o-z)
  dimension g(18)

  g(1) = a + b
  g(2) = a + a
  g(3) = b + b
  g(4) = g(1) + g(1)
  g(5) = a / g(1)
  g(6) = b / g(1)
  g(7) = g(2) * g(2)
  g(8) = g(3) * g(3)
  g(9) = g(4) * g(4)
  g(10)= g(2) * g(5)
  g(11)= g(3) * g(6)
  g(12)= g(5) * g(5)
  g(13)= g(6) * g(6)
  g(14)= g(5) * g(6)
  g(15)= Sqrt(a)
  g(15)= Sqrt(g(15) * g(15) * g(15))
  g(16)= Sqrt(b)
  g(16)= Sqrt(g(16) * g(16) * g(16))
  g(17)= g(2) * g(4)
  g(18)= Sqrt(a * b)

  Return
end subroutine exponval

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the powers of 1/(4Gamma) in A(i,r,u), A(j,s,v), A(k,t,w)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine EpsPower(MaxL,Gamma,Epsp)
  use Constant, only : Quarter, One
  implicit real(kind=8) (a-h,o-z)
  dimension Epsp(*)

  Epsp(1)=One
  if(MaxL > 0)Epsp(2)=Quarter/Gamma
  do I=3,MaxL+1
    Epsp(I)=Epsp(I-1)*Epsp(2)
  end do

  Return
End Subroutine EpsPower

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the scaling factor Zc*(2pi/Gamma)*exp(-aa*bb/gg*Rab^2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine ZExABG(aa,bb,gg,Zc,AB,fe,tmp)
  use Constant, only : One, Two
  implicit real(kind=8) (a-h,o-z)
  dimension AB(3),tmp(*)

  tmp(1) = Two * acos(-One)
  tmp(2) = AB(1) * AB(1) + AB(2) * AB(2) + AB(3) * AB(3)
  tmp(3) = exp(- aa * bb * tmp(2) / gg)
  fe = Zc * tmp(1) * tmp(3) / gg

  Return
End Subroutine ZExABG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates F_nu(z) integrals with precision better than 1.d-12.
!
!     z ... argument
!    r0 ... nuclear RMS charge radii in a.u.
!   fnu ... array of F integrals from nu=0 to nu=24
!           Note, i = nu+1
!           This must be more than enough for calculations
!  faux ... an auxiliary array
! maxnu ... maximal value of (nu+1)
!   tmp ... scratch
!
!  For z < 100: uses a downward recursion starting with numax=maxnu+24 and approximating f(numax) by its value for z=0.
!  F(nu,0)=1/(2*nu+1)
!  F(nu,z)=(z*F(nu+1,z)+Exp(-z)/2)/(nu+1/2)
!
!  For z >=100: uses an upward recursion starting with F(0,z).
!  F(nu+1,z)=(nu+1/2)*F(nu,z)/z
!
!  Note: exponent is neglected here!
!
!  Uses an external function derf to calculate F(0,z).
!  F(0,z)=1/2 * Sqrt(Pi/z) * Erf(Sqrt(z))
!
!  Finally, F(nu,z)=F(nu,z)/(1+gamma*r0**2)**(nu+1/2)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Fintgl(z,gamma,r0,fnu,maxnu,faux,tmp)
  use Constant, only : One,Half
  implicit real(kind=8) (a-h,o-z)
  parameter(LimD=55)
  dimension :: fnu(*),faux(maxnu+LimD),tmp(7)

  tmp(4)=Sqrt(Atan(One))    ! Sqrt(Pi)/2
  tmp(5)=One/(One+gamma*r0*r0)
  tmp(2)=sqrt(tmp(5))

  if (abs(z) <= EpsFun(25)) then
!   z=0: F(nu,0)=1/(2*nu+1), here i=nu+1. Then F(nu,0)/(1+gamma*r0**2)**(nu+1/2) --> F(nu,0)
    fnu(1)=tmp(2)
    tmp(3)=tmp(2)
    do i=2,maxnu
     tmp(3)=tmp(3)*tmp(5)
     fnu(i)=tmp(3)/dble(i+i-1)
    end do
  else
    tmp(6)=sqrt(z)
    fnu(1)=tmp(4)/tmp(6)*derf(tmp(6))
    if (abs(z) < 1.d2) then
      tmp(7)=Half*exp(-z)
!     Starting downward recursion
      lim=maxnu+LimD
      faux(lim)=One/(dble(2*lim)+One)
      do i=1,lim-1
        faux(lim-i)=(z*faux(lim+1-i)+tmp(7))/(dble(lim-1-i)+Half)
      end do
!     Factor for renormalization
      tmp(3)=tmp(2)*fnu(1)/faux(1)
    else
!     Starting upward recursion
      faux(1)=fnu(1)
      tmp(1)=one/z
      do i=1,maxnu
        faux(i+1)=(dble(i-1)+Half)*faux(i)*tmp(1)
      end do
      tmp(3)=tmp(2)
    end if
!   Renormalization
    fnu(1)=fnu(1)*tmp(2)
    if (maxnu >= 2) then
     do i=2,maxnu
       tmp(3)=tmp(3)*tmp(5)
       fnu(i)=faux(i)*tmp(3)
     end do
    end if
  end if

  Return
End Subroutine Fintgl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate coefficients F_k(l1,l2,a,b)
!   e0:       a small value which is close to zero.
!   facdat:   it saves (N)!
!
! NOTE:
!   1) k, l1 and l2 start from 0.
!   2) Except for k=l1+l2, 0^0 is unity.
!   3) It is assumed that l1 >= l2.
!
! Working formula:
! Fk = SUM(j=0,l2) [l1! * l2! / ( (k-j)! (l1-k+j)! j! (l2-j)! )] * a^(l1-k+j) * b^(l2-j)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine Fkd0(Fk,e0,facdat,pab,l1,l2,k,LQSum,x)
  use Constant, only : Zero, One
  implicit real(kind=8) (a-h,o-z)
  dimension :: pab(LQSum+1,2),facdat(*),x(*)

! Special case k=l1+l2, Fk = 1
  if(k==l1+l2) then
    Fk=One
  else
!   First check if l1 >= l2. If not => swap a and b.
    if(l1 < l2) then
      l1x=l2
      l2x=l1
      ia=2
      ib=1
    else
      l1x=l1
      l2x=l2
      ia=1
      ib=2
    end if
!   Now starts summation
    x(1)=Zero

    x(2)=facdat(l1x+1)*facdat(l2x+1)
    do j=0,l2x
      ind1=l1x-k+j
      ind2=l2x-j
      ind3=k-j
!     Check for negative index
      idx=(isign(1,ind1)+1)*(isign(1,ind2)+1)*(isign(1,ind3)+1)
      if (idx.ne.0) then
        x(3)=x(2)/(facdat(ind3+1)*facdat(ind1+1)*facdat(ind2+1)*facdat(j+1))
        x(1)=x(1)+x(3)*pab(ind1+1,ia)*pab(ind2+1,ib)
      end if
    end do

    Fk=x(1)
!   if(abs(Fk) < e0)Fk=e0
  end if

  Return
End Subroutine Fkd0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine calculates DDD integrals in primitive Cartesian Gaussian functions.
! The basic method was described in
! H. Taketa, S. Huzinaga, and K. O-ohata, J. Phys. Soc. Jap. 21(11), 2313, 1966.
!
! Input:
!   natom          : number of atoms
!   nshell         : number of primitive shells
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   NTT            : nbasis*(nbasis+1)/2. nbasis: number of primitive Cartesian functions
!
! Output:
!   DMat(NTT) : L.T. matrices in primitive Cartesian basis functions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ddd_drv(natom,nshell,mxla,nshlla,gtoexp,NTT,DMat)
  use Constant, only : maxlq, Zero
  implicit real(kind=8) (a-h,o-z)
  dimension :: mxla(natom),nshlla(0:maxlq,natom),gtoexp(*), DMat(NTT)
  allocatable :: K(:)

  ! ==================
  ! | initialization |
  ! ==================

  KOff = 1
  KMap = KOff + NSHELL
  KLqn = KMap + NSHELL
  LenK = KLqn + NSHELL

  allocate(K(LenK))
  K = 0

  ! setup indices of primitive shells
  call SetShl(natom,nshell,mxla,nshlla,K(KOff),K(KMap),K(KLqn))

  DMat(:) = Zero

  call DDDmain(NSHELL,NTT,K(KOff),K(KMap),K(KLqn),gtoexp,DMat)

  deallocate(K)

  Return
End subroutine ddd_drv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main subroutine to calculate the S, dV/dF and dW/dF matrices in primitive Cartesian Gaussian functions.
!
! Input:
!   NSHELL    Number of primitive shells
!   NTT       NTT=NBAS*(NBAS+1)/2. NBAS is the number of primitive Cartesian Gaussian functions.
!   NOff      offset of each primitive shell
!   MapA      map relationship of each shell
!   LQnm      L-quantum number of each shell
!   EXX       primitive exponents
!
! Output:
!   SMat      DDD
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DDDmain(NSHELL,NTT,NOff,MapA,LQnm,EXX,SMat)
  implicit real(kind=8) (a-h,o-z)
  dimension NOff(*),MapA(*),LQnm(*),EXX(*),SMat(NTT)

  loop_ISHELL: DO ISHELL=1,NSHELL
    ! IA center to which ISHELL is attached
    IA=MapA(ISHELL)
    ITYPE = LQnm(ISHELL)
    IST = NOff(ISHELL)
    IStart=(ITYPE+0)*(ITYPE+1)*(ITYPE+2)/6+1
    IEnd=(ITYPE+1)*(ITYPE+2)*(ITYPE+3)/6

    loop_JSHELL: DO JSHELL=1,ISHELL
      ! JA center to which JSHELL is attached
      JA=MapA(JSHELL)
      if(IA == JA)cycle

      JTYPE = LQnm(JSHELL)
      JST = NOff(JSHELL)
      JStart=(JTYPE+0)*(JTYPE+1)*(JTYPE+2)/6+1
      JEnd=(JTYPE+1)*(JTYPE+2)*(JTYPE+3)/6

      ! calculate SMat matrices
      Do I = IStart,IEnd
        MU=IST+I-IStart+1
        MUS=(MU*(MU-1))/2
        JND=JEnd
        If(ISHELL-JSHELL == 0) JND=I
        Do J =JStart,JND
          NU=JST+J-JStart+1
          MUNU=MUS+NU
          SMat(MUNU)= EXX(ISHELL) / (EXX(ISHELL) + EXX(JSHELL))
        end do
      end do

    end do loop_JSHELL
  end do loop_ISHELL

  Return
End subroutine DDDmain

!--- END
