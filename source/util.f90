!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Print title
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine prt_title(iout)
 use Constant, only : Ver, Date
 implicit real(kind=8) (a-h,o-z)

 write(iout,"(1x,77('='),/, &
 ' ==  RelEPro: a program to calculate relativistic one-electron properties.  ==',/, &
 ' ==',73x,'==',/, &
 ' ==',44x,'Ver. ',a5,',  ',a12,4x,'==',/, &
 1x,77('='),/ )") Ver, Date

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! count the number of point charges (npchg = 0) or read point charges (npchg > 0)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine rd_pchar(iinp,iout,npchg,pchg,pxyz,ctmp)
 use Constant, only : au2ang
 implicit real(kind=8) (a-h,o-z)
 dimension          :: pchg(*), pxyz(3,*)
 character*200      :: ctmp
 namelist/PChar/unit

 rewind(iinp)

 if(npchg > 0) then
   unit = 0.0d0
   read(iinp,PChar)
   do i = 1, npchg
     read(iinp,*) pchg(i), pxyz(:,i)
   end do
   if(nint(unit) == 1) pxyz(:,1:npchg) = pxyz(:,1:npchg)/au2ang
   call prtpcoord(iout,npchg,pchg,pxyz)
 else if(npchg == 0) then
   read(iinp,PChar,end=1000,err=1000)
   do while(.true.)
     read(iinp,"(a200)",end=200,err=1200) ctmp
     if(len_trim(ctmp) ==0) then
       exit
     else if(len_trim(ctmp) > 6) then
       read(ctmp,*,err=1200) pchg(1), pxyz(:,1)
       npchg = npchg + 1
     else
       go to 1200
     end if
   end do
   200  continue
 end if

 1000  continue
 return
 1200  write(iout,"(/,' The format of point charge (',i3,') is wrong!')") npchg+1
 call estop(1)
end Subroutine rd_pchar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read nuclear quadrupole moments (Q in millibarn).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine rd_nqm(iinp,iout,natom,qdat,ctmp)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension          :: qdat(2,natom)
 character*200      :: ctmp
 namelist/nqmdat/dum

 qdat = Zero

 rewind(iinp)

 read(iinp,nqmdat,end=1000,err=1000)
 do while(.true.)
   read(iinp,"(a200)",end=1000,err=1200) ctmp
   if(len_trim(ctmp) ==0) exit

   read(ctmp,*,err=1200) ia, q
   if(ia < 1 .or. ia > natom) then
     goto 1100
   else
     qdat(1,ia) = One
     qdat(2,ia) = q
   end if
 end do

 1000  continue
 return
 1100  write(iout,"(/,' The atomic index ',i3,' of NQM data is wrong!')") ia
 call estop(1)
 1200  write(iout,"(/,' The format of NQM data is wrong!')")
 call estop(1)
end Subroutine rd_nqm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read input and open data file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine open_dfile(iinp,iout,ifham,ifpro,ifdbg,irel,NIOP,IOP,NJob,JobTyp,imod,fnam)
 implicit real(kind=8) (a-h,o-z)
 logical            :: ifham, ifpro, ifdbg
 dimension          :: IOP(NIOP), JobTyp(NJob)
 character*200      :: fnam
 character*61       :: string
 character*20       :: prop
 character*1        :: cha1(0:1)
 integer            :: popu
 namelist/Contrl/prop,ifmt,iham,iden,minz,nucl,loca,popu
 namelist/QCData/fnam

 IOP=0
 JobTyp=0

 cha1 = (/' ','x'/)
 prop = " "
 ifmt = 0
 iham = 0
 iden = 0
 minz = 0
 nucl = 0
 loca = 0
 popu = 0
 irel = 1  ! 1c / 2c

 ! the first rewind is not needed by ifort:
 ! warning (499): REWIND is a no-op on terminal device
 !rewind(iinp)

 read(iinp,Contrl,end=100,err=10)
 goto 100
 10  write(iout,"(/,' $Contrl is wrong!')")
 call estop(1)

 100  continue
 ifham = .false.  ! whether compute Hamiltonian
 ifpro = .false.  ! whether compute properties
 if(len_trim(prop) < 1) then
   ! default: do nothing
   !JobTyp(1) = 1
   !ifham = .true.
 else
   call charl2u(prop)

   if(index(prop,'+DIP') > 0) then
     JobTyp(1) = 1
     ifham = .true.
     ifpro = .true.
   end if
   if(index(prop,'+APT') > 0 .and. ifdbg) then
     JobTyp(2) = 1
     ifham = .true.
     ifpro = .true.
   end if
   if(index(prop,'+DDD') > 0 .and. ifdbg) then
     JobTyp(3) = 1
     ifham = .true.
     ifpro = .true.
   end if
   if(index(prop,'+ED') > 0) then
     JobTyp(11) = 1
     ifham = .true.
     ifpro = .true.
   end if
   if(index(prop,'+CD') > 0) then    ! Hamiltonian is not needed for n.r. CD
     JobTyp(12) = 1
     ifpro = .true.
   end if
   if(index(prop,'+EFG') > 0) then
     JobTyp(13) = 1
     ifham = .true.
     ifpro = .true.
   end if

   if(isum(NJob,JobTyp) < 1) then       ! default: DIP calculation
     JobTyp(1) = 1
     ifham = .true.
     ifpro = .true.
   end if
 end if

 call hostnm( string )
 write(iout,"(/,' One-electron properties to be computed on the computer ',a,/)") trim( string )
 write(iout,"(2x,'( ',a1,' ) electric dipole moment (DIP)')") cha1(JobTyp(1))
 if(ifdbg) then
   write(iout,"(2x,'( ',a1,' ) atomic polar tensor (APT), no CPHF procedure')") cha1(JobTyp(2))
   write(iout,"(2x,'( ',a1,' ) atomic DDD charges (DDD)')") cha1(JobTyp(3))
 end if
 write(iout,"(2x,'( ',a1,' ) effective contact density (ED)')") cha1(JobTyp(11))
 write(iout,"(2x,'( ',a1,' ) contact density with PC effects (CD)')") cha1(JobTyp(12))
 write(iout,"(2x,'( ',a1,' ) electric field gradient (EFG)')") cha1(JobTyp(13))
 write(iout,*)

 ! if properties are not computed, the Hamiltonian is not required.
 if(.not. ifpro) ifham = .false.
 ! molden, fchk, et. al.
 if(ifmt /= 1 .and. ifmt /= 2) ifmt = 0
 ! DLX, AX, NR, et. al.
 if(iham < 1 .or. iham > 6) iham = 4
 ! X2C Hamiltonian is not needed for n.r. calculations
 if(iham > 1) ifham = .true.
 ! APT is computed by NR formula
 if(JobTyp(2) == 1 .and. iham /= 1) then
   write(iout,"(/,' Error! APT has been programmed for the NR Hamiltonian only.',/)")
   call estop(1)
 end if
 ! DDD is computed by NR formula
 if(JobTyp(3) == 1 .and. iham /= 1) then
   write(iout,"(/,' Error! DDD has been programmed for the NR Hamiltonian only.',/)")
   call estop(1)
 end if
 ! ZA or min(ZA) for ED/CD/EFG
 if(minz == 0) minz = 1
 ! Finite nuclear charge density distribution
 if(nucl /= 1) nucl = 1
 ! atomic block of dNAI
 if(loca /= 2 .and. loca /= 3) loca = 1
 if(.NOT. ifdbg) loca = 1
 ! population
 if(ifham) then
   if(popu /= 1 .or. loca /= 1) popu = 0
 else
   popu = 0
 end if
 IOP(1) = ifmt
 IOP(2) = iham
 IOP(3) = minz
 IOP(4) = nucl
 IOP(5) = loca

 fnam=" "
 rewind(iinp)
 read(iinp,QCData,end=200,err=110)
 goto 200
 110  write(iout,"(/,' $QCData is wrong!')")
 call estop(1)

 200  continue
 lstr=nonspace(fnam)
 lend=LEN_TRIM(fnam)
 if(lend == 0)then                 ! use default file name
   if(IOP(1) == 1) then
     lstr=1
     lend=6
     fnam='MOLDEN'
   end if
 end if
 open(imod,file=fnam(lstr:lend),status='old',err=210)
 goto 300
 210  write(iout,"(' ### Wrong! The file ',a,' does not exist!')") trim(fnam)
 call estop(1)

 300  continue
 if(IOP(1) == 0) then
   call dfiletest(imod,ifmt,string)
   if(ifmt == 1 .or. ifmt == 2) IOP(1) = ifmt
 end if

 if(IOP(1) == 1) then
   write(iout,"(/,' The file ',a,' in the MOLDEN format has been found.',/)") trim(fnam)
   iden =-1
 else if(IOP(1) == 2) then
   write(iout,"(/,' The file ',a,' in the FCHK format has been found.',/)") trim(fnam)
   if(iden < 0) then
     iden =-1
   else if(iden /= 2) then
     iden = 1
   end if
   ! 1c or 2c?
   call reltest(imod,irel,string)
   if(irel == 2 .and. iden == 2) then
     write(iout,"(/,' Two-component post-HF calculation is not supported.',/)")
     call estop(1)
   end if
 else
   write(iout,"(/,' The file ',a,' has been found, but its format is not known.',/)") trim(fnam)
   call estop(1)
 end if

 if(iden == 2) popu = 0

 IOP(6) = popu
 IOP(7) = iden

 return
end subroutine open_dfile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! check the format of a data file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine dfiletest(imod,ifmt,string)
 implicit real(kind=8) (a-h,o-z)
 character*61 :: string

 ifmt = 0  ! 1 (molden) / 2 (fchk)
 rewind(imod)

 read(imod,"(a20)",err=1000,end=1000) string
 call charl2u(string)
 if(index(string,'[MOLDEN FORMAT]') > 0) then
   ifmt = 1
 else
   do while(.true.)
     read(imod,"(a61)",err=1000,end=1000) string
     if(index(string,'IOpCl') == 1) then
       ifmt = 2
       exit
     end if
   end do
 end if

 1000 continue

 return
end Subroutine dfiletest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! check 1c or 2c Hamiltonian in fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine reltest(imod,irel,string)
 implicit real(kind=8) (a-h,o-z)
 character*61 :: string

 rewind(imod)

 do while(.true.)
   read(imod,"(a61)",err=1000,end=1000) string
   if(index(string,'IOpCl') == 1) then
     read(string(45:61),*) i
     if(i == 6) irel = 2
     exit
   end if
 end do

 1000 continue

 return
end Subroutine reltest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! length of a string without the first and last spaces.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonspace(string)
 implicit real(kind=8) (a-h,o-z)
 character*(*)     :: string
 character*1       :: space

 space=' '
 length=LEN_TRIM(string)
 if(length == 0) then
  i=1
 else
  do i=1,length
    if(string(i:i) /= space) goto 20
  end do
 endif

 20  nonspace=i

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! b(mob+1:mob+ml,*) = a(moa+1:moa+ml,*)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine xcopy(ma,mb,n,moa,mob,ml,a,b)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: a(ma,n),b(mb,n)

 b(mob+1:mob+ml,1:n) = a(moa+1:moa+ml,1:n)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! b(*) = a(*)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine acopy(n,a,b)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: a(n),b(n)

 b = a

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! CHA --> cha
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charu2l(cha)
 implicit real(kind=8) (a-h,o-z)
 character*(*) :: cha
 character*1  :: U2L

 do i=1,len_trim(cha)
   cha(i:i)=U2L(cha(i:i))
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! cha --> CHA
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charl2u(cha)
 implicit real(kind=8) (a-h,o-z)
 character*(*) :: cha
 character*1  :: L2U

 do i=1,len_trim(cha)
   cha(i:i)=L2U(cha(i:i))
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! L --> l
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U2L(letter)
 implicit real(kind=8) (a-h,o-z)
 character*1 :: letter,U2L

 if((ichar(letter) >= 65).and.(ichar(letter) <= 90))then
   U2L=char(ichar(letter)+32)
 else
   U2L=letter
 endif

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! l --> L
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L2U(letter)
 implicit real(kind=8) (a-h,o-z)
 character*1 :: letter,L2U

 if( ichar(letter) >= 97 .and. ichar(letter) <= 122 )then
   L2U=char(ichar(letter)-32)
 else
   L2U=letter
 endif

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! remove duplicate primitive functions for a given L quantum number of an atom
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Rmdupfun(maxcn2, ngto, ncsh, idel, iptr, extmp, cotmp)
 use Constant, only : exptol,Zero
 implicit real(kind=8) (a-h,o-z)
 dimension         :: idel(ngto), iptr(2,ncsh), extmp(ngto), cotmp(maxcn2,ncsh)

 idel = 0

 do i=1,ngto-1
   if(idel(i) > 0) cycle
   do j=i+1,ngto
     if(idel(j) > 0) cycle
     if(abs(extmp(i)-extmp(j)) < exptol) idel(j) = i
   end do
 end do

 do i=1,ngto
   if(idel(i) /= 0) then
     k=idel(i)
     do j=1,ncsh
       if(abs(cotmp(i,j)) > Zero .and. abs(cotmp(k,j)) < exptol) cotmp(k,j) = cotmp(i,j)
     end do
   end if
 end do

 ngto1=0
 do i=1,ngto
   if(idel(i) == 0) then
     ngto1=ngto1+1
     extmp(ngto1) = extmp(i)
     cotmp(ngto1,:) = cotmp(i,:)
   end if
 end do

 ngto = ngto1
 iptr(1,:) = 1
 iptr(2,:) = ngto

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! check basis functions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine checkbas(iout,ia,lq,maxcn2,ngto,ncsh,extmp, cotmp)
 use Constant, only : exptol
 implicit real(kind=8) (a-h,o-z)
 dimension         :: extmp(ngto), cotmp(maxcn2,ncsh)

 ! check exponents
 do i=1,ngto
   if(extmp(i) < exptol) then
     write(iout,"(' ### Error: zero exponent ',d16.6,' of IA = ',i4,', LQ = ',i2)") extmp(i),ia,lq
     call estop(1)
   end if
 end do

 ! check contraction coefficients
 do i=1,ncsh
   x = abssum(ngto,cotmp(1,i))
   if(x < exptol) then
     write(iout,"(' ### Error: zero contraction coefficients in CGTO = ',i2,' of IA = ',i4,', LQ = ',i2)") i,ia,lq
     call estop(1)
   end if
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Do contraction (icase=0) or uncontraction (/=0) for scalar orbitals (irelc=1) or 2-component spinors (irelc=2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine drv_mocontra(icase,irelc,natom,maxlq,nc,np,nmo,mxla,nshlla,iptcon,gtocon,cmo,pmo,cmat)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: mxla(natom), nshlla(0:maxlq,natom,2), iptcon(0:maxlq,natom),gtocon(*),cmo(nc*irelc,nmo*irelc),  &
              pmo(np*irelc,nmo*irelc),cmat(np,nc,*)

 cmat(:,:,1) = Zero

 icoff = 0
 ipoff = 0
 do ia = 1,natom
   do iq = 0,mxla(ia)
     iml   = iq + iq + 1
     ishlp = nshlla(iq,ia,1)
     ishlc = nshlla(iq,ia,2)
     if(ishlp < 1 .or. ishlc < 1) cycle

     ic    = iptcon(iq,ia)
     call conmat(np,iq,ishlp,ishlc,iml,gtocon(ic),ipoff,cmat(1,icoff+1,1))

     icoff = icoff + ishlc*iml
     ipoff = ipoff + ishlp*iml
   end do
 end do

 if(irelc == 1) then

   if(icase == 0) then
     do imo = 1, nmo
       call MatMult(1,1,np,nc,One,Zero,pmo(1,imo),cmat,cmo(1,imo))
     end do
   else
     do imo = 1, nmo
       call MatMult(3,1,nc,np,One,Zero,cmo(1,imo),cmat,pmo(1,imo))
     end do
   end if

 else if(irelc == 2) then

   npd = np + np
   ncd = nc + nc
   
   call DblMat(np,nc,cmat(1,1,1),cmat(1,1,2))

   if(icase == 0) then
     do imo = 1, nmo+nmo
       call MatMult(1,1,npd,ncd,One,Zero,pmo(1,imo),cmat(1,1,2),cmo(1,imo))
     end do
   else
     do imo = 1, nmo+nmo
       call MatMult(3,1,ncd,npd,One,Zero,cmo(1,imo),cmat(1,1,2),pmo(1,imo))
     end do
   end if

 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Do contraction (icase=0) or uncontraction (/=0) for scalar density matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine drv_contra(icase,natom,maxlq,nc,np,mxla,nshlla,iptcon,gtocon,pc,pp,cmat,sc)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension :: mxla(natom), nshlla(0:maxlq,natom,2), iptcon(0:maxlq,natom),gtocon(*),pc(nc,nc),pp(np,np),cmat(np,nc),sc(*)

 cmat = Zero

 icoff = 0
 ipoff = 0
 do ia = 1,natom
   do iq = 0,mxla(ia)
     iml   = iq + iq + 1
     ishlp = nshlla(iq,ia,1)
     ishlc = nshlla(iq,ia,2)
     if(ishlp < 1 .or. ishlc < 1) cycle

     ic    = iptcon(iq,ia)
     call conmat(np,iq,ishlp,ishlc,iml,gtocon(ic),ipoff,cmat(1,icoff+1))

     icoff = icoff + ishlc*iml
     ipoff = ipoff + ishlp*iml
   end do
 end do

 if(icase == 0) then
   call docontr(np,nc,cmat,pc,pp,sc)
 else
   call douncon(np,nc,cmat,pc,pp,sc)
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Do contraction (icase=0) or uncontraction (/=0) for 2-component density matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine drv_contra2(icase,natom,maxlq,nc,np,mxla,nshlla,iptcon,gtocon,pc,pp,cmat,sc)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension :: mxla(natom), nshlla(0:maxlq,natom,2), iptcon(0:maxlq,natom),gtocon(*),pc(4*nc*nc,2),pp(4*np*np,2),  &
              cmat(4*np*nc),sc(np,*)

 sc(:,1:nc) = Zero

 icoff = 0
 ipoff = 0
 do ia = 1,natom
   do iq = 0,mxla(ia)
     iml   = iq + iq + 1
     ishlp = nshlla(iq,ia,1)
     ishlc = nshlla(iq,ia,2)
     if(ishlp < 1 .or. ishlc < 1) cycle

     ic    = iptcon(iq,ia)
     call conmat(np,iq,ishlp,ishlc,iml,gtocon(ic),ipoff,sc(1,icoff+1))

     icoff = icoff + ishlc*iml
     ipoff = ipoff + ishlp*iml
   end do
 end do

 call DblMat(np,nc,sc,cmat)

 if(icase == 0) then
   call docontr(np*2,nc*2,cmat,pc(1,1),pp(1,1),sc)
   call docontr(np*2,nc*2,cmat,pc(1,2),pp(1,2),sc)
 else
   call douncon(np*2,nc*2,cmat,pc(1,1),pp(1,1),sc)
   call douncon(np*2,nc*2,cmat,pc(1,2),pp(1,2),sc)
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Do contraction for a density matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DoContr(NP,NC,Con,Pc,Pp,Sc)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: Con(NP,NC),Pc(*),Pp(*),Sc(*)

 call dgemm("t","n",NC,NP,NP,One,Con,NP,Pp, NP,Zero,Sc, NC)
 call dgemm("n","n",NC,NC,NP,One,Sc, NC,Con,NP,Zero,Pc, NC)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Do uncontraction for a density matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DoUnCon(NP,NC,Con,Pc,Pp,Sc)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: Con(NP,NC),Pc(*),Pp(*),Sc(*)

 call dgemm("n","n",NP,NC,NC,One,Con,NP,Pc, NC,Zero,Sc, NP)
 call dgemm("n","t",NP,NP,NC,One,Sc, NP,Con,NP,Zero,Pp, NP)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! make contraction coefficient matrix for a given L-type
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine conmat(np,lq,nshlp,nshlc,nml,coeff,ioff,con)
 implicit real(kind=8) (a-h,o-z)
 dimension :: coeff(nshlp,nshlc),con(np,nml,nshlc)

 do i = 1, nshlc
   m = 0
   do j = 1, nshlp
     n = 0
     do k = -lq, lq
       m = m + 1
       n = n + 1
       con(ioff+m,n,i) = coeff(j,i)
     end do
   end do
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! fill in the C2S (icase=0, cart2pure) or S2C (icase/=0, pure2cart) coefficient array
! mpatt : pattern mode of basis functions. 0: Gaussian, -1: Molden
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine FilC2S(icase,mpatt,MaxL,C2SCoef)
 implicit real(kind=8) (a-h,o-z)
 dimension :: C2SCoef(*)

 ioff=0
 do L=0,MaxL
   NC = (L+1)*(L+2)/2
   if(icase == 0) then
     if(mpatt == 0) then
       call CoefC2S_GAUSS(L, NC, C2SCoef(ioff+1))
     else if(mpatt == -1) then
       call CoefC2S_MOLDEN(L, NC, C2SCoef(ioff+1))
     else
       write(*,"(/,' ### Error: unknown MPATT!')")
       call estop(1)
     end if
   else
     if(mpatt == 0) then
       call CoefS2C_GAUSS(L, NC, C2SCoef(ioff+1))
     else if(mpatt == -1) then
       call CoefS2C_MOLDEN(L, NC, C2SCoef(ioff+1))
     else
       write(*,"(/,' ### Error: unknown MPATT!')")
       call estop(1)
     end if
   end if
   ioff =(NC*(L+3)*(L*3+2))/6
 enddo

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! transformation coefficients from normalized Cartesian functions to spherical functions.
!
! Input
!   Lq      angular momentum
!   NC      = (Lq+1)*(Lq+2)/2
!
! Output
!   C2S     transformation coefficients
!
! NOTE
!   MOLDEN ordering in Cartesian & spherical functions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CoefC2S_MOLDEN(Lq, NC, C2S)
  implicit real(kind=8) (a-h,o-z)
  dimension :: C2S(NC,1+Lq+Lq)

  C2S = 0.0d0

  select case (Lq)
    case (0)
      C2S( 1, 1) = 1.0D0
    case (1)
      C2S( 1, 1) = 1.0D0
      C2S( 2, 2) = 1.0D0
      C2S( 3, 3) = 1.0D0
    case (2)
      C2S( 1, 1) =-0.5D+00
      C2S( 2, 1) =-0.5D+00
      C2S( 3, 1) = 1.0D0
      C2S( 5, 2) = 1.0D0
      C2S( 6, 3) = 1.0D0
      C2S( 1, 4) = 0.8660254037844387D+00
      C2S( 2, 4) =-0.8660254037844387D+00
      C2S( 4, 5) = 1.0D0
    case (3)
      C2S( 3, 1) = 1.0D0
      C2S( 6, 1) =-0.6708203932499369D+00
      C2S( 9, 1) =-0.6708203932499369D+00
      C2S( 1, 2) =-0.6123724356957946D+00
      C2S( 4, 2) =-0.2738612787525831D+00
      C2S( 7, 2) = 0.1095445115010332D+01
      C2S( 2, 3) =-0.6123724356957946D+00
      C2S( 5, 3) =-0.2738612787525831D+00
      C2S( 8, 3) = 0.1095445115010332D+01
      C2S( 6, 4) = 0.8660254037844387D+00
      C2S( 9, 4) =-0.8660254037844387D+00
      C2S(10, 5) = 1.0D0
      C2S( 1, 6) = 0.7905694150420949D+00
      C2S( 4, 6) =-0.1060660171779821D+01
      C2S( 2, 7) =-0.7905694150420949D+00
      C2S( 5, 7) = 0.1060660171779821D+01
    case (4)
      C2S( 1, 1) = 0.375D+00
      C2S( 2, 1) = 0.375D+00
      C2S( 3, 1) = 1.0D0
      C2S(10, 1) = 0.2195775164134200D+00
      C2S(11, 1) =-0.8783100656536799D+00
      C2S(12, 1) =-0.8783100656536799D+00
      C2S( 5, 2) =-0.8964214570007952D+00
      C2S( 8, 2) = 0.1195228609334394D+01
      C2S(14, 2) =-0.4008918628686366D+00
      C2S( 7, 3) =-0.8964214570007952D+00
      C2S( 9, 3) = 0.1195228609334394D+01
      C2S(13, 3) =-0.4008918628686366D+00
      C2S( 1, 4) =-0.5590169943749475D+00
      C2S( 2, 4) = 0.5590169943749475D+00
      C2S(11, 4) = 0.9819805060619657D+00
      C2S(12, 4) =-0.9819805060619657D+00
      C2S( 4, 5) =-0.4225771273642583D+00
      C2S( 6, 5) =-0.4225771273642583D+00
      C2S(15, 5) = 0.1133893419027682D+01
      C2S( 5, 6) = 0.7905694150420949D+00
      C2S(14, 6) =-0.1060660171779821D+01
      C2S( 7, 7) =-0.7905694150420949D+00
      C2S(13, 7) = 0.1060660171779821D+01
      C2S( 1, 8) = 0.7395099728874520D+00
      C2S( 2, 8) = 0.7395099728874520D+00
      C2S(10, 8) =-0.1299038105676658D+01
      C2S( 4, 9) = 0.1118033988749895D+01
      C2S( 6, 9) =-0.1118033988749895D+01
    case (5)
      C2S( 1, 1) = 1.0D0
      C2S( 3, 1) =-0.1091089451179962D+01
      C2S( 5, 1) = 0.625D+00
      C2S(12, 1) =-0.1091089451179962D+01
      C2S(14, 1) = 0.3659625273557000D+00
      C2S(19, 1) = 0.625D+00
      C2S( 7, 2) = 0.1290994448735806D+01
      C2S( 9, 2) =-0.5669467095138407D+00
      C2S(11, 2) = 0.1613743060919757D+00
      C2S(16, 2) =-0.1267731382092775D+01
      C2S(18, 2) = 0.2112885636821291D+00
      C2S(21, 2) = 0.4841229182759271D+00
      C2S( 2, 3) = 0.1290994448735806D+01
      C2S( 4, 3) =-0.1267731382092775D+01
      C2S( 6, 3) = 0.4841229182759271D+00
      C2S(13, 3) =-0.5669467095138407D+00
      C2S(15, 3) = 0.2112885636821291D+00
      C2S(20, 3) = 0.1613743060919757D+00
      C2S( 3, 4) =-0.1118033988749895D+01
      C2S( 5, 4) = 0.8539125638299665D+00
      C2S(12, 4) = 0.1118033988749895D+01
      C2S(19, 4) =-0.8539125638299665D+00
      C2S( 8, 5) = 0.1290994448735806D+01
      C2S(10, 5) =-0.6454972243679028D+00
      C2S(17, 5) =-0.6454972243679028D+00
      C2S( 9, 6) =-0.1224744871391589D+01
      C2S(11, 6) = 0.5229125165837972D+00
      C2S(16, 6) = 0.9128709291752769D+00
      C2S(18, 6) = 0.2282177322938192D+00
      C2S(21, 6) =-0.5229125165837972D+00
      C2S( 4, 7) =-0.9128709291752769D+00
      C2S( 6, 7) = 0.5229125165837972D+00
      C2S(13, 7) = 0.1224744871391589D+01
      C2S(15, 7) =-0.2282177322938192D+00
      C2S(20, 7) =-0.5229125165837972D+00
      C2S( 5, 8) = 0.7395099728874520D+00
      C2S(14, 8) =-0.1299038105676658D+01
      C2S(19, 8) = 0.7395099728874520D+00
      C2S(10, 9) =-0.1118033988749895D+01
      C2S(17, 9) = 0.1118033988749895D+01
      C2S(11,10) = 0.1169267933366857D+01
      C2S(18,10) =-0.1530931089239487D+01
      C2S(21,10) = 0.7015607600201140D+00
      C2S( 6,11) = 0.7015607600201140D+00
      C2S(15,11) =-0.1530931089239487D+01
      C2S(20,11) = 0.1169267933366857D+01
    case default
      write(*,"(/,' CoefC2S_MOLDEN is available only for Lq = 0 ~ 5!')")
      call estop(1)
  end select

!  do i=1,nc
!    write(77,"(100f10.5)")C2S(i,:)
!  end do

  return
end subroutine CoefC2S_MOLDEN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! transformation coefficients from normalized spherical functions to Cartesian functions.
!
! Input
!   Lq      angular momentum
!   NC      = (Lq+1)*(Lq+2)/2
!
! Output
!   S2C     transformation coefficients
!
! NOTE
!   MOLDEN ordering in Cartesian & spherical functions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CoefS2C_MOLDEN(Lq, NC, S2C)
  implicit real(kind=8) (a-h,o-z)
  dimension :: S2C(NC,1+Lq+Lq)

  S2C = 0.0d0

  select case (Lq)
    case (0)
      S2C( 1, 1) = 1.0D0
    case (1)
      S2C( 1, 1) = 1.0D0
      S2C( 2, 2) = 1.0D0
      S2C( 3, 3) = 1.0D0
    case (2)
      S2C( 1, 1) =-0.3333333333333333D+00
      S2C( 2, 1) =-0.3333333333333333D+00
      S2C( 3, 1) = 0.6666666666666666D+00
      S2C( 5, 2) = 1.0D0
      S2C( 6, 3) = 1.0D0
      S2C( 1, 4) = 0.5773502691896257D+00
      S2C( 2, 4) =-0.5773502691896257D+00
      S2C( 4, 5) = 1.0D0
    case (3)
      S2C( 3, 1) = 0.4D0
      S2C( 6, 1) =-0.4472135954999578D+00
      S2C( 9, 1) =-0.4472135954999578D+00
      S2C( 1, 2) =-0.2449489742783178D+00
      S2C( 4, 2) =-0.1825741858350554D+00
      S2C( 7, 2) = 0.7302967433402215D+00
      S2C( 2, 3) =-0.2449489742783178D+00
      S2C( 5, 3) =-0.1825741858350554D+00
      S2C( 8, 3) = 0.7302967433402214D+00
      S2C( 6, 4) = 0.5773502691896257D+00
      S2C( 9, 4) =-0.5773502691896257D+00
      S2C(10, 5) = 1.0D0
      S2C( 1, 6) = 0.3162277660168379D+00
      S2C( 4, 6) =-0.7071067811865475D+00
      S2C( 2, 7) =-0.3162277660168379D+00
      S2C( 5, 7) = 0.7071067811865475D+00
    case (4)
      S2C( 1, 1) = 0.8571428571428570D-01
      S2C( 2, 1) = 0.8571428571428572D-01
      S2C( 3, 1) = 0.2285714285714286D+00
      S2C(10, 1) = 0.9759000729485329D-01
      S2C(11, 1) =-0.3903600291794133D+00
      S2C(12, 1) =-0.3903600291794133D+00
      S2C( 5, 2) =-0.3585685828003183D+00
      S2C( 8, 2) = 0.4780914437337576D+00
      S2C(14, 2) =-0.2672612419124243D+00
      S2C( 7, 3) =-0.3585685828003183D+00
      S2C( 9, 3) = 0.4780914437337576D+00
      S2C(13, 3) =-0.2672612419124243D+00
      S2C( 1, 4) =-0.1277753129999880D+00
      S2C( 2, 4) = 0.1277753129999880D+00
      S2C(11, 4) = 0.4364357804719848D+00
      S2C(12, 4) =-0.4364357804719847D+00
      S2C( 4, 5) =-0.1690308509457035D+00
      S2C( 6, 5) =-0.1690308509457035D+00
      S2C(15, 5) = 0.7559289460184542D+00
      S2C( 5, 6) = 0.3162277660168382D+00
      S2C(14, 6) =-0.7071067811865474D+00
      S2C( 7, 7) =-0.3162277660168382D+00
      S2C(13, 7) = 0.7071067811865475D+00
      S2C( 1, 8) = 0.1690308509457033D+00
      S2C( 2, 8) = 0.1690308509457032D+00
      S2C(10, 8) =-0.5773502691896257D+00
      S2C( 4, 9) = 0.4472135954999581D+00
      S2C( 6, 9) =-0.4472135954999581D+00
    case (5)
      S2C( 1, 1) = 0.1269841269841268D+00
      S2C( 3, 1) =-0.2909571869813233D+00
      S2C( 5, 1) = 0.1428571428571428D+00
      S2C(12, 1) =-0.2909571869813233D+00
      S2C(14, 1) = 0.1626500121580887D+00
      S2C(19, 1) = 0.1428571428571429D+00
      S2C( 7, 2) = 0.2950844454253272D+00
      S2C( 9, 2) =-0.2519763153394846D+00
      S2C(11, 2) = 0.3688555567816591D-01
      S2C(16, 2) =-0.3380617018914065D+00
      S2C(18, 2) = 0.5634361698190118D-01
      S2C(21, 2) = 0.6147592613027647D-01
      S2C( 2, 3) = 0.2950844454253272D+00
      S2C( 4, 3) =-0.3380617018914065D+00
      S2C( 6, 3) = 0.6147592613027650D-01
      S2C(13, 3) =-0.2519763153394846D+00
      S2C(15, 3) = 0.5634361698190118D-01
      S2C(20, 3) = 0.3688555567816593D-01
      S2C( 3, 4) =-0.2981423969999721D+00
      S2C( 5, 4) = 0.1951800145897066D+00
      S2C(12, 4) = 0.2981423969999720D+00
      S2C(19, 4) =-0.1951800145897067D+00
      S2C( 8, 5) = 0.5163977794943222D+00
      S2C(10, 5) =-0.2581988897471612D+00
      S2C(17, 5) =-0.2581988897471612D+00
      S2C( 9, 6) =-0.5443310539518175D+00
      S2C(11, 6) = 0.1195228609334393D+00
      S2C(16, 6) = 0.2434322477800738D+00
      S2C(18, 6) = 0.6085806194501837D-01
      S2C(21, 6) =-0.6640158940746627D-01
      S2C( 4, 7) =-0.2434322477800739D+00
      S2C( 6, 7) = 0.6640158940746628D-01
      S2C(13, 7) = 0.5443310539518175D+00
      S2C(15, 7) =-0.6085806194501844D-01
      S2C(20, 7) =-0.1195228609334393D+00
      S2C( 5, 8) = 0.1690308509457033D+00
      S2C(14, 8) =-0.5773502691896257D+00
      S2C(19, 8) = 0.1690308509457032D+00
      S2C(10, 9) =-0.4472135954999579D+00
      S2C(17, 9) = 0.4472135954999580D+00
      S2C(11,10) = 0.2672612419124244D+00
      S2C(18,10) =-0.4082482904638632D+00
      S2C(21,10) = 0.8908708063747461D-01
      S2C( 6,11) = 0.8908708063747466D-01
      S2C(15,11) =-0.4082482904638631D+00
      S2C(20,11) = 0.2672612419124244D+00
    case default
      write(*,"(/,' CoefS2C_MOLDEN is available only for Lq = 0 ~ 5!')")
      call estop(1)
  end select

!  do i=1,nc
!    write(77,"(100f10.5)")S2C(i,:)
!  end do

  return
end subroutine CoefS2C_MOLDEN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! transformation coefficients from normalized Cartesian functions to spherical functions.
!
! Input
!   Lq      angular momentum
!   NC      = (Lq+1)*(Lq+2)/2
!
! Output
!   C2S     transformation coefficients
!
! NOTE
!   Gaussian ordering in Cartesian & spherical functions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CoefC2S_GAUSS(Lq, NC, C2S)
  implicit real(kind=8) (a-h,o-z)
  dimension :: C2S(NC,1+Lq+Lq)

  C2S = 0.0d0

  select case (Lq)
    case (0)
      C2S( 1, 1) = 1.0D0
    case (1)
      C2S( 1, 1) = 1.0D0
      C2S( 2, 2) = 1.0D0
      C2S( 3, 3) = 1.0D0
    case (2)
      C2S( 1, 1) =-0.5D+00
      C2S( 2, 1) =-0.5D+00
      C2S( 3, 1) = 1.0D0
      C2S( 5, 2) = 1.0D0
      C2S( 6, 3) = 1.0D0
      C2S( 1, 4) = 0.8660254037844387D+00
      C2S( 2, 4) =-0.8660254037844387D+00
      C2S( 4, 5) = 1.0D0
    case (3)
      C2S( 3, 1) = 1.0D0
      C2S( 6, 1) =-0.6708203932499369D+00
      C2S( 9, 1) =-0.6708203932499369D+00
      C2S( 1, 2) =-0.6123724356957946D+00
      C2S( 4, 2) =-0.2738612787525831D+00
      C2S( 7, 2) = 0.1095445115010332D+01
      C2S( 2, 3) =-0.6123724356957946D+00
      C2S( 5, 3) =-0.2738612787525831D+00
      C2S( 8, 3) = 0.1095445115010332D+01
      C2S( 6, 4) = 0.8660254037844387D+00
      C2S( 9, 4) =-0.8660254037844387D+00
      C2S(10, 5) = 1.0D0
      C2S( 1, 6) = 0.7905694150420949D+00
      C2S( 4, 6) =-0.1060660171779821D+01
      C2S( 2, 7) =-0.7905694150420949D+00
      C2S( 5, 7) = 0.1060660171779821D+01
    case (4)
      C2S( 1, 1) = 1.0D0
      C2S( 3, 1) =-0.8783100656536799D+00
      C2S( 5, 1) = 0.375D+00
      C2S(10, 1) =-0.8783100656536799D+00
      C2S(12, 1) = 0.2195775164134200D+00
      C2S(15, 1) = 0.375D+00
      C2S( 6, 2) = 0.1195228609334394D+01
      C2S( 8, 2) =-0.4008918628686366D+00
      C2S(13, 2) =-0.8964214570007952D+00
      C2S( 2, 3) = 0.1195228609334394D+01
      C2S( 4, 3) =-0.8964214570007952D+00
      C2S(11, 3) =-0.4008918628686366D+00
      C2S( 3, 4) =-0.9819805060619657D+00
      C2S( 5, 4) = 0.5590169943749475D+00
      C2S(10, 4) = 0.9819805060619657D+00
      C2S(15, 4) =-0.5590169943749475D+00
      C2S( 7, 5) = 0.1133893419027682D+01
      C2S( 9, 5) =-0.4225771273642583D+00
      C2S(14, 5) =-0.4225771273642583D+00
      C2S( 8, 6) =-0.1060660171779821D+01
      C2S(13, 6) = 0.7905694150420949D+00
      C2S( 4, 7) =-0.7905694150420949D+00
      C2S(11, 7) = 0.1060660171779821D+01
      C2S( 5, 8) = 0.7395099728874520D+00
      C2S(12, 8) =-0.1299038105676658D+01
      C2S(15, 8) = 0.7395099728874520D+00
      C2S( 9, 9) =-0.1118033988749895D+01
      C2S(14, 9) = 0.1118033988749895D+01
    case (5)
      C2S( 1, 1) = 1.0D0
      C2S( 3, 1) =-0.1091089451179962D+01
      C2S( 5, 1) = 0.625D+00
      C2S(12, 1) =-0.1091089451179962D+01
      C2S(14, 1) = 0.3659625273557000D+00
      C2S(19, 1) = 0.625D+00
      C2S( 7, 2) = 0.1290994448735806D+01
      C2S( 9, 2) =-0.5669467095138409D+00
      C2S(11, 2) = 0.1613743060919757D+00
      C2S(16, 2) =-0.1267731382092775D+01
      C2S(18, 2) = 0.2112885636821291D+00
      C2S(21, 2) = 0.4841229182759271D+00
      C2S( 2, 3) = 0.1290994448735806D+01
      C2S( 4, 3) =-0.1267731382092775D+01
      C2S( 6, 3) = 0.4841229182759271D+00
      C2S(13, 3) =-0.5669467095138409D+00
      C2S(15, 3) = 0.2112885636821291D+00
      C2S(20, 3) = 0.1613743060919757D+00
      C2S( 3, 4) =-0.1118033988749895D+01
      C2S( 5, 4) = 0.8539125638299665D+00
      C2S(12, 4) = 0.1118033988749895D+01
      C2S(19, 4) =-0.8539125638299665D+00
      C2S( 8, 5) = 0.1290994448735806D+01
      C2S(10, 5) =-0.6454972243679028D+00
      C2S(17, 5) =-0.6454972243679028D+00
      C2S( 9, 6) =-0.1224744871391589D+01
      C2S(11, 6) = 0.5229125165837972D+00
      C2S(16, 6) = 0.9128709291752769D+00
      C2S(18, 6) = 0.2282177322938192D+00
      C2S(21, 6) =-0.5229125165837972D+00
      C2S( 4, 7) =-0.9128709291752769D+00
      C2S( 6, 7) = 0.5229125165837972D+00
      C2S(13, 7) = 0.1224744871391589D+01
      C2S(15, 7) =-0.2282177322938192D+00
      C2S(20, 7) =-0.5229125165837972D+00
      C2S( 5, 8) = 0.7395099728874520D+00
      C2S(14, 8) =-0.1299038105676658D+01
      C2S(19, 8) = 0.7395099728874520D+00
      C2S(10, 9) =-0.1118033988749895D+01
      C2S(17, 9) = 0.1118033988749895D+01
      C2S(11,10) = 0.1169267933366857D+01
      C2S(18,10) =-0.1530931089239486D+01
      C2S(21,10) = 0.7015607600201140D+00
      C2S( 6,11) = 0.7015607600201140D+00
      C2S(15,11) =-0.1530931089239486D+01
      C2S(20,11) = 0.1169267933366857D+01
    case (6)
      C2S( 1, 1) = 1.0D0
      C2S( 3, 1) =-0.1305582419667734D+01
      C2S( 5, 1) = 0.9791868147508004D+00
      C2S( 7, 1) =-0.3125D+00
      C2S(14, 1) =-0.1305582419667734D+01
      C2S(16, 1) = 0.5733530903673287D+00
      C2S(18, 1) =-0.1631978024584667D+00
      C2S(23, 1) = 0.9791868147508004D+00
      C2S(25, 1) =-0.1631978024584667D+00
      C2S(28, 1) =-0.3125D+00
      C2S( 8, 2) = 0.1381698559415514D+01
      C2S(10, 2) =-0.7537783614444091D+00
      C2S(12, 2) = 0.2878538665448990D+00
      C2S(19, 2) =-0.1685499656158105D+01
      C2S(21, 2) = 0.3768891807222046D+00
      C2S(26, 2) = 0.8635615996346969D+00
      C2S( 2, 3) = 0.1381698559415514D+01
      C2S( 4, 3) =-0.1685499656158105D+01
      C2S( 6, 3) = 0.8635615996346969D+00
      C2S(15, 3) =-0.7537783614444091D+00
      C2S(17, 3) = 0.3768891807222046D+00
      C2S(24, 3) = 0.2878538665448990D+00
      C2S( 3, 4) =-0.1261312447773783D+01
      C2S( 5, 4) = 0.1261312447773783D+01
      C2S( 7, 4) =-0.4528555233184200D+00
      C2S(14, 4) = 0.1261312447773783D+01
      C2S(18, 4) =-0.7883202798586142D-01
      C2S(23, 4) =-0.1261312447773783D+01
      C2S(25, 4) = 0.7883202798586142D-01
      C2S(28, 4) = 0.4528555233184200D+00
      C2S( 9, 5) = 0.1456438162508838D+01
      C2S(11, 5) =-0.9534625892455924D+00
      C2S(13, 5) = 0.2730821554704072D+00
      C2S(20, 5) =-0.9534625892455924D+00
      C2S(22, 5) = 0.2665008954445131D+00
      C2S(27, 5) = 0.2730821554704072D+00
      C2S(10, 6) =-0.1430193883868389D+01
      C2S(12, 6) = 0.8192464664112216D+00
      C2S(19, 6) = 0.1066003581778052D+01
      C2S(21, 6) = 0.3575484709670972D+00
      C2S(26, 6) =-0.8192464664112216D+00
      C2S( 4, 7) =-0.1066003581778052D+01
      C2S( 6, 7) = 0.8192464664112216D+00
      C2S(15, 7) = 0.1430193883868389D+01
      C2S(17, 7) =-0.3575484709670971D+00
      C2S(24, 7) =-0.8192464664112216D+00
      C2S( 5, 8) = 0.8635615996346969D+00
      C2S( 7, 8) =-0.4960783708246108D+00
      C2S(16, 8) =-0.1516949690542295D+01
      C2S(18, 8) = 0.4317807998173485D+00
      C2S(23, 8) = 0.8635615996346969D+00
      C2S(25, 8) = 0.4317807998173485D+00
      C2S(28, 8) =-0.4960783708246108D+00
      C2S(11, 9) =-0.1305582419667734D+01
      C2S(13, 9) = 0.5982930264130992D+00
      C2S(20, 9) = 0.1305582419667734D+01
      C2S(27, 9) =-0.5982930264130992D+00
      C2S(12,10) = 0.1169267933366857D+01
      C2S(21,10) =-0.1530931089239487D+01
      C2S(26,10) = 0.7015607600201140D+00
      C2S( 6,11) = 0.7015607600201140D+00
      C2S(17,11) =-0.1530931089239487D+01
      C2S(24,11) = 0.1169267933366857D+01
      C2S( 7,12) =-0.6716932893813962D+00
      C2S(18,12) = 0.1753901900050285D+01
      C2S(25,12) =-0.1753901900050285D+01
      C2S(28,12) = 0.6716932893813962D+00
      C2S(13,13) = 0.1215138880951474D+01
      C2S(22,13) =-0.1976423537605237D+01
      C2S(27,13) = 0.1215138880951474D+01
    case default
      write(*,"(/,' CoefC2S_GAUSS has been programmed only for Lq = 0 ~ 6!')")
      call estop(1)
  end select

!  do i=1,nc
!    write(77,"(100f10.5)")C2S(i,:)
!  end do

  return
end subroutine CoefC2S_GAUSS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! transformation coefficients from normalized spherical functions to Cartesian functions.
!
! Input
!   Lq      angular momentum
!   NC      = (Lq+1)*(Lq+2)/2
!
! Output
!   S2C     transformation coefficients
!
! NOTE
!   Gaussian ordering in Cartesian & spherical functions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CoefS2C_GAUSS(Lq, NC, S2C)
  implicit real(kind=8) (a-h,o-z)
  dimension :: S2C(NC,1+Lq+Lq)

  S2C = 0.0d0

  select case (Lq)
    case (0)
      S2C( 1, 1) = 1.0D0
    case (1)
      S2C( 1, 1) = 1.0D0
      S2C( 2, 2) = 1.0D0
      S2C( 3, 3) = 1.0D0
    case (2)
      S2C( 1, 1) =-0.3333333333333333D+00
      S2C( 2, 1) =-0.3333333333333333D+00
      S2C( 3, 1) = 0.6666666666666666D+00
      S2C( 5, 2) = 1.0D0
      S2C( 6, 3) = 1.0D0
      S2C( 1, 4) = 0.5773502691896257D+00
      S2C( 2, 4) =-0.5773502691896257D+00
      S2C( 4, 5) = 1.0D0
    case (3)
      S2C( 3, 1) = 0.4D0
      S2C( 6, 1) =-0.4472135954999578D+00
      S2C( 9, 1) =-0.4472135954999578D+00
      S2C( 1, 2) =-0.2449489742783178D+00
      S2C( 4, 2) =-0.1825741858350554D+00
      S2C( 7, 2) = 0.7302967433402215D+00
      S2C( 2, 3) =-0.2449489742783178D+00
      S2C( 5, 3) =-0.1825741858350554D+00
      S2C( 8, 3) = 0.7302967433402214D+00
      S2C( 6, 4) = 0.5773502691896257D+00
      S2C( 9, 4) =-0.5773502691896257D+00
      S2C(10, 5) = 1.0D0
      S2C( 1, 6) = 0.3162277660168379D+00
      S2C( 4, 6) =-0.7071067811865475D+00
      S2C( 2, 7) =-0.3162277660168379D+00
      S2C( 5, 7) = 0.7071067811865475D+00
    case (4)
      S2C( 1, 1) = 0.2285714285714285D+00
      S2C( 3, 1) =-0.3903600291794133D+00
      S2C( 5, 1) = 0.8571428571428572D-01
      S2C(10, 1) =-0.3903600291794133D+00
      S2C(12, 1) = 0.9759000729485334D-01
      S2C(15, 1) = 0.8571428571428574D-01
      S2C( 6, 2) = 0.4780914437337576D+00
      S2C( 8, 2) =-0.2672612419124244D+00
      S2C(13, 2) =-0.3585685828003180D+00
      S2C( 2, 3) = 0.4780914437337576D+00
      S2C( 4, 3) =-0.3585685828003181D+00
      S2C(11, 3) =-0.2672612419124244D+00
      S2C( 3, 4) =-0.4364357804719847D+00
      S2C( 5, 4) = 0.1277753129999880D+00
      S2C(10, 4) = 0.4364357804719847D+00
      S2C(15, 4) =-0.1277753129999880D+00
      S2C( 7, 5) = 0.7559289460184544D+00
      S2C( 9, 5) =-0.1690308509457033D+00
      S2C(14, 5) =-0.1690308509457033D+00
      S2C( 8, 6) =-0.7071067811865477D+00
      S2C(13, 6) = 0.3162277660168379D+00
      S2C( 4, 7) =-0.3162277660168379D+00
      S2C(11, 7) = 0.7071067811865477D+00
      S2C( 5, 8) = 0.1690308509457033D+00
      S2C(12, 8) =-0.5773502691896257D+00
      S2C(15, 8) = 0.1690308509457034D+00
      S2C( 9, 9) =-0.4472135954999580D+00
      S2C(14, 9) = 0.4472135954999580D+00
    case (5)
      S2C( 1, 1) = 0.1269841269841270D+00
      S2C( 3, 1) =-0.2909571869813230D+00
      S2C( 5, 1) = 0.1428571428571429D+00
      S2C(12, 1) =-0.2909571869813229D+00
      S2C(14, 1) = 0.1626500121580890D+00
      S2C(19, 1) = 0.1428571428571429D+00
      S2C( 7, 2) = 0.2950844454253270D+00
      S2C( 9, 2) =-0.2519763153394848D+00
      S2C(11, 2) = 0.3688555567816585D-01
      S2C(16, 2) =-0.3380617018914067D+00
      S2C(18, 2) = 0.5634361698190108D-01
      S2C(21, 2) = 0.6147592613027641D-01
      S2C( 2, 3) = 0.2950844454253269D+00
      S2C( 4, 3) =-0.3380617018914067D+00
      S2C( 6, 3) = 0.6147592613027646D-01
      S2C(13, 3) =-0.2519763153394848D+00
      S2C(15, 3) = 0.5634361698190109D-01
      S2C(20, 3) = 0.3688555567816588D-01
      S2C( 3, 4) =-0.2981423969999720D+00
      S2C( 5, 4) = 0.1951800145897066D+00
      S2C(12, 4) = 0.2981423969999720D+00
      S2C(19, 4) =-0.1951800145897066D+00
      S2C( 8, 5) = 0.5163977794943222D+00
      S2C(10, 5) =-0.2581988897471611D+00
      S2C(17, 5) =-0.2581988897471611D+00
      S2C( 9, 6) =-0.5443310539518175D+00
      S2C(11, 6) = 0.1195228609334393D+00
      S2C(16, 6) = 0.2434322477800739D+00
      S2C(18, 6) = 0.6085806194501838D-01
      S2C(21, 6) =-0.6640158940746632D-01
      S2C( 4, 7) =-0.2434322477800739D+00
      S2C( 6, 7) = 0.6640158940746634D-01
      S2C(13, 7) = 0.5443310539518175D+00
      S2C(15, 7) =-0.6085806194501847D-01
      S2C(20, 7) =-0.1195228609334393D+00
      S2C( 5, 8) = 0.1690308509457033D+00
      S2C(14, 8) =-0.5773502691896257D+00
      S2C(19, 8) = 0.1690308509457034D+00
      S2C(10, 9) =-0.4472135954999580D+00
      S2C(17, 9) = 0.4472135954999580D+00
      S2C(11,10) = 0.2672612419124245D+00
      S2C(18,10) =-0.4082482904638630D+00
      S2C(21,10) = 0.8908708063747484D-01
      S2C( 6,11) = 0.8908708063747484D-01
      S2C(15,11) =-0.4082482904638629D+00
      S2C(20,11) = 0.2672612419124245D+00
    case (6)
      S2C( 1, 1) = 0.6926406926406940D-01
      S2C( 3, 1) =-0.1989458925207974D+00
      S2C( 5, 1) = 0.1492094193905983D+00
      S2C( 7, 1) =-0.2164502164502157D-01
      S2C(14, 1) =-0.1989458925207974D+00
      S2C(16, 1) = 0.1698823971458754D+00
      S2C(18, 1) =-0.2486823656509965D-01
      S2C(23, 1) = 0.1492094193905982D+00
      S2C(25, 1) =-0.2486823656509963D-01
      S2C(28, 1) =-0.2164502164502158D-01
      S2C( 8, 2) = 0.1754537853226046D+00
      S2C(10, 2) =-0.2010075630518427D+00
      S2C(12, 2) = 0.6579516949597686D-01
      S2C(19, 2) =-0.2696799449852970D+00
      S2C(21, 2) = 0.1005037815259212D+00
      S2C(26, 2) = 0.1096586158266282D+00
      S2C( 2, 3) = 0.1754537853226046D+00
      S2C( 4, 3) =-0.2696799449852971D+00
      S2C( 6, 3) = 0.1096586158266281D+00
      S2C(15, 3) =-0.2010075630518426D+00
      S2C(17, 3) = 0.1005037815259212D+00
      S2C(24, 3) = 0.6579516949597686D-01
      S2C( 3, 4) =-0.1921999920417193D+00
      S2C( 5, 4) = 0.1921999920417193D+00
      S2C( 7, 4) =-0.3136661633374332D-01
      S2C(14, 4) = 0.1921999920417192D+00
      S2C(18, 4) =-0.1201249950260741D-01
      S2C(23, 4) =-0.1921999920417192D+00
      S2C(25, 4) = 0.1201249950260746D-01
      S2C(28, 4) = 0.3136661633374332D-01
      S2C( 9, 5) = 0.3329001514305916D+00
      S2C(11, 5) =-0.2542566904654913D+00
      S2C(13, 5) = 0.3467709910735329D-01
      S2C(20, 5) =-0.2542566904654912D+00
      S2C(22, 5) = 0.4264014327112210D-01
      S2C(27, 5) = 0.3467709910735328D-01
      S2C(10, 6) =-0.3813850356982369D+00
      S2C(12, 6) = 0.1872563351797078D+00
      S2C(19, 6) = 0.1705605730844883D+00
      S2C(21, 6) = 0.9534625892455936D-01
      S2C(26, 6) =-0.1040312973220598D+00
      S2C( 4, 7) =-0.1705605730844883D+00
      S2C( 6, 7) = 0.1040312973220599D+00
      S2C(15, 7) = 0.3813850356982370D+00
      S2C(17, 7) =-0.9534625892455920D-01
      S2C(24, 7) =-0.1872563351797079D+00
      S2C( 5, 8) = 0.1315903389919539D+00
      S2C( 7, 8) =-0.3436040663720241D-01
      S2C(16, 8) =-0.4494665749754946D+00
      S2C(18, 8) = 0.6579516949597697D-01
      S2C(23, 8) = 0.1315903389919539D+00
      S2C(25, 8) = 0.6579516949597695D-01
      S2C(28, 8) =-0.3436040663720241D-01
      S2C(11, 9) =-0.3481553119113957D+00
      S2C(13, 9) = 0.7597371763975860D-01
      S2C(20, 9) = 0.3481553119113958D+00
      S2C(27, 9) =-0.7597371763975858D-01
      S2C(12,10) = 0.2672612419124243D+00
      S2C(21,10) =-0.4082482904638632D+00
      S2C(26,10) = 0.8908708063747472D-01
      S2C( 6,11) = 0.8908708063747472D-01
      S2C(17,11) =-0.4082482904638631D+00
      S2C(24,11) = 0.2672612419124243D+00
      S2C( 7,12) =-0.4652421051992352D-01
      S2C(18,12) = 0.2672612419124243D+00
      S2C(25,12) =-0.2672612419124243D+00
      S2C(28,12) = 0.4652421051992350D-01
      S2C(13,13) = 0.1543033499620918D+00
      S2C(22,13) =-0.3162277660168381D+00
      S2C(27,13) = 0.1543033499620918D+00
    case default
      write(*,"(/,' CoefS2C_GAUSS has been programmed only for Lq = 0 ~ 6!')")
      call estop(1)
  end select

!  do i=1,nc
!    write(77,"(100f10.5)")S2C(i,:)
!  end do

  return
end subroutine CoefS2C_GAUSS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Transform MO coefficients from Cartesian to spherical dimension.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MOCar2Sph(natom,mxla,nshlla,ncbas,nsbas,nmo,cmo,smo,C2S)
 use Constant, only : maxlq, One, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension  :: mxla(natom), nshlla(0:maxlq,natom), cmo(ncbas,nmo), smo(nsbas,nmo), C2S(*)

 do imo=1,nmo
     ioffs=0
     ioffc=0
     do iatom=1,natom
        icoef=0
        do li=0,mxla(iatom)
           nprimi=nshlla(li,iatom)
           nsph=li+li+1
           ncar=(li+1)*(li+2)/2
           do ish=1,nprimi
              call MatMult(1,1,ncar,nsph,One,Zero,cmo(ioffc+1,imo),C2S(icoef+1),smo(ioffs+1,imo))
              ioffs=ioffs+nsph
              ioffc=ioffc+ncar
           end do ! ish
           icoef=((li+1)*(li+2)*(li+3)*(li*3+2))/12
        enddo ! li
     enddo ! iatom
 enddo ! imo

 return
end

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!
!! Transform MO coefficients from Cartesian to spherical dimension. This is a 2-component version of subroutine MOCar2Sph.
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine MOCar2Sph2C(natom,mxla,nshlla,ncbas,nsbas,nmo,cmo,smo,C2S)
! use Constant, only : maxlq, One, Zero
! implicit real(kind=8) (a-h,o-z)
! dimension  :: mxla(natom), nshlla(0:maxlq,natom), cmo(4,ncbas,nmo), smo(4,nsbas,nmo), C2S(*)
! allocatable:: src(:,:), srs(:,:)
!
! allocate(src(ncbas,nmo),srs(nsbas,nmo))
!
! do idim = 1, 4
!   do i = 1, nmo
!     do j = 1, ncbas
!       src(j,i) = cmo(idim,j,i)
!     end do
!   end do
!
!   srs = Zero
!
!   do imo=1,nmo
!       ioffs=0
!       ioffc=0
!       do iatom=1,natom
!          icoef=0
!          do li=0,mxla(iatom)
!             nprimi=nshlla(li,iatom)
!             nsph=li+li+1
!             ncar=(li+1)*(li+2)/2
!             do ish=1,nprimi
!                call MatMult(1,1,ncar,nsph,One,Zero,src(ioffc+1,imo),C2S(icoef+1),srs(ioffs+1,imo))
!                ioffs=ioffs+nsph
!                ioffc=ioffc+ncar
!             end do ! ish
!             icoef=((li+1)*(li+2)*(li+3)*(li*3+2))/12
!          enddo ! li
!       enddo ! iatom
!   enddo ! imo
!
!   do i = 1, nmo
!     do j = 1, nsbas
!       smo(idim,j,i) = srs(j,i)
!     end do
!   end do
! end do
!
! deallocate(src,srs)
!
! return
!end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Transform a matrix from Cartesian to spherical dimension (IOP=0) or spherical to Cartesian dimension (IOP/=0).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Car2Sph(iop,natom,mxla,nshlla,NBAS,np,DC,DS,C2S,SCR1,SCR2)
 use Constant, only : maxlq
 implicit real(kind=8) (a-h,o-z)
 dimension  :: mxla(natom),nshlla(0:maxlq,natom), DC(NBAS,NBAS),DS(np,np),C2S(*),SCR1(*),SCR2(*)

 ioffs=0
 ioffc=0
 do iatom=1,natom
    icoef=0
    do li=0,mxla(iatom)
       nprimi=nshlla(li,iatom)
       isph=li+li+1
       icar=(li+1)*(li+2)/2
       do ish=1,nprimi

          joffs=0
          joffc=0
          do jatom=1,natom
             jcoef=0
             do lj=0,mxla(jatom)
                nprimj=nshlla(lj,jatom)
                jsph=lj+lj+1
                jcar=(lj+1)*(lj+2)/2
                do jsh=1,nprimj
                   if(iop == 0) then
                     call Car2SphSl(icar,isph,C2S(icoef+1),jcar,jsph,C2S(jcoef+1),NBAS,np,DC,DS,ioffs,ioffc,joffs,joffc,SCR1,SCR2)
                   else
                     call Sph2CarSl(icar,isph,C2S(icoef+1),jcar,jsph,C2S(jcoef+1),NBAS,np,DC,DS,ioffs,ioffc,joffs,joffc,SCR1,SCR2)
                   end if
                   joffs=joffs+jsph
                   joffc=joffc+jcar
                end do ! jsh
                jcoef=((lj+1)*(lj+2)*(lj+3)*(lj*3+2))/12
             enddo ! lj
          enddo ! jatom

          ioffs=ioffs+isph
          ioffc=ioffc+icar
       end do ! ish
       icoef=((li+1)*(li+2)*(li+3)*(li*3+2))/12
    enddo ! li
 enddo ! iatom

! do i=1,NBAS
!   write(77,"(100f12.6)")DC(:,i)
! end do
! write(77,*)
! do i=1,np
!   write(77,"(100f12.6)")DS(:,i)
! end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! transform a matrix from Cartesian to spherical dimension for a give block.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Car2SphSl(ncl,nsl,C2SL,ncr,nsr,C2SR,nc,np,DC,DS,ioffs,ioffc,joffs,joffc,SCR1,SCR2)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: C2SL(ncl,nsl),C2SR(ncr,nsr),DC(nc,nc),DS(np,np),SCR1(*),SCR2(*)
 ! SCR1(ncl,ncr),SCR2(nsl,nsr),SCR3(*)

 !SCR1 = DC(ioffc+1:ioffc+ncl,joffc+1:joffc+ncr)
 call xcopy(nc,ncl,ncr,ioffc,0,ncl,DC(1,joffc+1),SCR1)

 if(ncl == nsl) then
   call acopy(ncl*ncr,SCR1,SCR2)
 else
   call dgemm("t","n",nsl,ncr,ncl,One,C2SL,ncl,SCR1,ncl,Zero,SCR2,nsl)
 end if

 if(ncr == nsr) then
   call acopy(nsl*nsr,SCR2,SCR1)
 else
   call dgemm("n","n",nsl,nsr,ncr,One,SCR2,nsl,C2SR,ncr,Zero,SCR1,nsl)
 end if

 !DS(ioffs+1:ioffs+nsl,joffs+1:joffs+nsr) = SCR2
 call xcopy(nsl,np,nsr,0,ioffs,nsl,SCR1,DS(1,joffs+1))

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! transform a matrix from spherical to Cartesian dimension for a give block.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Sph2CarSl(ncl,nsl,C2SL,ncr,nsr,C2SR,nc,np,DC,DS,ioffs,ioffc,joffs,joffc,SCR1,SCR2)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: C2SL(ncl,nsl),C2SR(ncr,nsr),DC(nc,nc),DS(np,np),SCR1(*),SCR2(*)
 ! SCR1(ncl,ncr),SCR2(nsl,nsr),SCR3(*)

 !SCR2 = DS(ioffs+1:ioffs+nsl,joffs+1:joffs+nsr)
 call xcopy(np,nsl,nsr,ioffs,0,nsl,DS(1,joffs+1),SCR1)

 if(ncl == nsl) then
   call acopy(nsl*nsr,SCR1,SCR2)
 else
   call dgemm("n","n",ncl,nsr,nsl,One,C2SL,ncl,SCR1,nsl,Zero,SCR2,ncl)
 end if

 if(ncr == nsr) then
   call acopy(ncl*ncr,SCR2,SCR1)
 else
   call dgemm("n","t",ncl,ncr,nsr,One,SCR2,ncl,C2SR,ncr,Zero,SCR1,ncl)
 end if

 !DC(ioffc+1:ioffc+ncl,joffc+1:joffc+ncr) = SCR1
 call xcopy(ncl,nc,ncr,0,ioffc,ncl,SCR1,DC(1,joffc+1))

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read fchk: 1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_fchk_1(iout,imod,ctmp,natom,nmo)
 implicit real(kind=8) (a-h,o-z)
 character*(*)     :: ctmp

 rewind(imod)

 ! the third line should be "Number of atoms ..."
 read(imod,"(//,a61)",err=1000,end=1000) ctmp
 if(index(ctmp,'Number of atoms                            I') == 0) goto 1000
 read(ctmp(45:61),*) natom

 nbs  = 0
 nmo  = 0
 do while(.true.)
   read(imod,"(a61)",err=1100,end=1100) ctmp
   if(index(ctmp,'Number of basis functions                  I') > 0) then
     read(ctmp(45:61),*) nbs
   else if(index(ctmp,'Alpha Orbital Energies                     R   N=') > 0) then
     read(ctmp(50:61),*) nmoa
     nmo = nmoa
   else if(index(ctmp,'Beta Orbital Energies                      R   N=') > 0) then
     read(ctmp(50:61),*) nmob
     nmo = nmo + nmob
   else if(index(ctmp,'Alpha MO coefficients                      R   N=') > 0) then
     exit
   end if
 end do
 if(nbs /= nmo .and. nbs*2 /= nmo) &
   write(iout,"(' ### Warning in rd_fchk_1: #bs=',i5,' and #mo',i5,' do not match.')") nbs, nmo

 return

 1000 write(iout,"(' ### Error: this is not a standardized fchk file!')")
 call estop(1)
 1100 write(iout,"(' ### Error in rd_fchk_1.')")
 call estop(1)
 return
end subroutine rd_fchk_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read fchk: 2
!
! Count the numbers of contracted Cartesian (ncbas), contracted spherical (nsbas), primitive Cartesian (npcar), and primitive
! spherical basis functions (npsph).
! The MOs may be in Cartesian (ispher=0) or spherical (ispher=1) functions, although spherical functions should always be used
! in calculations.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_fchk_2(iout,imod,igto,irelc,ctmp, natom, nmo, maxtyp, ispher, ncbas, nsbas, npcar, npsph, nexp, ncoeff)
 use Constant, only : maxlq_G, maxcn, Zero
 implicit real(kind=8) (a-h,o-z)
 character*(*)     :: ctmp
 allocatable       :: ityp(:),npsh(:),imap(:),exp1(:),con1(:),extmp(:,:),cotmp(:,:,:),ncsh(:),iptr(:,:,:),idel(:)

 rewind(imod)

 itag  = 0
 maxtyp= 0
 ncbas = 0
 nsbas = 0
 npcar = 0
 npsph = 0
 nexp  = 0
 ncoeff= 0
 ispher= 0
 do while(.true.)
   read(imod,"(a80)",err=1100) ctmp
   if     (index(ctmp,'Highest angular momentum                   I') > 0) then
     read(ctmp(45:61),*) maxtyp
     if(maxtyp > maxlq_G) then
       write(iout,"(' ### Error: please increase maxlq_G!')")
       call estop(1)
     end if
     itag = itag + 1

   else if(index(ctmp,'Largest degree of contraction              I') > 0) then
     read(ctmp(45:61),*) ncon
     if(ncon < 1 .or. ncon > maxcn) then
       write(iout,"(' ### Error: unreasonable NCon',i4)") ncon
       call estop(1)
     end if
     itag = itag + 1

   else if(index(ctmp,'Number of contracted shells                I') > 0) then
     read(ctmp(45:61),*) nshl
     if(nshl < 1) then
       write(iout,"(' ### Error: unreasonable NShl',i4)") nshl
       call estop(1)
     end if
     allocate(ityp(nshl),npsh(nshl),imap(nshl))
     itag = itag + 1

   else if(index(ctmp,'Shell types                                I   N=') > 0) then
     ! 0: S; 1: P; -1: SP; -L / +L (L > 1): spherical / Cartesian L
     read(imod,*,err=1100) (ityp(i), i=1,nshl)
     do i = 1, nshl
       if(ityp(i) == -1) then
         write(iout,"(' ### Error: the SP shell of basis function is not supported!')")
         call estop(1)
       end if
     end do
     itag = itag + 1

   else if(index(ctmp,'Number of primitives per shell             I   N=') > 0) then
     read(imod,*,err=1100) (npsh(i), i=1,nshl)
     ndat = sum( npsh(1:nshl) )
     allocate(exp1(ndat),con1(ndat))
     itag = itag + 1

   else if(index(ctmp,'Shell to atom map                          I   N=') > 0) then
     read(imod,*,err=1100) (imap(i), i=1,nshl)
     do i = 2, nshl
       if(imap(i) < imap(i-1)) then
         write(iout,"(' ### Error: basis functions of each atom must be saved continuously.')")
         call estop(1)
       end if
     end do
     itag = itag + 1

   else if(index(ctmp,'Primitive exponents                        R   N=') > 0) then
     read(ctmp(50:61),*) i
     if(ndat /= i) then
       write(iout,"(' ### Error: unreasonable number of exponents',i4)") i
       call estop(1)
     end if
     read(imod,*,err=1100) (exp1(i), i=1,ndat)
     itag = itag + 1

   else if(index(ctmp,'Contraction coefficients                   R   N=') > 0) then
     read(ctmp(50:61),*) i
     if(ndat /= i) then
       write(iout,"(' ### Error: unreasonable number of coefficients',i4)") i
       call estop(1)
     end if
     read(imod,*,err=1100) (con1(i), i=1,ndat)
     itag = itag + 1

   end if

   if(itag == 8) exit
 end do

 maxcn2=maxcn*maxcn
 allocate(extmp(maxcn2,0:maxlq_G),cotmp(maxcn2,maxcn,0:maxlq_G),ncsh(0:maxlq_G),iptr(2,maxcn,0:maxlq_G),idel(maxcn2), stat=ierr)
   if(ierr /= 0) then
     write(iout,"(' ### Insufficient Memory in rd_fchk_2!')")
     call estop(1)
   end if
 do ia = 1, natom
   mlqa = 0
   ncsh = 0
   iptr = 0
   ilst = 0
   extmp= Zero
   cotmp= Zero
   do is = 1, nshl
     ilst = ilst + npsh(is)
     if(imap(is) /= ia) then
       if(ncsh(1) == 0) then
         cycle
       else
         exit
       end if
     end if
     lq = abs(ityp(is))
     if(mlqa < lq) then
       mlqa = lq
     else if(mlqa > lq) then
       write(iout,"(' ### Error: basis functions must be sorted in ascending order according to LQ.')")
       call estop(1)
     end if

     ncsh(lq) = ncsh(lq) + 1
     if(ncsh(lq) == 1) then
       iptr(1,ncsh(lq),lq) = 1
       iptr(2,ncsh(lq),lq) = npsh(is)
     else
       iptr(1,ncsh(lq),lq) = iptr(2,ncsh(lq)-1,lq)+1
       iptr(2,ncsh(lq),lq) = iptr(2,ncsh(lq)-1,lq)+npsh(is)
     end if

     extmp(iptr(1,ncsh(lq),lq):iptr(2,ncsh(lq),lq),lq) = exp1(ilst+1-npsh(is):ilst)
     cotmp(iptr(1,ncsh(lq),lq):iptr(2,ncsh(lq),lq),ncsh(lq),lq) = con1(ilst+1-npsh(is):ilst)
   end do  ! is

   ! remove duplicate primitive functions
   do lq = 0, mlqa
     if(ncsh(lq) < 1) cycle

     ngto = iptr(2,ncsh(lq),lq)
     call Rmdupfun(maxcn2, ngto, ncsh(lq), idel, iptr(1,1,lq), extmp(1,lq), cotmp(1,1,lq))
     ngto = iptr(2,ncsh(lq),lq)
     call checkbas(iout,ia,lq,maxcn2,ngto,ncsh(lq),extmp(1,lq), cotmp(1,1,lq))
     ncbas = ncbas + (lq+1)*(lq+2)*ncsh(lq)/2
     nsbas = nsbas + (lq+lq+1)*ncsh(lq)
     npcar = npcar + (lq+1)*(lq+2)*ngto/2
     npsph = npsph + (lq+lq+1)*ngto
     nexp  = nexp  + ngto
     ncoeff= ncoeff+ ngto*ncsh(lq)
   end do

   ! save basis functions to scratch file igto in a MOLCAS-like format
   write(igto,"(2i5)")ia, mlqa
   write(igto,"(4x,10i5)") (iptr(2,max(ncsh(lq),1),lq), lq=0,mlqa)
   write(igto,"(4x,10i5)") (ncsh(lq), lq=0,mlqa)
   do lq = 0, mlqa
     if(ncsh(lq) < 1) cycle

     ngto = iptr(2,ncsh(lq),lq)
     ! write(igto,"(4x,2i5)") ngto, ncsh(lq)
     write(igto,"(d20.12)") extmp(1:ngto,lq)
     do j = 1,ngto
       write(igto,"(100d20.12)") cotmp(j,1:ncsh(lq),lq)
     end do
   end do

 end do  ! ia

 ! RHF & ROHF: nbas = nmo; UHF & GHF: nbas*2 = nmo
 if(ncbas*irelc == nmo .and. maxtyp > 1) then
   ispher=0
 else if(irelc == 1 .and. ncbas*2 == nmo .and. maxtyp > 1) then    ! UHF, UKS
   ispher=0
 ! else if(nsbas*irelc == nmo) then  ! Gaussian may eliminate some highest-lying virtual MOs.
 else if(nsbas*irelc >= nmo) then
   ispher=1
 else if(irelc == 1 .and. nsbas*2 == nmo) then    ! UHF, UKS
   ispher=1
 else
   write(iout,*) nsbas, irelc, nmo
   write(iout,"(' ### Error: Cartesian and spherical functions cannot be mixed together.')")
   call estop(1)
 end if

 deallocate(ityp,npsh,imap,exp1,con1,extmp,cotmp,iptr,ncsh,idel)

 return

 1100 write(iout,"(' ### Error in rd_fchk_2 (1).')")
 call estop(1)
 return
end subroutine rd_fchk_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read fchk 3. Cartesian coordinates
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_fchk_3(iout,imod,ctmp,natom,ifecp,iza,za,xyz)
 implicit real(kind=8) (a-h,o-z)
 character*(*)     :: ctmp
 dimension         :: iza(natom), za(natom), xyz(3*natom)
 logical           :: ifecp

 rewind(imod)
 do while(.true.)
   read(imod,"(a80)",err=1100,end=1100) ctmp
   if     (index(ctmp,'Atomic numbers                             I   N=') > 0) then
     read(imod,*,err=1200,end=1200) (iza(i), i=1,natom)
   else if(index(ctmp,'Nuclear charges                            R   N=') > 0) then
     read(imod,*,err=1200,end=1200) (za(i), i=1,natom)
   else if(index(ctmp,'Current cartesian coordinates              R   N=') > 0) then
     read(imod,*,err=1200,end=1200) (xyz(i), i=1,3*natom)
     exit
   end if
 end do
 do i = 1, natom
   if(nint(za(i)) /= iza(i)) then
     !write(iout,"(' ### Error: PP cannot be used for all-electron relativistic calculation.')")
     !call estop(1)
     ifecp = .true.
     exit
   end if
 end do
 call prtcoord(iout,natom,iza,xyz)

 return

 1100 write(iout,"(' ### Error in rd_fchk_3 (1).')")
 call estop(1)
 1200 write(iout,"(' ### Error in rd_fchk_3 (2).')")
 call estop(1)
 return
end subroutine rd_fchk_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read fchk: 4. read MOs.
! In the case of 2c calculation, nbas <= 4*nbas for alpha MO coefficients
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_fchk_4(iout,imod,irelc,nbas,nmo,occtot,occ,ene,ispn,symm,ctmp,cmo,den)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 character*10      :: symm
 character*(*)     :: ctmp
 dimension         :: occ(nmo), ene(nmo), ispn(nmo), symm(nmo), cmo(nbas,*), den(*)
 logical           :: uscf

 uscf = .false.
 rewind(imod)

 do while(.true.)
   read(imod,"(a80)",err=1100,end=1100) ctmp
   if     (index(ctmp,'Number of alpha electrons                  I') > 0) then
     read(ctmp(45:61),*) na
   else if(index(ctmp,'Number of beta electrons                   I') > 0) then
     read(ctmp(45:61),*) nb
   else if(index(ctmp,'Alpha Orbital Energies                     R   N=') > 0) then
     read(ctmp(50:61),*) nmoa
     read(imod,*,err=1200,end=1200) (ene(i), i=1,nmoa)
   else if(index(ctmp,'Beta Orbital Energies                      R   N=') > 0) then
     read(ctmp(50:61),*) nmob
     if(nmob /= nmoa) then
       write(iout,"(' ### Error: #MO_alpha .ne. #MO_beta.')")
       call estop(1)
     end if
     read(imod,*,err=1200,end=1200) (ene(i), i=nmoa+1,nmoa+nmob)
     uscf = .true.
   else if(index(ctmp,'Alpha MO coefficients                      R   N=') > 0) then
     read(ctmp(50:61),*) i
     if(i /= nbas*nmoa) goto 1300
     read(imod,*,err=1200,end=1200) ((cmo(i,j), i=1,nbas), j=1,nmoa)
   else if(index(ctmp,'Beta MO coefficients                       R   N=') > 0) then
     read(ctmp(50:61),*) i
     if(i /= nbas*nmoa) goto 1300
     read(imod,*,err=1200,end=1200) ((cmo(i,j), i=1,nbas), j=nmoa+1,nmoa+nmob)
   else if(index(ctmp,'Total SCF Density                          R   N=') > 0) then
     read(ctmp(50:61),*) i
     if(irelc == 1) then
       if(i /= nbas*(nbas+1)/2) goto 1400
     else
       ! nbas*4 --> nbas
       if(i /= nbas*(nbas/2+1)/2) goto 1400
     end if
     read(imod,*,err=1200,end=1200) (den(j), j=1,i)
     exit
   end if
 end do

 ispn(1:nmoa) = 1
 occ = Zero
 if(irelc == 1) then  ! 1c
   symm(1:nmo) = 'A'
   occ(1:na) = One
   if(uscf) then
     if(nmoa+nmob /= nmo) goto 1500
     occ(nmoa+1:nmoa+nb) = One
     ispn(nmoa+1:nmo) = 2
   else
     if(nmoa /= nmo) goto 1500
     do i = 1, nb
       occ(i) = occ(i) + One
     end do
   end if
 else if(irelc == 2) then  ! 2c
   symm(1:nmo) = 'E_1/2'
   if(nmoa /= nmo) goto 1500
   occ(1:na+nb) = One
 end if
 occtot = sum(occ)

 return

 1100 write(iout,"(' ### Error in rd_fchk_4 (1).')")
 call estop(1)
 1200 write(iout,"(' ### Error in rd_fchk_4 (2).')")
 call estop(1)
 1300 write(iout,"(' ### Error in rd_fchk_4 (3). The number of MO coefficients is wrong!')")
 call estop(1)
 1400 write(iout,"(' ### Error in rd_fchk_4 (4). The number of elements in density matrix is wrong!')")
 call estop(1)
 1500 write(iout,"(' ### Error in rd_fchk_4 (5). The number of MOs is wrong!')")
 call estop(1)
 return
end subroutine rd_fchk_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read post-HF density matrix from fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_fchk_den(iout,imod,iden,irelc,nbas,ctmp,den)
 implicit real(kind=8) (a-h,o-z)
 character*(*)     :: ctmp
 dimension         :: den(*)

 if(iden /= 2) then
   if(iden == 1) write(iout,"(/,' The SCF density matrix will be used.')")
   return
 else if(irelc == 2) then
   write(iout,"(/,' The two-component post-HF calculation is not supported.',/)")
   call estop(1)
 end if

 ! rewind(imod)

 do while(.true.)
   read(imod,"(a80)",err=1000,end=1000) ctmp
   if     (index(ctmp,'Total MP2 Density                          R   N=') > 0) then
     write(iout,"(/,' The MP2 density matrix will be used.')")
     exit
   else if(index(ctmp,'Total MP3 Density                          R   N=') > 0) then
     write(iout,"(/,' The MP3 density matrix will be used.')")
     exit
   else if(index(ctmp,'Total MP4 Density                          R   N=') > 0) then
     write(iout,"(/,' The MP4(SDQ) density matrix will be used.')")
     exit
   else if(index(ctmp,'Total CI Density                           R   N=') > 0) then
     write(iout,"(/,' The CI density matrix will be used.')")
     exit
   else if(index(ctmp,'Total CC Density                           R   N=') > 0) then
     write(iout,"(/,' The CC density matrix will be used.')")
     exit
   end if
 end do

 read(ctmp(50:61),*) i
 if(i /= nbas*(nbas+1)/2) goto 1100
 read(imod,*,err=1200,end=1200) (den(j), i=1,i)

 return

 1000 write(iout,"(/,' No post-HF density matrix found. The SCF density matrix will be used.',/)")
 return

 1100 write(iout,"(' ### Error in rd_fchk_den (1). The number of elements in density matrix is wrong!')")
 call estop(1)
 1200 write(iout,"(' ### Error in rd_fchk_den (2).')")
 call estop(1)
 return
end subroutine rd_fchk_den

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read molden: 1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_molden_1(iout,imod,ctmp,natom,nmo)
 implicit real(kind=8) (a-h,o-z)
 character*(*)     :: ctmp

 rewind(imod)

 ! the first line should be [Molden Format]
 read(imod,"(a200)",err=1000,end=1000) ctmp
 if(index(ctmp,'[Molden Format]') == 0) goto 1000
 ! the second line should be [Porgram] BDF or [Title]
 read(imod,"(a200)",err=1000,end=1000) ctmp
 if(index(ctmp,'[Porgram] BDF') /= 0) then
   ! do nothing
 else if(index(ctmp,'[Title]') /= 0) then
   ! there should be "Molden2AIM" or "Multiwfn" in the third or fifth line
   do while(.true.)
     read(imod,"(a200)",err=1000,end=1000) ctmp
     if(index(ctmp,'Host name:') > 0 .or. index(ctmp,'File name:') > 0) cycle
     call charl2u(ctmp)
     if(index(ctmp,'MOLDEN2AIM') == 0 .and. index(ctmp,'MULTIWFN') == 0) goto 1000
     exit
   end do
 else
   goto 1000
 end if

 natom = 0
 nmo   = 0
 do while(.true.)
   read(imod,"(a200)",err=1100,end=1100) ctmp
   call charl2u(ctmp)
   if(index(ctmp,'[PSEUDO]') > 0 .or. index(ctmp,'[CORE]') > 0) then
     write(iout,"(' ### Error: PP cannot be used for all-electron relativistic calculation.')")
     call estop(1)
   else if(index(ctmp,'[ATOMS]') > 0) then
     ! count the number of atoms
     do while(.true.)
       read(imod,"(a200)",err=1200,end=1200) ctmp
       if(len_trim(ctmp) == 0) exit
       if(index(ctmp,'[') > 0 .and. index(ctmp,']') > 0) exit
       natom = natom + 1
     end do
   else if(index(ctmp,'[MO]') > 0) then
     ! count the number of MOs
     do while(.true.)
       read(imod,"(a200)",err=1300,end=2000) ctmp
       if(len_trim(ctmp) == 0) exit
       if(index(ctmp,'[') > 0 .and. index(ctmp,']') > 0) exit
       call charl2u(ctmp)
       if(index(ctmp,'OCCUP') > 0) nmo = nmo + 1
     end do
     2000  continue
   end if
   if(natom > 0 .and. nmo > 0) exit
 end do

 return

 1000 write(iout,"(' ### Error: this is not a standardized molden file!')")
 call estop(1)
 1100 write(iout,"(' ### Error in rd_molden_1 (1).')")
 call estop(1)
 1200 write(iout,"(' ### Error in rd_molden_1 (2).')")
 call estop(1)
 1300 write(iout,"(' ### Error in rd_molden_1 (3).')")
 call estop(1)
 return
end subroutine rd_molden_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read molden: 2
!
! Count the numbers of contracted Cartesian (ncbas), contracted spherical (nsbas), primitive Cartesian (npcar), and primitive
! spherical basis functions (npsph).
! The MOs may be in Cartesian (ispher=0) or spherical (ispher=1) functions, although spherical functions should always be used
! in calculations.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_molden_2(iout,imod,igto,ctmp,clqdat, natom, maxtyp, ispher, ncbas, nsbas, npcar, npsph, nexp, ncoeff)
 use Constant, only : maxlq, maxcn, Zero
 implicit real(kind=8) (a-h,o-z)
 character*(*)     :: ctmp
 character*1       :: clqdat(0:maxlq), clq
 logical           :: ifok
 allocatable       :: extmp(:,:),cotmp(:,:,:),iptr(:,:,:),ncsh(:),idel(:)

 clqdat(0:maxlq) = (/'S','P','D','F','G','H'/)

 rewind(imod)

 ifok=.false.
 maxcn2=maxcn*maxcn
 maxtyp= 0
 ncbas = 0
 nsbas = 0
 npcar = 0
 npsph = 0
 nexp  = 0
 ncoeff= 0
 ispher= 0
 do while(.true.)
   read(imod,"(a200)",err=1100) ctmp
   call charl2u(ctmp)
   if(index(ctmp,'[PSEUDO]') > 0 .or. index(ctmp,'[CORE]') > 0) then
     write(iout,"(' ### Error: PP cannot be used for all-electron relativistic calculation.')")
     call estop(1)
   else if(index(ctmp,'[GTO]') > 0) then
     allocate(extmp(maxcn2,0:maxlq),cotmp(maxcn2,maxcn,0:maxlq),iptr(2,maxcn,0:maxlq),ncsh(0:maxlq),idel(maxcn2), stat=ierr)
       if(ierr /= 0) then
         write(iout,"(' ### Insufficient Memory in rd_molden_2!')")
         call estop(1)
       end if

     do ia = 1, natom
       mlqa = 0
       ncsh = 0
       iptr = 0
       extmp= Zero
       cotmp= Zero
       read(imod,*,err=1200,end=1200)
       do while(.true.)
         read(imod,"(a200)",err=1200,end=1200) ctmp
         if(len_trim(ctmp) == 0) exit
         read(ctmp,*) clq, ncon
         call charl2u(clq)
         do lq = 0, maxlq
           if(clq == clqdat(lq)) exit
           if(lq == maxlq) then
             write(iout,"(' ### Error: unknown angular momentum ',a1)") clq
             call estop(1)
           end if
         end do
         if(mlqa < lq) then
           mlqa = lq
         else if(mlqa > lq) then
           write(iout,"(' ### Error: basis functions must be sorted in ascending order according to LQ.')")
           call estop(1)
         end if
         if(ncon < 1 .or. ncon > maxcn) then
           write(iout,"(' ### Error: unreasonable NCon',i4)") ncon
           call estop(1)
         end if

         ncsh(lq) = ncsh(lq) + 1
         if(ncsh(lq) == 1) then
           iptr(1,ncsh(lq),lq) = 1
           iptr(2,ncsh(lq),lq) = ncon
         else
           iptr(1,ncsh(lq),lq) = iptr(2,ncsh(lq)-1,lq)+1
           iptr(2,ncsh(lq),lq) = iptr(2,ncsh(lq)-1,lq)+ncon
         end if

         do i = 0, ncon-1
           read(imod,*,err=1200,end=1200) extmp(iptr(1,ncsh(lq),lq)+i,lq), cotmp(iptr(1,ncsh(lq),lq)+i,ncsh(lq),lq)
         end do
       end do
       maxtyp = max(maxtyp,mlqa)
       !! <<< debug
       !write(77,"('ia=',i5,', mlqa=',i1)")ia, mlqa
       !do i = 0, mlqa
       !  ngto = iptr(2,ncsh(i),i)
       !  write(77,"(4x,'LQ=',i2,', ngto=',i5,', ncsh=',i5)") i, ngto, ncsh(i)
       !  do j = 1,ngto
       !    write(77,"(4x,d16.8,100d16.6)") extmp(j,i), cotmp(j,1:ncsh(i),i)
       !  end do
       !end do
       !! >>>

       ! remove duplicate primitive functions
       do i = 0, mlqa
         if(ncsh(i) < 1) cycle

         ngto = iptr(2,ncsh(i),i)
         call Rmdupfun(maxcn2, ngto, ncsh(i), idel, iptr(1,1,i), extmp(1,i), cotmp(1,1,i))
         ngto = iptr(2,ncsh(i),i)
         call checkbas(iout,ia,i,maxcn2,ngto,ncsh(i),extmp(1,i), cotmp(1,1,i))
         ncbas = ncbas + (i+1)*(i+2)*ncsh(i)/2
         nsbas = nsbas + (i+i+1)*ncsh(i)
         npcar = npcar + (i+1)*(i+2)*ngto/2
         npsph = npsph + (i+i+1)*ngto
         nexp  = nexp  + ngto
         ncoeff= ncoeff+ ngto*ncsh(i)
       end do

       ! save basis functions to scratch file igto in a MOLCAS-like format
       write(igto,"(2i5)")ia, mlqa
       write(igto,"(4x,10i5)") (iptr(2,max(ncsh(i),1),i), i=0,mlqa)
       write(igto,"(4x,10i5)") (ncsh(i), i=0,mlqa)
       do i = 0, mlqa
         if(ncsh(i) < 1) cycle

         ngto = iptr(2,ncsh(i),i)
         ! write(igto,"(4x,2i5)") ngto, ncsh(i)
         write(igto,"(d20.12)") extmp(1:ngto,i)
         do j = 1,ngto
           write(igto,"(100d20.12)") cotmp(j,1:ncsh(i),i)
         end do
       end do

     end do  ! ia
   else if(index(ctmp,'[MO]') > 0) then
     ! check whether Cartesian basis functions in MOLDEN file
     do while(.true.)
       read(imod,"(a200)",err=1300,end=1300) ctmp
       if(len_trim(ctmp) == 0) exit
       if(index(ctmp,'[') > 0 .and. index(ctmp,']') > 0) exit
     
       call charl2u(ctmp)
       if(index(ctmp,'OCCUP') > 0) exit
     end do
     nmo = 0
     do while(.true.)
       read(imod,"(a200)",err=1300,end=2000) ctmp
       if(len_trim(ctmp) == 0 .or. index(ctmp,'=') > 0) exit
       if(index(ctmp,'[') > 0 .and. index(ctmp,']') > 0) exit
       nmo = nmo + 1
       read(ctmp,*,err=1400,end=1400) i, atmp
     end do
     2000  continue
     if(ncbas == nmo .and. maxtyp > 1) then
       ifok = .true.
       ispher=0
     else if(nsbas == nmo) then
       ifok = .true.
       ispher=1
     else
       write(iout,"(' ### Error: MOs should be in the dimension of either Cartesian or spherical functions.')")
       call estop(1)
     end if
   end if
   if(ifok) exit
 end do

 deallocate(extmp,cotmp,iptr,ncsh,idel)

 return

 1100 write(iout,"(' ### Error in rd_molden_2 (1).')")
 call estop(1)
 1200 write(iout,"(' ### Error in rd_molden_2 (2).')")
 call estop(1)
 1300 write(iout,"(' ### Error in rd_molden_2 (3).')")
 call estop(1)
 1400 write(iout,"(' ### Error in rd_molden_2 (4).')")
 call estop(1)
 return
end subroutine rd_molden_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read molden: 3. Cartesian coordinates
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_molden_3(iout,imod,ctmp,natom,iza,za,xyz)
 use Constant, only : One, au2ang
 implicit real(kind=8) (a-h,o-z)
 character*(*)     :: ctmp
 logical           :: ifangs
 dimension         :: iza(natom), za(natom), xyz(3,natom)

 rewind(imod)
 do while(.true.)
   read(imod,"(a200)",err=1100,end=1100) ctmp
   call charl2u(ctmp)
   if(index(ctmp,'[ATOMS]') > 0) then
     ifangs = (index(ctmp,'ANGS') > 0)
     exit
   end if
 end do

 do i = 1, natom
   read(imod,*,err=1200,end=1200) ctmp, j, iza(i), xyz(:,i)
   call ElemZA(0,ctmp,j)
   if(iza(i) /= j) then
     write(iout,"(' ### Error: PP cannot be used for all-electron relativistic calculation.')")
     call estop(1)
   end if
   za(i) = dble(iza(i))
 end do
 if(ifangs) xyz = (One/au2ang) * xyz

 call prtcoord(iout,natom,iza,xyz)

 return

 1100 write(iout,"(' ### Error in rd_molden_3 (1).')")
 call estop(1)
 1200 write(iout,"(' ### Error in rd_molden_3 (2).')")
 call estop(1)
 return
end subroutine rd_molden_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read molden: 4
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_molden_4(iout,imod,nbas,nmo,occtot,occ,ene,ispn,symm,ctmp,cmo)
 implicit real(kind=8) (a-h,o-z)
 character*10      :: symm,stmp
 character*(*)     :: ctmp
 dimension         :: occ(nmo), ene(nmo), ispn(nmo), symm(nmo), cmo(nbas,*)

 !============================================================================== read MOs
 rewind(imod)
 do while(.true.)
   read(imod,"(a200)",err=1100,end=1100) ctmp
   call charl2u(ctmp)
   if(index(ctmp,'[MO]') > 0) exit
 end do

 stmp = 'A'
 etmp = 0.0d0
 itmp = 1
 imo  = 0
 do while(.true.)
   read(imod,"(a200)",err=1100,end=2000) ctmp
   if(len_trim(ctmp) == 0) exit
   if(index(ctmp,'[') > 0 .and. index(ctmp,']') > 0) exit

   call charl2u(ctmp)
   if(index(ctmp,'SYM') > 0) then
     ieq = index(ctmp,'=')
     read(ctmp(ieq+1:),*) stmp
   else if(index(ctmp,'ENE') > 0) then
     ieq = index(ctmp,'=')
     if(index(ctmp(ieq+1:),'***') == 0) read(ctmp(ieq+1:),*) etmp
   else if(index(ctmp,'SPIN') > 0) then
     ieq = index(ctmp,'=')
     if(index(ctmp(ieq+1:),'BETA') > 0) itmp = 2
   else if(index(ctmp,'OCCUP') > 0) then
     imo = imo + 1
     ieq = index(ctmp,'=')
     read(ctmp(ieq+1:),*) occ(imo)
     ene(imo) = etmp
     ispn(imo) = itmp
     symm(imo) = stmp
     stmp = 'A'
     etmp = 0.0d0
     itmp = 1
     do i = 1, nbas
       read(imod,*,err=1200,end=1200) j, cmo(i,imo)
     end do
   end if
 end do

 2000  if(imo /= nmo) goto 1300
 occtot = sum(occ)
 return

 1100 write(iout,"(' ### Error in rd_molden_4 (1).')")
 call estop(1)
 1200 write(iout,"(' ### Error in rd_molden_4 (2).')")
 call estop(1)
 1300 write(iout,"(' ### Error in rd_molden_4 (3).')")
 call estop(1)
 return
end subroutine rd_molden_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read basis functions from igto in a MOLCAS-like format
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rd_basfun(iout,igto,natom,mxla,nshlla,ialbs,ialbc,iptexp,iptcon,gtoexp,gtocon,scr)
 use Constant, only : maxlq
 implicit real(kind=8) (a-h,o-z)
 dimension         :: mxla(natom), nshlla(0:maxlq,natom,2), ialbs(0:maxlq,natom), ialbc(0:maxlq,natom), iptexp(0:maxlq,natom),  &
                      iptcon(0:maxlq,natom), gtoexp(*), gtocon(*), scr(*)

 rewind(igto)
 mxla   = 0
 nshlla = 0
 ialbs  = 0
 ialbc  = 0
 iptexp = 0
 iptcon = 0
 ngto_tmp = 0
 ncon_tmp = 0
 nbas_tmp = 1
 nbac_tmp = 1
 do ia = 1, natom
   read(igto,*,err=1300,end=1300) i, mxla(ia)
   read(igto,*,err=1300,end=1300) nshlla(0:mxla(ia),ia,1)
   read(igto,*,err=1300,end=1300) nshlla(0:mxla(ia),ia,2)

   do i = 0, mxla(ia)
     if(nshlla(i,ia,1) < 1) cycle

     ialbs(i,ia)  = nbas_tmp
     ialbc(i,ia)  = nbac_tmp
     nbas_tmp = nbas_tmp + (i+i+1)*nshlla(i,ia,1)
     nbac_tmp = nbac_tmp + (i+1)*(i+2)*nshlla(i,ia,1)/2
     iptexp(i,ia) = ngto_tmp + 1
     iptcon(i,ia) = ncon_tmp + 1
     ngto_tmp = ngto_tmp + nshlla(i,ia,1)
     icon = nshlla(i,ia,1)*nshlla(i,ia,2)
     ncon_tmp = ncon_tmp + icon
     read(igto,*,err=1300,end=1300) gtoexp(iptexp(i,ia) : iptexp(i,ia)+nshlla(i,ia,1)-1)
     read(igto,*,err=1300,end=1300) scr(1:icon)
     call Transp(nshlla(i,ia,2),nshlla(i,ia,1),scr,gtocon(iptcon(i,ia)))
   end do
 end do

 return

 1300 write(iout,"(' ### Error in rd_basfun.')")
 call estop(1)
 return
end subroutine rd_basfun

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Search for unique atoms according to their atomic charges and primitive basis functions
!
! In mapatm, positive indices correspond to unique atoms. For example, if the atoms are
!   A1, B1, A2, C1, D1, A3, C2, ...
! mapatm will be
!    1   2  -1   3   4  -1  -3, ...
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine uniatm(natom,nuqatm,mapatm,iza,mxla,nshlla,iptexp,gtoexp)
 use Constant, only : exptol, maxlq
 implicit real(kind=8) (a-h,o-z)
 real(kind=8), intent(in)    :: gtoexp(*)
 integer, intent(out)        :: nuqatm, mapatm(natom)
 integer, intent(in)         :: iza(natom), mxla(natom), nshlla(0:maxlq,natom), iptexp(0:maxlq,natom)

 nuqatm = 0
 mapatm = 0

 do i=1,natom
   if(mapatm(i) /= 0) cycle
   nuqatm = nuqatm + 1
   mapatm(i) = nuqatm

   jsearch : do j=i+1,natom
     if(mapatm(j) /= 0) cycle jsearch
     if(abs(iza(i) - iza(j)) > 0) cycle jsearch
     if(abs(mxla(i) - mxla(j)) > 0) cycle jsearch
     do k=0,mxla(i)
       if(abs(nshlla(k,i) - nshlla(k,j)) > 0) cycle jsearch
     end do
     do k=0,mxla(i)
       ioff=iptexp(k,i)-1
       joff=iptexp(k,j)-1
       do l=1,nshlla(k,i)
         if(abs(gtoexp(ioff+l) - gtoexp(joff+l)) > exptol) cycle jsearch
       end do
     end do

     mapatm(j) = -mapatm(i)
   end do jsearch
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! get atomic index (List(1)) of the Iuq-th unique atom and its equivalents (with -Iuq)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetEqAtm(NAtm,Iuq,IAuq,NeqA,List)
 implicit real(kind=8) (a-h,o-z)
 dimension :: IAuq(*), List(*)

 NeqA = 1
 do i=1,NAtm
   if(IAuq(i) == Iuq) then
     List(1) = i
   else if(-IAuq(i) == Iuq) then
     NeqA = NeqA + 1
     List(NeqA) = i
   end if
 end do

 return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! set nuclear RMS charge radii in au
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine setrms(natom,iza,rms)
 use Constant, only : Two, Three, atom_nucrad
 implicit real(kind=8) (a-h,o-z)
 dimension :: iza(natom), rms(natom)

 q23=sqrt(Two/Three)

 do i=1,natom
   rms(i) = q23 * atom_nucrad(iza(i))
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! check minza
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ChkMinZa(iout,iham,natom,iza,minza,itype,ncalc)
 use Constant, only : pronam
 implicit real(kind=8) (a-h,o-z)
 dimension :: iza(natom)

 ncalc = 0
 if(minza < 0) then
   do iatom = 1, natom
     if(iza(iatom) == abs(minza)) ncalc = ncalc + 1
   end do
   if(ncalc == 0) then
     write(iout,"(/,' No atoms with Za = ',i3)") abs(minza)
     return
   else
     if(iham == 1) then
       write(iout,"(/,' Non-relativistic ',a,' for the atoms with Za = ',i3)") trim(pronam(itype)), abs(minza)
     else
       write(iout,"(/,' Relativistic ',a,' for the atoms with Za = ',i3)") trim(pronam(itype)), abs(minza)
     end if
   end if
 else
   do iatom = 1, natom
     if(iza(iatom) >= minza) ncalc = ncalc + 1
   end do
   if(ncalc == 0) then
     write(iout,"(/,' No atoms with Za > ',i3)") minza-1
     return
   else
     if(iham == 1) then
       write(iout,"(/,' Non-relativistic ',a,' for the atoms with Za > ',i3)") trim(pronam(itype)), minza-1
     else
       write(iout,"(/,' Relativistic ',a,' for the atoms with Za > ',i3)") trim(pronam(itype)), minza-1
     end if
   end if
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! check density matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine checkden(iout,natom,nbas,irelc,nchag,occtot,za,smat,pmat,sden)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension :: za(natom), smat(*), pmat(*), sden(*)

 zatot = sum(za)

 if(irelc == 1) then
   elen=traceSS(nbas,smat,pmat)
   elei=Zero
 else if(irelc == 2) then
   call SglMat(nbas,nbas,pmat,sden)
   elen=traceSS(nbas,smat,sden)
   call SglMat(nbas,nbas,pmat(nbas*nbas*4+1),sden)
   elei=traceSS(nbas,smat,sden)
 end if

 write(iout,"(/,' Checking density matrix')")
 write(iout,"('   Sum of nuclear charges',22x,': ',f9.4)") zatot
 write(iout,"('   Sum of MO occupations',23x,': ',f9.4)") occtot
 if(irelc == 1) then
   write(iout,"('   Analytically integrated number of electrons : ',f9.4)") elen
 else if(irelc == 2) then
   if(elei < 0.0d0) then
     write(iout,"('   Analytically integrated number of electrons : ',f9.4,' - i * ',f6.4)") elen, abs(elei)
   else
     write(iout,"('   Analytically integrated number of electrons : ',f9.4,' + i * ',f6.4)") elen, abs(elei)
   end if
 end if
 nchag = nint(zatot-occtot)
 write(iout,"('   Net charge',34x,': ',f9.4)") dble(nchag)

 if(abs(zatot - occtot) > 6.0d0) write(iout,"(3x,'### Warning in checkden: is ECP used?')")
 if(abs( dble( nint(elen) ) - elen ) > 1.0d-4) then
   write(iout,"(/,' ### Error in checkden: fractional charge is not allowed.')")
   call estop(1)
 end if
 if(abs(elen - occtot) > 1.0d-4) then
   write(iout,"(/,' ### Error in checkden: occtot and elen do not agree.')")
   call estop(1)
 end if
 if(abs(elei) > 1.0d-5) then
   write(iout,"(/,' ### Error in checkden: the imaginary part is non-zero.')")
   call estop(1)
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Delete atomic-overlap blocks in the density matrix. Used for LOCA = 3.
!
! ic can be 1 (scalar) or 2 (two-component).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine abdens(ic,natom,np,mxla,ialbs,nshlla,pden,scr)
 use Constant, only : maxlq, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension     :: mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom), pden(ic,np,ic,np,ic), scr(np,np,*)

 nmat = ic*ic*ic
 isc1 = nmat + 1
 isc2 = isc1 + 1

 ! split pden into scr
 k = 0
 do i3 = 1, ic
   do i2 = 1, ic
     do i1 = 1, ic
       k = k + 1

       do j = 1, np
         do i = 1, np
           scr(i,j,k) = pden(i1,i,i2,j,i3)
         end do
       end do

     end do
   end do
 end do

 do k = 1, nmat
   scr(:,:,isc1) = Zero
   do iatm = 1, natom
     I0 = ialbs(0,iatm)
     Iw = 0
     do lq = 0, mxla(iatm)
       Iw = Iw + (lq+lq+1)*nshlla(lq,iatm)
     end do
     call TakDBlk(np,I0,Iw,scr(1,1,k),scr(1,1,isc2))
     call PutDBlk(np,I0,Iw,scr(1,1,isc1),scr(1,1,isc2))
   end do
   scr(:,:,k) = scr(:,:,isc1)
 end do

 ! save scr to pden
 k = 0
 do i3 = 1, ic
   do i2 = 1, ic
     do i1 = 1, ic
       k = k + 1

       do j = 1, np
         do i = 1, np
           pden(i1,i,i2,j,i3) = scr(i,j,k)
         end do
       end do

     end do
   end do
 end do

 return
end subroutine abdens

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! print Cartesian coordinates of point charges
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine prtpcoord(iout,npchg,pchg,pxyz)
 use Constant, only : au2ang
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: pchg(npchg), pxyz(3,npchg)

  write(iout,"(/,' Cartesian coordinates of point charges (in Angstrom)',/,1x,70('-'),/,3x, &
    'No.       Charge',17x,'X             Y             Z',/,1x,70('-'))")

  do i=1,npchg
    write(iout,"(i6,f13.5,8x,3f14.8)") i,pchg(i),(pxyz(j,i)*au2ang,j=1,3)
  end do

  write(iout,"(1x,70('-'))")

  return
end subroutine prtpcoord

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! print molecular geometry
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine prtcoord(iout,NAtm,IZA,XYZ)
 use Constant, only : au2ang
 implicit real(kind=8) (a-h,o-z)
 integer :: IZA(NAtm)
 real(kind=8) :: XYZ(3,NAtm)
 character*3 :: Elm

  write(iout,"(/,' Cartesian coordinates of atoms (in Angstrom)',/,1x,70('-'),/,3x, &
    'No.   Atom    ZA                 X             Y             Z',/,1x,70('-'))")

  do i=1,NAtm
    call ElemZA(1,Elm,IZA(i))
    write(iout,"(i6,4x,a3,1x,i5,8x,3f14.8)") i,Elm,IZA(i),(XYZ(j,i)*au2ang,j=1,3)
  end do

  write(iout,"(1x,70('-'))")

  return
end subroutine prtcoord

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! calculate mass-unweighted principal axes and moments of inertia
!
! Eigenvalues and Eigenvectors are saved in scr.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RotCons(iout,NAtom,XYZ,cxyz,scr)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: XYZ(3,NAtom),cxyz(3,NAtom),scr(21)

 ! XYZ in CMCS --> cxyz
 call GeomCent(NAtom,XYZ,cxyz,scr)

 ! principal moment of inertia;
 ! Only 3+9+9 elements in of scr will be used.
 IE = 1          ! Eigenvalues
 IR = IE + 3     ! Eigenvectors
 IW = IR + 9
 call MIner(iout,NAtom,cxyz,scr(IR),scr(IE),scr(IW))

 return

 contains

 ! center of geometry
 subroutine GeomCent(NAtm,XYZ,XYZCM,CM)
  use Constant, only : Zero
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: XYZ(3,*),XYZCM(3,*),CM(3)
 
  CM = Zero
  do i=1,NAtm
    CM(1)=CM(1)+XYZ(1,i)
    CM(2)=CM(2)+XYZ(2,i)
    CM(3)=CM(3)+XYZ(3,i)
  end do
  CM = CM / dble(NAtm)
 
  do i=1,NAtm
    XYZCM(1,i)=XYZ(1,i)-CM(1)
    XYZCM(2,i)=XYZ(2,i)-CM(2)
    XYZCM(3,i)=XYZ(3,i)-CM(3)
  end do
 
  return
 end subroutine GeomCent
 
 ! mass-unweighted principal moment of inertia
 subroutine MIner(iout,NAtm,XYZCM,RI,E,WORK)
  use Constant, only : Zero
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: XYZCM(3,NAtm),RI(3,3),E(3),WORK(9)
 
  RI = Zero
  do i=1,NAtm
    X=XYZCM(1,i)
    Y=XYZCM(2,i)
    Z=XYZCM(3,i)
    RI(1,1) = RI(1,1) + Y * Y + Z * Z
    RI(2,2) = RI(2,2) + X * X + Z * Z
    RI(3,3) = RI(3,3) + X * X + Y * Y
    RI(1,2) = RI(1,2) - X * Y
    RI(1,3) = RI(1,3) - X * Z
    RI(2,3) = RI(2,3) - Y * Z
  end do
  RI(2,1) = RI(1,2)
  RI(3,1) = RI(1,3)
  RI(3,2) = RI(2,3)
 
  Call DSYEV('V','L',3,RI,3,E,WORK,9,INFO)
  if(INFO /= 0) then
    write(iout,"(/,' Error in sub. MIner: Diagnolization failed.')")
    call estop(1)
  end if
 
  return
 end subroutine MIner

end subroutine RotCons

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Mode = 0 : returns nuclear charge iza for an element symbol "el".
!     /= 0 : returns element symbol "el" for nuclear charge iza.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ElemZA(Mode,el,iza)
  implicit real(kind=8) (a-h,o-z)
  parameter (maxza=120)
  character*3 :: el,atomlib(maxza)
  data (atomlib(i),i=1,maxza) / &
   'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ',   'NA ','MG ','AL ','SI ','P  ','S  ','CL ','AR ','K  ','CA ', &
   'SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',   'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ', &
   'NB ','MO ','TC ','RU ','RH ','PD ','AG ','CD ','IN ','SN ',   'SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ', &
   'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ',   'LU ','HF ','TA ','W  ','RE ','OS ','IR ','PT ','AU ','HG ', &
   'TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',   'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ', &
   'MD ','NO ','LR ','RF ','DB ','SG ','BH ','HS ','MT ','DS ',   'RG ','CN ','NH ','FL ','MC ','LV ','TS ','OG ','UUE','UBN'/
  save atomlib

  if (Mode == 0) then

    call charl2u(el)
    iza = 0
    do i=1,maxza
      if(index(el,atomlib(i)) /= 0)then
        iza = i
        exit
      end if
    end do

  else

    el = "???"
    if(iza > 0 .and. iza <= maxza) el = adjustl(atomlib(iza))
    call charu2l(el(2:3))

  end if

  return
end subroutine ElemZA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! nuclear quadrupole moments (Q in millibarn). From
! P. Pyykko. Year-2017 nuclear quadrupole moments. Mol. Phys. 116, 1328 (2018).
!
! Update.
!
! Bi-209.
! J.-P. Dognon and P. Pyykko, Determining nuclear quadrupole moments of Bi and Sb from molecular data, Phys. Chem. Chem. Phys.
! 25, 2758-2761 (2023).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NQMlib(iz,ndat,iso,eqv)
 implicit real(kind=8) (a-h,o-z)
 integer      :: iso(*)
 real(kind=8) :: eqv(*)

 ndat = 0
 select case(iz)
   case(01)
     ndat = 1
     iso(1:ndat) = (/2/)
     eqv(1:ndat) = (/2.85783d0/)
   case(03)
     ndat = 2
     iso(1:ndat) = (/6,7/)
     eqv(1:ndat) = (/-0.808d0,-40.1d0/)
   case(04)
     ndat = 1
     iso(1:ndat) = (/9/)
     eqv(1:ndat) = (/52.88d0/)
   case(05)
     ndat = 2
     iso(1:ndat) = (/10,11/)
     eqv(1:ndat) = (/84.59d0,40.59d0/)
   case(06)
     ndat = 1
     iso(1:ndat) = (/11/)
     eqv(1:ndat) = (/33.27d0/)
   case(07)
     ndat = 1
     iso(1:ndat) = (/14/)
     eqv(1:ndat) = (/20.44d0/)
   case(08)
     ndat = 1
     iso(1:ndat) = (/17/)
     eqv(1:ndat) = (/-25.58d0/)
   case(09)
     ndat = 1
     iso(1:ndat) = (/19/)
     eqv(1:ndat) = (/-94.2d0/)
   case(10)
     ndat = 1
     iso(1:ndat) = (/21/)
     eqv(1:ndat) = (/101.55d0/)
   case(11)
     ndat = 1
     iso(1:ndat) = (/23/)
     eqv(1:ndat) = (/104.0d0/)
   case(12)
     ndat = 1
     iso(1:ndat) = (/25/)
     eqv(1:ndat) = (/199.4d0/)
   case(13)
     ndat = 1
     iso(1:ndat) = (/27/)
     eqv(1:ndat) = (/146.6d0/)
   case(16)
     ndat = 2
     iso(1:ndat) = (/33,35/)
     eqv(1:ndat) = (/-69.4d0,48.3d0/)
   case(17)
     ndat = 2
     iso(1:ndat) = (/35,37/)
     eqv(1:ndat) = (/-81.12d0,-63.93d0/)
   case(19)
     ndat = 3
     iso(1:ndat) = (/39,40,41/)
     eqv(1:ndat) = (/60.3d0,-75.0d0,73.4d0/)
   case(20)
     ndat = 2
     iso(1:ndat) = (/41,43/)
     eqv(1:ndat) = (/-66.5d0,-40.8d0/)
   case(21)
     ndat = 1
     iso(1:ndat) = (/45/)
     eqv(1:ndat) = (/-220.0d0/)
   case(22)
     ndat = 2
     iso(1:ndat) = (/47,49/)
     eqv(1:ndat) = (/302.0d0,247.0d0/)
   case(23)
     ndat = 2
     iso(1:ndat) = (/50,51/)
     eqv(1:ndat) = (/210.0d0,-52.0d0/)
   case(24)
     ndat = 1
     iso(1:ndat) = (/53/)
     eqv(1:ndat) = (/-150.0d0/)
   case(25)
     ndat = 1
     iso(1:ndat) = (/55/)
     eqv(1:ndat) = (/330.0d0/)
   case(26)
     ndat = 1
     iso(1:ndat) = (/57/)
     eqv(1:ndat) = (/160.0d0/)
   case(27)
     ndat = 1
     iso(1:ndat) = (/59/)
     eqv(1:ndat) = (/420.0d0/)
   case(28)
     ndat = 1
     iso(1:ndat) = (/61/)
     eqv(1:ndat) = (/162.0d0/)
   case(29)
     ndat = 2
     iso(1:ndat) = (/63,65/)
     eqv(1:ndat) = (/-220.0d0,-204.0d0/)
   case(30)
     ndat = 1
     iso(1:ndat) = (/67/)
     eqv(1:ndat) = (/122.0d0/)
   case(31)
     ndat = 2
     iso(1:ndat) = (/69,71/)
     eqv(1:ndat) = (/171.0d0,107.0d0/)
   case(32)
     ndat = 1
     iso(1:ndat) = (/73/)
     eqv(1:ndat) = (/-196.0d0/)
   case(33)
     ndat = 1
     iso(1:ndat) = (/75/)
     eqv(1:ndat) = (/311.0d0/)
   case(34)
     ndat = 1
     iso(1:ndat) = (/77/)
     eqv(1:ndat) = (/760.0d0/)
   case(35)
     ndat = 2
     iso(1:ndat) = (/79,81/)
     eqv(1:ndat) = (/308.7d0,257.9d0/)
   case(36)
     ndat = 1
     iso(1:ndat) = (/83/)
     eqv(1:ndat) = (/259.0d0/)
   case(37)
     ndat = 2
     iso(1:ndat) = (/85,87/)
     eqv(1:ndat) = (/276.0d0,133.5d0/)
   case(38)
     ndat = 1
     iso(1:ndat) = (/87/)
     eqv(1:ndat) = (/305.0d0/)
   case(39)
     ndat = 1
     iso(1:ndat) = (/90/)
     eqv(1:ndat) = (/-125.0d0/)
   case(40)
     ndat = 1
     iso(1:ndat) = (/91/)
     eqv(1:ndat) = (/-176.0d0/)
   case(41)
     ndat = 1
     iso(1:ndat) = (/93/)
     eqv(1:ndat) = (/-320.0d0/)
   case(42)
     ndat = 2
     iso(1:ndat) = (/95,97/)
     eqv(1:ndat) = (/-22.0d0,255.0d0/)
   case(43)
     ndat = 1
     iso(1:ndat) = (/99/)
     eqv(1:ndat) = (/-129.0d0/)
   case(44)
     ndat = 2
     iso(1:ndat) = (/99,101/)
     eqv(1:ndat) = (/79.0d0,457.0d0/)
   case(45)
     ndat = 1
     iso(1:ndat) = (/100/)
     eqv(1:ndat) = (/153.0d0/)
   case(46)
     ndat = 1
     iso(1:ndat) = (/105/)
     eqv(1:ndat) = (/660.0d0/)
   case(49)
     ndat = 2
     iso(1:ndat) = (/113,115/)
     eqv(1:ndat) = (/761.0d0,772.0d0/)
   case(50)
     ndat = 1
     iso(1:ndat) = (/119/)
     eqv(1:ndat) = (/-132.0d0/)
   case(51)
     ndat = 2
     iso(1:ndat) = (/121,123/)
     eqv(1:ndat) = (/-543.0d0,-692.0d0/)
   case(53)
     ndat = 1
     iso(1:ndat) = (/127/)
     eqv(1:ndat) = (/-688.22d0/)
   case(54)
     ndat = 1
     iso(1:ndat) = (/131/)
     eqv(1:ndat) = (/-114.6d0/)
   case(55)
     ndat = 1
     iso(1:ndat) = (/133/)
     eqv(1:ndat) = (/-3.43d0/)
   case(56)
     ndat = 2
     iso(1:ndat) = (/135,137/)
     eqv(1:ndat) = (/153.0d0,236.0d0/)
   case(57)
     ndat = 2
     iso(1:ndat) = (/138,139/)
     eqv(1:ndat) = (/450.0d0,206.0d0/)
   case(59)
     ndat = 1
     iso(1:ndat) = (/141/)
     eqv(1:ndat) = (/-58.9d0/)
   case(60)
     ndat = 2
     iso(1:ndat) = (/143,145/)
     eqv(1:ndat) = (/-630.0d0,-330.0d0/)
   case(61)
     ndat = 1
     iso(1:ndat) = (/147/)
     eqv(1:ndat) = (/740.0d0/)
   case(62)
     ndat = 2
     iso(1:ndat) = (/147,149/)
     eqv(1:ndat) = (/-259.0d0,75.0d0/)
   case(63)
     ndat = 2
     iso(1:ndat) = (/151,153/)
     eqv(1:ndat) = (/903.0d0,2412.0d0/)
   case(64)
     ndat = 2
     iso(1:ndat) = (/155,157/)
     eqv(1:ndat) = (/1270.0d0,1350.0d0/)
   case(65)
     ndat = 1
     iso(1:ndat) = (/159/)
     eqv(1:ndat) = (/1432.0d0/)
   case(66)
     ndat = 2
     iso(1:ndat) = (/161,163/)
     eqv(1:ndat) = (/2507.0d0,2648.0d0/)
   case(67)
     ndat = 1
     iso(1:ndat) = (/165/)
     eqv(1:ndat) = (/3580.0d0/)
   case(68)
     ndat = 1
     iso(1:ndat) = (/167/)
     eqv(1:ndat) = (/3565.0d0/)
   case(70)
     ndat = 1
     iso(1:ndat) = (/173/)
     eqv(1:ndat) = (/2800.0d0/)
   case(71)
     ndat = 2
     iso(1:ndat) = (/175,176/)
     eqv(1:ndat) = (/3490.0d0,4970.0d0/)
   case(72)
     ndat = 2
     iso(1:ndat) = (/177,179/)
     eqv(1:ndat) = (/3365.0d0,3793.0d0/)
   case(73)
     ndat = 1
     iso(1:ndat) = (/181/)
     eqv(1:ndat) = (/3170.0d0/)
   case(75)
     ndat = 2
     iso(1:ndat) = (/185,187/)
     eqv(1:ndat) = (/2180.0d0,2070.0d0/)
   case(76)
     ndat = 1
     iso(1:ndat) = (/189/)
     eqv(1:ndat) = (/856.0d0/)
   case(77)
     ndat = 2
     iso(1:ndat) = (/191,193/)
     eqv(1:ndat) = (/816.0d0,751.0d0/)
   case(79)
     ndat = 1
     iso(1:ndat) = (/197/)
     eqv(1:ndat) = (/547.0d0/)
   case(80)
     ndat = 1
     iso(1:ndat) = (/201/)
     eqv(1:ndat) = (/387.0d0/)
   case(82)
     ndat = 1
     iso(1:ndat) = (/209/)
     eqv(1:ndat) = (/-269.0d0/)
   case(83)
     ndat = 1
     iso(1:ndat) = (/209/)
     ! eqv(1:ndat) = (/-516.0d0/)
     ! new NQM value by Pyykko. see PCCP 25, 2758 (2023).
     eqv(1:ndat) = (/-422.0d0/)
   case(86)
     ndat = 1
     iso(1:ndat) = (/209/)
     eqv(1:ndat) = (/311.0d0/)
   case(87)
     ndat = 1
     iso(1:ndat) = (/211/)
     eqv(1:ndat) = (/-210.0d0/)
   case(89)
     ndat = 1
     iso(1:ndat) = (/227/)
     eqv(1:ndat) = (/1700.0d0/)
   case(90)
     ndat = 1
     iso(1:ndat) = (/229/)
     eqv(1:ndat) = (/3110.0d0/)
   case(91)
     ndat = 1
     iso(1:ndat) = (/231/)
     eqv(1:ndat) = (/-1720.0d0/)
   case(92)
     ndat = 2
     iso(1:ndat) = (/233,235/)
     eqv(1:ndat) = (/3663.0d0,4936.0d0/)
   case(93)
     ndat = 1
     iso(1:ndat) = (/237/)
     eqv(1:ndat) = (/3886.0d0/)
   case(94)
     ndat = 1
     iso(1:ndat) = (/241/)
     eqv(1:ndat) = (/5600.0d0/)
   case(95)
     ndat = 2
     iso(1:ndat) = (/241,243/)
     eqv(1:ndat) = (/4200.0d0,4200.0d0/)
   case(99)
     ndat = 1
     iso(1:ndat) = (/253/)
     eqv(1:ndat) = (/6700.0d0/)
 end select

 return
end subroutine NQMlib

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! print L number of real matrices in A(M,N,L) in the format of nD14.6.
!
! mxcol: the max number of column to be printed
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine writemat(m,n,l,a,mxcol)
 implicit real(kind=8) (a-h,o-z)
 dimension     :: a(m,n,l)
 character*10  :: frmt

 !frmt='(   d14.6)'
 frmt='(   f14.4)'
 write(frmt(2:4),"(i3)") max(1,min(mxcol,n))

 do im = 1, l
   write(*,"('im = ',i2)") im
   do i=1,m
     write(*,trim(frmt)) a(i,:,im)
   end do
 end do

 return
end subroutine writemat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! print L number of complex matrices in A(M,N,L) in the format of nD14.6.
!
! mxcol: the max number of column to be printed
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine writecmat(m,n,l,a,mxcol)
 implicit real(kind=8) (a-h,o-z)
 dimension     :: a(2,m,n,l)
 character*10  :: frmt

 frmt='(   d14.6)'
 write(frmt(2:4),"(i3)") max(1,min(mxcol,n))

 write(*,"('real')")
 do im = 1, l
   write(*,"('im = ',i2)") im
   do i=1,m
     write(*,trim(frmt)) a(1,i,:,im)
   end do
 end do

 write(*,"('imag')")
 do im = 1, l
   write(*,"('im = ',i2)") im
   do i=1,m
     write(*,trim(frmt)) a(2,i,:,im)
   end do
 end do

 return
end subroutine writecmat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read an <ENTER> and (if imod /= 0) stop.
! stop can trigger a message "Warning: ieee_inexact is signaling" by pgf90.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine estop(imod)
 implicit real(kind=8) (a-h,o-z)

 !write(*,"(//,' Press <ENTER> to exit',/)")
 !read(*,*)

 if(imod /= 0) stop

 return
end

!--- END
