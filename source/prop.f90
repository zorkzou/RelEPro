!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Effective density calculation (local version): in tr(P * dNAI), only calculate atomic blocks of dNAI
!
! Input:
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   irelc          : scalar or 2-component calculation
!   natom          : number of atoms
!   iham           : NR or X2C Hamiltonians
!   minza          : min. atomic nuclear charge for ED calculation
!   loca           : 1CA approximation (see the manual)
!   iza,za         : atomic nuclear charges
!   rms            : nuclear RMS charge radii
!   npsph          : number of primitive spherical functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   ialbs          : starting position of each L in primitive spherical basis functions
!   iptexp         : starting positions of GTO exponents for each L
!   gtoexp         : GTO exponents
!   c2scoef        : Car. to sph. coefficients
!   pv,pw          : transformed density matrices for V and W
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine edcalclo(iout,mpatt,irelc,natom,iham,minza,loca,iza,za,rms,npsph,mxla,nshlla,ialbs,iptexp,gtoexp,c2scoef,pv,pw,  &
 scr,lscr)
 use Constant, only : maxlq, dpi, au2ang, f_4cc, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension     :: iza(natom), za(natom), rms(natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom),  &
                  iptexp(0:maxlq,natom), gtoexp(*), c2scoef(*), pv(npsph*npsph,*), pw(npsph*npsph,*), scr(lscr)
 character*3   :: Elm
 allocatable   :: xyz0(:)

 if(iham == 1) then
   if(loca > 1) write(iout,"(/,' The Atomic block approximation is turned on in P.')")
 else
   if(loca == 2) then
     write(iout,"(/,' The Atomic block approximation is turned on in Pv and Pw.')")
   else if(loca == 3) then
     write(iout,"(/,' The Atomic block approximation is turned on in P.')")
   end if
 end if

 fac1 = sqrt(1.5d0) * 1.0d5 * au2ang
 facg = fac1 * sqrt(2.0d0/3.0d0)
 allocate(xyz0(3))
 xyz0 = Zero

 ! check minza
 call ChkMinZa(iout,iham,natom,iza,minza,4,ncalc)
 if(ncalc == 0) return

 ! single-value density matrices
 if(irelc == 1) then
   lden = 0
 else if(irelc == 2) then
   npss = npsph*npsph
   lden = npss
   ipv  = 1
   ! Only the real part of Pv should be taken
   call SglMat(npsph,npsph,pv,scr(ipv))
   ! Split Pw into Pwo, Pwx, Pwy, and Pwz
   if(iham > 1) then
     lden = npss*5
     ipwo = ipv  + npss
     ipwx = ipwo + npss
     ipwy = ipwx + npss
     ipwz = ipwy + npss
     call SGLCOMPW(npsph,pw(1,1),pw(1,5),scr(ipwo),scr(ipwx),scr(ipwy),scr(ipwz))
   end if
 end if

 write(iout,"(1x,63('-'),/, '   No.   Atom   ZA     RMS (fm)      RG (fm)          Rho (a.u.)',/, 1x,63('-'))")
 isc1 = lden +   1
 do iatom = 1, natom
   if(minza < 0) then
     if(iza(iatom) /= abs(minza)) cycle
   else
     if(iza(iatom) < minza) cycle
   end if

   nshell0 = 0
   I0 = ialbs(0,iatom)
   Ic = 0
   Ip = 0
   do i = 0, mxla(iatom)
     ngto = nshlla(i,iatom)
     nshell0  = nshell0  + ngto
     Ic = Ic + (i+1)*(i+2)*ngto/2
     Ip = Ip + (i+i+1)*ngto
   end do
   ntt = Ic*(Ic+1)/2
   nss = Ic*Ic
   lnai = ntt*irelc*irelc
   iscr = isc1 + lnai
   isc2 = isc1 + nss
   isc3 = isc2 + nss
   isc4 = isc3 + nss
   if(iham > 1 .and. irelc == 2) then
     isc5 = isc4 + nss
     isc6 = isc5 + nss
     isc7 = isc6 + nss
     isc8 = isc7 + nss
     isc9 = isc8 + nss
     isc10= isc9 + nss
   end if
   call pnai_drv(iout,0,0,0,1,mpatt,1,za(iatom),xyz0,mxla(iatom),nshell0,Ic,mxla(iatom),nshlla(0,iatom),  &
     gtoexp(iptexp(0,iatom)),rms(iatom),scr(isc1),scr(iscr),lscr-lnai-lden)
   call lt2sqr(0,Ic,scr(isc1),scr(isc2))
   ! dV/dR0 --> sc1
   call car2sph(0,1,mxla(iatom),nshlla(0,iatom),Ic,Ip,scr(isc2),scr(isc1),c2scoef,scr(isc3),scr(isc4))
   if(irelc == 1) then
     call TakDBlk(npsph,I0,Ip,pv,scr(isc3))
   else if(irelc == 2) then
     call TakDBlk(npsph,I0,Ip,scr(ipv),scr(isc3))
   end if
   val_cd = traceSS(Ip,scr(isc3),scr(isc1))

   if(iham > 1) then
     call pnai_drv(iout,irelc,0,0,1,mpatt,1,za(iatom),xyz0,mxla(iatom),nshell0,Ic,mxla(iatom),nshlla(0,iatom),  &
       gtoexp(iptexp(0,iatom)),rms(iatom),scr(isc1),scr(iscr),lscr-lnai-lden)
     if(irelc == 1) then
       ! dW/dR0 --> sc1 (1c)
       call lt2sqr(0,Ic,scr(isc1),scr(isc2))
       call car2sph(0,1,mxla(iatom),nshlla(0,iatom),Ic,Ip,scr(isc2),scr(isc1),c2scoef,scr(isc3),scr(isc4))
       call TakDBlk(npsph,I0,Ip,pw,scr(isc3))
       val_cd = val_cd + traceSS(Ip,scr(isc3),scr(isc1)) / f_4cc
     else if(irelc == 2) then
       ! dW/dR0 --> sc1-sc4 (2c)
       do ipvp = 1, 4
         j = (ipvp+1)/3
         jsc1 = isc1 + ntt * (ipvp - 1)
         jsc5 = isc5 + nss * (ipvp - 1)
         call lt2sqr(j,Ic,scr(jsc1),scr(jsc5))
       end do
       do ipvp = 1, 4
         jsc1 = isc1 + nss * (ipvp - 1)
         jsc5 = isc5 + nss * (ipvp - 1)
         call car2sph(0,1,mxla(iatom),nshlla(0,iatom),Ic,Ip,scr(jsc5),scr(jsc1),c2scoef,scr(isc9),scr(isc10))
       end do
       call TakDBlk(npsph,I0,Ip,scr(ipwo),scr(isc5))
       call TakDBlk(npsph,I0,Ip,scr(ipwx),scr(isc6))
       call TakDBlk(npsph,I0,Ip,scr(ipwy),scr(isc7))
       call TakDBlk(npsph,I0,Ip,scr(ipwz),scr(isc8))
       val_cd = val_cd + traceSS(Ip,scr(isc5),scr(isc1)) / f_4cc
       ! * -1 due to antisymmetry
       val_cd = val_cd - traceSS(Ip,scr(isc6),scr(isc2)) / f_4cc
       val_cd = val_cd - traceSS(Ip,scr(isc7),scr(isc3)) / f_4cc
       val_cd = val_cd - traceSS(Ip,scr(isc8),scr(isc4)) / f_4cc
     end if
   end if

   call ElemZA(1,Elm,iza(iatom),Elm)
   val_cd = val_cd / (rms(iatom) * za(iatom) * dpi)
   write(iout,"(i6,4x,a3,i5,2f13.5,f20.5)") iatom, Elm, iza(iatom), rms(iatom)*fac1, rms(iatom)*facg, val_cd
 end do
 write(iout,"(1x,63('-'))")
 deallocate(xyz0)

 return
end subroutine edcalclo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Effective density calculation
!
! Input:
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   irelc          : scalar or 2-component calculation
!   natom          : number of atoms
!   iham           : NR or X2C Hamiltonians
!   minza          : min. atomic nuclear charge for ED calculation
!   ipop           : population if =1
!   iza,za         : atomic nuclear charges
!   rms            : nuclear RMS charge radii
!   xyz            : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   npcar          : number of primitive Cartesian functions
!   npsph          : number of primitive spherical functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   pv,pw          : transformed density matrices for V and W
!   nmo,occ,ene,ispn,symm
!                  : information of MOs
!   cr,cu          : transformed MOs
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine edcalc(iout,mpatt,irelc,natom,iham,minza,ipop,iza,za,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,  &
 c2scoef,pv,pw, nmo,occ,ene,ispn,symm,cr,cu, scr,lscr)
 use Constant, only : maxlq, dpi, au2ang, f_4cc, One, Zero, occtol
 implicit real(kind=8) (a-h,o-z)
 character*10  :: symm
 dimension     :: iza(natom), za(natom), rms(natom), xyz(3,natom), mxla(natom), nshlla(0:maxlq,natom), gtoexp(*), c2scoef(*),  &
                  pv(npsph*npsph,*), pw(npsph*npsph,*),  &
                  occ(nmo), ene(nmo), ispn(nmo), symm(nmo), cr(npsph*irelc*nmo,irelc), cu(npsph*irelc*nmo,irelc), scr(lscr)
 character*3   :: Elm
 allocatable   :: I0(:),IW(:),popden(:,:)

 fac1 = sqrt(1.5d0) * 1.0d5 * au2ang
 facg = fac1 * sqrt(2.0d0/3.0d0)
 ntt  = npcar*(npcar+1)/2
 nss  = npcar*npcar
 npss = npsph*npsph
 lnai = ntt*irelc*irelc

 ! single-value density matrices
 if(irelc == 1) then
   lden = 0
 else if(irelc == 2) then
   ipv  = 1
   lden = npss
   ! Only the real part of Pv should be taken
   call SglMat(npsph,npsph,pv,scr(ipv))
   ! Split Pw into Pwo, Pwx, Pwy, and Pwz
   if(iham > 1) then
     ipwo = ipv  + npss
     ipwx = ipwo + npss
     ipwy = ipwx + npss
     ipwz = ipwy + npss
     lden = lden + npss*4
     call SGLCOMPW(npsph,pw(1,1),pw(1,5),scr(ipwo),scr(ipwx),scr(ipwy),scr(ipwz))
   end if
 end if

 isc1 = lden +   1
 iscr = isc1 + lnai
 isc2 = isc1 + nss
 isc3 = isc2 + nss
 isc4 = isc3 + nss
 if(irelc == 2) then
   isc5 = isc4 + nss
   isc6 = isc5 + nss
   isc7 = isc6 + nss
   isc8 = isc7 + nss
   isc9 = isc8 + nss
   isc10= isc9 + nss
   isc11= isc10+ nss
   isc12= isc11+ nss
 end if

 ! check minza
 call ChkMinZa(iout,iham,natom,iza,minza,4,ncalc)
 if(ncalc == 0) return

 ! starting position and width of primitive sperical functions for each atom
 if(ipop == 1) then
   allocate(I0(natom), IW(natom), popden(natom,nmo))

   do iatom = 1, natom
     if(iatom == 1) then
       I0(iatom) = 1
     else
       I0(iatom) = I0(iatom-1) + Iw(iatom-1)
     end if
     Iw(iatom) = 0
     do i = 0, mxla(iatom)
       Iw(iatom) = Iw(iatom) + (i+i+1)*nshlla(i,iatom)
     end do
   end do
   matom = min(4,natom)

   write(iout,"(1x,132('-'),/,'   No.   Atom   ZA     RMS (fm)      RG (fm)          Rho (a.u.)',50x,  &
     'Population on atoms',/,1x,132('-'))")
 else
   write(iout,"(1x,63('-'),/, '   No.   Atom   ZA     RMS (fm)      RG (fm)          Rho (a.u.)',/, 1x,63('-'))")
 end if

 do iatom = 1, natom
   if(minza < 0) then
     if(iza(iatom) /= abs(minza)) cycle
   else
     if(iza(iatom) < minza) cycle
   end if

   call pnai_drv(iout,0,0,0,iatom,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,scr(isc1),scr(iscr),  &
     lscr-lnai-lden)
   call lt2sqr(0,npcar,scr(isc1),scr(isc2))
   ! dV/dR0 --> sc1
   call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc2),scr(isc1),c2scoef,scr(isc3),scr(isc4))
   if(irelc == 1) then
     val_cd = traceSS(npsph,pv,scr(isc1))
   else if(irelc == 2) then
     val_cd = traceSS(npsph,scr(ipv),scr(isc1))
   end if
   ! population
   if(ipop == 1) then
     popden = Zero
     if(irelc == 1) then
       call cdpop(npsph,nmo,natom,One,I0,Iw,scr(isc1),occ,cr,popden,scr(isc2))
     else if(irelc == 2) then
       call DblMat(npsph,npsph,scr(isc1),scr(isc2))
       call cdpop2(npsph*2,nmo,natom,One,I0,Iw,scr(isc2),occ,cr(1,1),cr(1,1),popden,scr(isc6))
       call cdpop2(npsph*2,nmo,natom,One,I0,Iw,scr(isc2),occ,cr(1,2),cr(1,2),popden,scr(isc6))
     end if
   end if

   if(iham > 1) then
     call pnai_drv(iout,irelc,0,0,iatom,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,scr(isc1),scr(iscr),  &
       lscr-lnai-lden)
     if(irelc == 1) then
       ! dW/dR0 --> sc1
       call lt2sqr(0,npcar,scr(isc1),scr(isc2))
       call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc2),scr(isc1),c2scoef,scr(isc3),scr(isc4))
       val_cd = val_cd + traceSS(npsph,pw,scr(isc1)) / f_4cc
     else if(irelc == 2) then
       ! dW/dR0 --> sc1-sc4
       do ipvp = 1, 4
         j = (ipvp+1)/3
         jsc1 = isc1 + ntt * (ipvp - 1)
         jsc5 = isc5 + nss * (ipvp - 1)
         call lt2sqr(j,npcar,scr(jsc1),scr(jsc5))
       end do
       do ipvp = 1, 4
         jsc1 = isc1 + nss * (ipvp - 1)
         jsc5 = isc5 + nss * (ipvp - 1)
         call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(jsc5),scr(jsc1),c2scoef,scr(isc9),scr(isc10))
       end do
       val_cd = val_cd + traceSS(npsph,scr(ipwo),scr(isc1)) / f_4cc
       ! * -1 due to antisymmetry
       val_cd = val_cd - traceSS(npsph,scr(ipwx),scr(isc2)) / f_4cc
       val_cd = val_cd - traceSS(npsph,scr(ipwy),scr(isc3)) / f_4cc
       val_cd = val_cd - traceSS(npsph,scr(ipwz),scr(isc4)) / f_4cc
     end if
     ! population
     if(ipop == 1) then
       if(irelc == 1) then
         call cdpop(npsph,nmo,natom,(One/f_4cc),I0,Iw,scr(isc1),occ,cu,popden,scr(isc2))
       else if(irelc == 2) then
         ! Re(C)' * Re(W) * Re(C) + Im(C)' * Re(W) * Im(C)
         call WoaddWy(npsph,scr(isc1),scr(isc3),scr(isc5))
         call cdpop2(npsph*2,nmo,natom,( One/f_4cc),I0,Iw,scr(isc5),occ,cu(1,1),cu(1,1),popden,scr(isc9))
         call cdpop2(npsph*2,nmo,natom,( One/f_4cc),I0,Iw,scr(isc5),occ,cu(1,2),cu(1,2),popden,scr(isc9))
         ! Im(C)' * Im(W) * Re(C) - Re(C)' * Im(W) * Im(C)
         call WxaddWz(npsph,scr(isc2),scr(isc4),scr(isc5))
         call cdpop2(npsph*2,nmo,natom,( One/f_4cc),I0,Iw,scr(isc5),occ,cu(1,2),cu(1,1),popden,scr(isc9))
         call cdpop2(npsph*2,nmo,natom,(-One/f_4cc),I0,Iw,scr(isc5),occ,cu(1,1),cu(1,2),popden,scr(isc9))
       end if
     end if
   end if

   call ElemZA(1,Elm,iza(iatom),Elm)
   fac2 = rms(iatom) * za(iatom) * dpi
   val_cd = val_cd / fac2
   write(iout,"(i6,4x,a3,i5,2f13.5,f20.5)") iatom, Elm, iza(iatom), rms(iatom)*fac1, rms(iatom)*facg, val_cd
   if(ipop == 1) then
     popden = popden / fac2
     if(irelc == 1) then
       write(iout,"(/,11x,'#MO      Irrep  S          Ene        Occ          Rho sum         Rho pop')")
       do imo = 1, nmo
         if(abs(occ(imo)) < occtol) cycle
         write(iout,"(8x,i6,1x,a10,i3,f13.4,f11.6,1x,5f16.5)") imo, adjustr(symm(imo)), ispn(imo), ene(imo), occ(imo),  &
           sum(popden(:,imo)), popden(1:matom,imo)
         if(natom > 4) write(iout,"(69x,4f16.5)") popden(matom+1:natom,imo)
       end do
     else if(irelc == 2) then
       write(iout,"(/,11x,'#MO      Irrep             Ene        Occ          Rho sum         Rho pop')")
       do imo = 1, nmo
         if(abs(occ(imo)) < occtol) cycle
         write(iout,"(8x,i6,1x,a10,3x,f13.4,f11.6,1x,5f16.5)") imo, adjustr(symm(imo)), ene(imo), occ(imo),  &
           sum(popden(:,imo)), popden(1:matom,imo)
         if(natom > 4) write(iout,"(69x,4f16.5)") popden(matom+1:natom,imo)
       end do
     end if
     write(iout,*)
   end if
 end do  ! iatom

 if(ipop == 1) then
   write(iout,"(1x,132('-'))")
   deallocate(I0, IW, popden)
 else
   write(iout,"(1x,63('-'))")
 end if

 return
end subroutine edcalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! effective density population for each atom
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine cdpop(npsph,nmo,natom,fac,Ip,Iw,dnai,occ,cmoL,popden,cmoR)
 use Constant, only : occtol, One, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension         :: Ip(natom), Iw(natom), dnai(npsph,npsph), occ(nmo), cmoL(npsph,nmo), popden(natom,nmo), cmoR(npsph,nmo)

 call MatMult(1,npsph,npsph,nmo,One,Zero,dnai,cmoL,cmoR)

 do imo = 1, nmo
   if(abs(occ(imo)) < occtol) cycle

   do iatom = 1, natom
     popden(iatom,imo) = popden(iatom,imo) + fac * occ(imo) * dotx(Iw(iatom), cmoL(Ip(iatom),imo), cmoR(Ip(iatom),imo))
   end do
 end do

 return
end subroutine cdpop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! effective density population for each atom. 2-component version.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine cdpop2(npd,nmo,natom,fac,Ip,Iw,dnai,occ,cmoL,cmoR,popden,scr)
 use Constant, only : occtol, One, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension         :: Ip(natom), Iw(natom), dnai(npd,npd), occ(nmo), cmoL(npd,nmo), cmoR(npd,nmo), popden(natom,nmo), scr(npd,nmo)

 call MatMult(1,npd,npd,nmo,One,Zero,dnai,cmoR,scr)

 do imo = 1, nmo
   if(abs(occ(imo)) < occtol) cycle

   do iatom = 1, natom
     popden(iatom,imo) = popden(iatom,imo) + fac * occ(imo) * dotx(Iw(iatom)*2, cmoL(Ip(iatom)*2-1,imo), scr(Ip(iatom)*2-1,imo))
   end do
 end do

 return
end subroutine cdpop2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! contact density calculation.
!
! Input:
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   irelc          : scalar or 2-component calculation (NYI)
!   natom          : number of atoms
!   iham           : NR or X2C Hamiltonians
!   minza          : min. atomic nuclear charge for CD calculation
!   iza            : atomic nuclear charges
!   rms            : nuclear RMS charge radii
!   xyz            : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   npcar          : number of primitive Cartesian functions
!   npsph          : number of primitive spherical functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   pv             : P (n.r.) or P_LV (r.) density matrix in primitive spherical functions
!   pw             : P_LW (r.) density matrix in primitive spherical functions
!   xmat           : X matrix in primitive spherical functions (not used)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine cdcalc(iout,mpatt,irelc,natom,iham,minza,iza,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,  &
 pv,pw,xmat, scr,lscr)
 use Constant, only : maxlq, au2ang, Zero, One, clight, f_4cc, fpi
 implicit real(kind=8) (a-h,o-z)
 !<<< curve data
 parameter(delx=0.0d0,dely=0.0d0,delz=2.0d-2)
 ! parameter(npt=250)
 parameter(npt=1)
 !>>>
 character*10  :: symm
 dimension     :: iza(natom), rms(natom), xyz(3,natom), mxla(natom), nshlla(0:maxlq,natom), gtoexp(*), c2scoef(*),  &
                  pv(npsph*npsph,*), pw(npsph*npsph,*), xmat(npsph*npsph,*), scr(lscr)
 allocatable   :: aol(:,:), aos(:,:), px(:,:), LPat(:,:), fnorm(:), pnt6(:,:), rpa(:)
 character*3   :: Elm

 fac1 = sqrt(1.5d0) * 1.0d5 * au2ang

 isc0 = 1
 isc1 = isc0 + npt
 isc2 = isc1 + npsph*npt
 isc3 = isc2 + npsph*npt
 iend = isc3 + npsph*npt
 if(iend > lscr) then
   write(iout,"(' ### Error: the scratch space in sub. cdcalc is too small.')")
   call estop(1)
 end if

 ! check minza
 call ChkMinZa(iout,iham,natom,iza,minza,5,ncalc)
 if(ncalc == 0) return

 MaxSum=(MaxTyp+1)*(MaxTyp+2)*(MaxTyp+3)/6
 allocate(aol(npcar,npt), LPat(3,MaxSum), fnorm(MaxSum), pnt6(3,npt), rpa(3))
 if(irelc == 1) then
   if(iham > 1) allocate( aos(npcar,npt) )
 else if(irelc == 2) then
   ! N.A.
   write(iout,"(/,' Two-component CD is NYI.')")
   return
 end if

 ! Generate (x,y,z) patterns in MOLDEN or Gaussian mode
 call GnPatt(MaxTyp,mpatt,LPat)

 ! common normalization factors
 call NormFac(MaxTyp,LPat,fnorm)

 write(iout,"(1x,63('-'),/, '   No.   Atom   ZA     RMS (fm)                       Rho (a.u.)',/, 1x,63('-'))")
 do iatom = 1, natom
   if(minza < 0) then
     if(iza(iatom) /= abs(minza)) cycle
   else
     if(iza(iatom) < minza) cycle
   end if

   ! point-1 at the atomic center
   pnt6(:,1) = xyz(:,iatom)
   rpa = (/delx,dely,delz/)
   do ipt = 2, npt
     pnt6(:,ipt) = pnt6(:,ipt-1) + rpa
   end do

   ! large component
   iao = 0
   jshl = 0
   do iam = 1, natom
     do lq = 0, mxla(iam)
       noff = (lq+0)*(lq+1)*(lq+2)/6
       npatt= (lq+1)*(lq+2)/2
       do ishl = 1, nshlla(lq,iam)
         jshl = jshl + 1
         x1 = gtoexp(jshl)**(0.75d0)
         If(lq > 0) x1 = x1 * Sqrt(gtoexp(jshl)**lq)
         do ipatt=1, npatt
           iao = iao + 1
           call LPatt(lq, ipatt, LPat, Lx, Ly, Lz)
           do ipt = 1, npt
             rpa = pnt6(:,ipt) - xyz(:,iam)
             x2 = One
             If(lq > 0) x2 = rpa(1)**Lx * rpa(2)**Ly * rpa(3)**Lz
             x3 = gtoexp(jshl) * dotx(3,rpa,rpa)
             aol(iao,ipt) = fnorm(noff+ipatt) * x1 * x2 * exp(-x3)
           end do
         end do
       end do
     end do
   end do

   call MOCar2Sph(natom,mxla,nshlla,npcar,npsph,npt,aol,scr(isc1),c2scoef)

   ! small component
   if(iham > 1) then
     iao = 0
     jshl = 0
     do iam = 1, natom
       do lq = 0, mxla(iam)
         noff = (lq+0)*(lq+1)*(lq+2)/6
         npatt= (lq+1)*(lq+2)/2
         do ishl = 1, nshlla(lq,iam)
           jshl = jshl + 1
           x1 = gtoexp(jshl)**(0.75d0)
           If(lq > 0) x1 = x1 * Sqrt(gtoexp(jshl)**lq)
           do ipatt=1, npatt
             iao = iao + 1
             call LPatt(lq, ipatt, LPat, Lx, Ly, Lz)
             do ipt = 1, npt
               rpa = pnt6(:,ipt) - xyz(:,iam)
               x2 = rpa(1)**(Lx+1) * rpa(2)**Ly * rpa(3)**Lz
               x2 = rpa(1)**Lx * rpa(2)**(Ly+1) * rpa(3)**Lz + x2
               x2 = rpa(1)**Lx * rpa(2)**Ly * rpa(3)**(Lz+1) + x2
               x2 = (gtoexp(jshl) + gtoexp(jshl)) * x2
               x3 = Zero
               if(Lx > 0) x3 = dble(Lx) * rpa(1)**(Lx-1) * rpa(2)**Ly * rpa(3)**Lz + x3
               if(Ly > 0) x3 = dble(Ly) * rpa(1)**Lx * rpa(2)**(Ly-1) * rpa(3)**Lz + x3
               if(Lz > 0) x3 = dble(Lz) * rpa(1)**Lx * rpa(2)**Ly * rpa(3)**(Lz-1) + x3
               x4 = gtoexp(jshl) * dotx(3,rpa,rpa)
               aos(iao,ipt) = fnorm(noff+ipatt) * x1 * (x3 - x2) * exp(-x4)
             end do
           end do
         end do
       end do
     end do

     call MOCar2Sph(natom,mxla,nshlla,npcar,npsph,npt,aos,scr(isc2),c2scoef)
   end if

   call MatMult(1,npsph,npsph,npt,One,Zero,pv,scr(isc1),scr(isc3))
   noff = 0
   do ipt = 1, npt
     scr(ipt) = dotx(npsph,scr(isc1+noff),scr(isc3+noff))
     noff = noff + npsph
   end do

   if(iham > 1) then
     call MatMult(1,npsph,npsph,npt,One,Zero,pw,scr(isc2),scr(isc3))
     noff = 0
     do ipt = 1, npt
       scr(ipt) = scr(ipt) + dotx(npsph,scr(isc2+noff),scr(isc3+noff)) / f_4cc
       noff = noff + npsph
     end do
   end if

   val_cd = scr(1)

   ! save curve data
   if(npt > 1) then
     x3 = sqrt(delx*delx + dely*dely + delz*delz) * au2ang
     x1 = -x3
     do ipt = 1, npt
       x1 = x1 + x3
       !<<< radial density
       x2 = scr(ipt)
       ! radial distribution function
       !x2 = fpi * x1 * x1 * scr(ipt)
       !>>>
       write(88,"(f12.7,f22.7)") x1, x2
     end do
   end if

   call ElemZA(1,Elm,iza(iatom),Elm)
   write(iout,"(i6,4x,a3,i5,f13.5,13x,f20.5)") iatom, Elm, iza(iatom), rms(iatom)*fac1, val_cd
 end do
 write(iout,"(1x,63('-'))")

 deallocate(aol, LPat, fnorm, pnt6, rpa)
 if(iham > 1) deallocate(aos)

 return
end subroutine cdcalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! dipole momemt contribution from nuclear charges
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine nucdip(natom,za,xyz, npchg,pchg,pxyz, dip)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension za(natom), xyz(3,natom), pchg(npchg), pxyz(3,npchg), dip(3)

 ! C is fixed at O
 do i = 1, natom + npchg
   if(i <= natom) then
     do j=1,3
       dip(j) = dip(j) + za(i) * xyz(j,i)
     end do
   else
     do j=1,3
       dip(j) = dip(j) + pchg(i-natom) * pxyz(j,i-natom)
     end do
   endif
 end do

 Return
End subroutine nucdip

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! relativistic dipole moment calculation
!
! Input:
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   irelc          : scalar or 2-component calculation
!   natom          : number of atoms
!   iham           : NR or X2C Hamiltonians
!   za             : atomic nuclear charges
!   xyz            : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   npcar          : number of primitive Cartesian functions
!   npsph          : number of primitive spherical functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   pv,pw          : transformed density matrices for V and W
!   npchg,pchg,pxyz: point charges
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dipcalc(iout,mpatt,irelc,natom,iham,za,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,pv,pw, &
 npchg,pchg,pxyz, scr,lscr)
 use Constant, only : maxlq, f_4cc, au2deb, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension     :: za(natom), xyz(3,natom), mxla(natom), nshlla(0:maxlq,natom), gtoexp(*), c2scoef(*),  &
                  pv(npsph*npsph,*), pw(npsph*npsph,*), pchg(*), pxyz(*), scr(lscr)
 allocatable   :: dip(:), dmat(:,:)

 ntt  = npcar*(npcar+1)/2
 nss  = npcar*npcar
 npss = npsph*npsph
 lmat = 3
 if(iham > 1 .and. irelc == 2) lmat = 12

 ! single-value density matrices
 if(irelc == 1) then
   lden = 0
 else if(irelc == 2) then
   ipv  = 1
   lden = npss
   ! Only the real part of Pv should be taken
   call SglMat(npsph,npsph,pv,scr(ipv))
   ! Split Pw into Pwo, Pwx, Pwy, and Pwz
   if(iham > 1) then
     ipwo = ipv  + npss
     ipwx = ipwo + npss
     ipwy = ipwx + npss
     ipwz = ipwy + npss
     lden = lden + npss*4
     call SGLCOMPW(npsph,pw(1,1),pw(1,5),scr(ipwo),scr(ipwx),scr(ipwy),scr(ipwz))
   end if
 end if

 isc1 = lden +   1
 isc2 = isc1 + nss
 isc3 = isc2 + nss
 isc4 = isc3 + nss

 allocate(dip(3),dmat(ntt,lmat))

 ! S (for test)
 ! call pdip_drv(iout,0,0,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,ntt,dmat,scr(isc1),lscr-lden)
 ! call lt2sqr(0,npcar,dmat(1,1),scr(isc1))
 ! call writemat(npcar,npcar,1,scr(isc1),30)

 dip = Zero
 ! contribution from nuclear charges
 call nucdip(natom,za,xyz, npchg,pchg,pxyz, dip)

 ! dV/dF
 call pdip_drv(iout,1,0,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,ntt,dmat,scr(isc1),lscr-lden)
 do imat = 1, 3
   call lt2sqr(0,npcar,dmat(1,imat),scr(isc1))
   ! call writemat(npcar,npcar,1,scr(isc1),30)
   ! dV/dF --> sc2
   call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
   if(irelc == 1) then
     dip(imat) = dip(imat) - traceSS(npsph,pv,scr(isc2))
   else if(irelc == 2) then
     dip(imat) = dip(imat) - traceSS(npsph,scr(ipv),scr(isc2))
   end if
 end do

 if(iham > 1) then
   if(irelc == 1) then
     ! dW0/dF
     call pdip_drv(iout,2,0,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,ntt,dmat,scr(isc1),lscr-lden)
     do imat = 1, 3
       call lt2sqr(0,npcar,dmat(1,imat),scr(isc1))
       ! call writemat(npcar,npcar,1,scr(isc1),30)
       ! dW0/dF --> sc2
       call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
       dip(imat) = dip(imat) - traceSS(npsph,pw,scr(isc2)) / f_4cc
     end do
   else if(irelc == 2) then
     ! dW/dF
     call pdip_drv(iout,3,0,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,ntt,dmat,scr(isc1),lscr-lden)
     do imat = 1, 3
       call lt2sqr(0,npcar,dmat(1,imat),scr(isc1))
       ! write(iout,*)'dWo/dF',imat
       ! call writemat(npcar,npcar,1,scr(isc1),30)
       ! dW0/dF --> sc2
       call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
       dip(imat) = dip(imat) - traceSS(npsph,scr(ipwo),scr(isc2)) / f_4cc
     end do
     do imat = 4, 6
       call lt2sqr(1,npcar,dmat(1,imat),scr(isc1))
       ! write(iout,*)'dWx/dF',imat-3
       ! call writemat(npcar,npcar,1,scr(isc1),30)
       ! dWx/dF --> sc2
       call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
       ! * -1 due to antisymmetry
       dip(imat-3) = dip(imat-3) + traceSS(npsph,scr(ipwx),scr(isc2)) / f_4cc
     end do
     do imat = 7, 9
       call lt2sqr(1,npcar,dmat(1,imat),scr(isc1))
       ! write(iout,*)'dWy/dF',imat-6
       ! call writemat(npcar,npcar,1,scr(isc1),30)
       ! dWy/dF --> sc2
       call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
       ! * -1 due to antisymmetry
       dip(imat-6) = dip(imat-6) + traceSS(npsph,scr(ipwy),scr(isc2)) / f_4cc
     end do
     do imat = 10, 12
       call lt2sqr(1,npcar,dmat(1,imat),scr(isc1))
       ! write(iout,*)'dWz/dF',imat-9
       ! call writemat(npcar,npcar,1,scr(isc1),30)
       ! dWz/dF --> sc2
       call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
       ! * -1 due to antisymmetry
       dip(imat-9) = dip(imat-9) + traceSS(npsph,scr(ipwz),scr(isc2)) / f_4cc
     end do
   end if
 end if

 do imat = 1, 3
   if(abs(dip(imat)) < 1.0d-8) dip(imat) = Zero
 end do

 write(iout,"(/,4x,'|',2(f12.5,','),f12.5,'|  =',2x,f12.5,'   (a.u.)')") dip, sqrt(dotx(3,dip,dip))
 write(iout,"(/,4x,'|',2(f12.5,','),f12.5,'|  =',2x,f12.5,'   (Debye)')") dip*au2deb, sqrt(dotx(3,dip,dip))*au2deb

 deallocate(dip,dmat)

 return
end subroutine dipcalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! atomic polar tensor (APT) calculation (NR Hamiltonian only)
!
! Input:
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   irelc          : scalar or 2-component calculation
!   natom          : number of atoms
!   iham           : NR or X2C Hamiltonians
!   nchag          : net charge
!   iza,za         : atomic nuclear charges
!   xyz            : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   npcar          : number of primitive Cartesian functions
!   npsph          : number of primitive spherical functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   pv,pw          : transformed density matrices for V and W
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine aptcalc(iout,mpatt,irelc,natom,iham,nchag,iza,za,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,smat,  &
 pv,pw,scr,lscr)
 use Constant, only : maxlq, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension     :: iza(natom), za(natom), xyz(3,natom), mxla(natom), nshlla(0:maxlq,natom), gtoexp(*), c2scoef(*),  &
                  smat(npsph*npsph), pv(npsph*npsph,*), pw(npsph*npsph,*), scr(lscr)
 character*3   :: Elm
 character*1   :: cha1(3), l2u
 allocatable   :: apt(:,:,:), dmat(:,:,:), nabc(:,:)

 ! call writemat(npcar,npsph,1,pv,30)
 cha1 = (/'x','y','z'/)
 ntt  = npcar*(npcar+1)/2
 nss  = npcar*npcar
 npss = npsph*npsph

 ! single-value density matrices
 if(irelc == 1) then
   lden = 0
 else if(irelc == 2) then
   ipv  = 1
   lden = npss
   ! Only the real part of Pv should be taken
   call SglMat(npsph,npsph,pv,scr(ipv))
   ! Split Pw into Pwo, Pwx, Pwy, and Pwz
   if(iham > 1) then
     write(iout,"(/,' Error! APT has been programmed for the NR Hamiltonian only.',/)")
     call estop(1)
     ! ipwo = ipv  + npss
     ! ipwx = ipwo + npss
     ! ipwy = ipwx + npss
     ! ipwz = ipwy + npss
     ! lden = lden + npss*4
     ! call SGLCOMPW(npsph,pw(1,1),pw(1,5),scr(ipwo),scr(ipwx),scr(ipwy),scr(ipwz))
   end if
 end if

 isc1 = lden +   1
 isc2 = isc1 + nss
 isc3 = isc2 + nss
 isc4 = isc3 + nss

 allocate(apt(3,3,natom),dmat(ntt,3,3),nabc(2,natom))

 apt = Zero

 ! S
 ! call writemat(npsph,npsph,1,smat,30)

 ! dS/dX (for test)
 ! call pdip_drv(iout,0,1,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,ntt,dmat,scr(isc1),lscr-lden)
 ! do imat = 1, 3
 !   call lt2sqr(0,npcar,dmat(1,imat),scr(isc1))
 !   call writemat(npcar,npcar,1,scr(isc1),30)
 ! end do

 nabc(1,1) = 1
 do imat = 1, natom
   if(imat > 1) nabc(1,imat) = nabc(1,imat-1) + nabc(2,imat-1)
   nabc(2,imat) = 0
   do i = 0, mxla(imat)
     nabc(2,imat) = nabc(2,imat) + (i+i+1)*nshlla(i,imat)
   end do
 end do

 ! contribution from nuclear charges
 do imat = 1, natom
   do i = 1, 3
     apt(i,i,imat) = za(imat)
   end do
 end do

 ! ddV/dFdX
 call pdip_drv(iout,1,1,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,ntt,dmat,scr(isc1),lscr-lden)
 
 do imat = 1, 3      ! dX, dY, dZ
   do i = 1, 3       ! F_x, F_y, F_z
     call lt2sqr(0,npcar,dmat(1,i,imat),scr(isc1))
     !if(i==imat) then
     !  write(*,"('dF=',i2,', dX=',i2)") i, imat
     !  call writemat(npcar,npcar,1,scr(isc1),30)
     !end if
     ! d^2V/dFdX --> sc2
     call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
     if(irelc == 1) then
       call STAPT(natom,i,imat,npsph,apt,nabc,pv,smat,scr(isc2))
     else if(irelc == 2) then
       call STAPT(natom,i,imat,npsph,apt,nabc,scr(ipv),smat,scr(isc2))
     end if
   end do
 end do

 if(iham > 1) then
   write(iout,"(/,' Error! APT has been programmed for the NR Hamiltonian only.',/)")
   call estop(1)
 end if

 ! shift
 do i = 1, 3      ! dX, dY, dZ
   do j = 1, 3    ! F_x, F_y, F_z
     scr(isc1) = Zero
     do imat = 1, natom
       scr(isc1) = scr(isc1) + apt(j,i,imat)
     end do
     if(i == j) then
       scr(isc1) = (scr(isc1) - dble(nchag)) / dble(natom)
     else
       scr(isc1) = scr(isc1) / dble(natom)
     end if
     do imat = 1, natom
       apt(j,i,imat) = apt(j,i,imat) - scr(isc1)
       if(abs(apt(j,i,imat)) < 1.0d-8) apt(j,i,imat) = Zero
     end do
   end do
 end do

 do imat = 1, natom
   call ElemZA(1,Elm,iza(imat),Elm)
   write(iout,"(/,i4,3x,a3,3x,3(11x,'dF_',a1),6x,'(a.u.)',/,1x,69('-'))") imat, Elm, cha1(1:3)
   do i = 1, 3
     write(iout,"(12x,'d',a1,4x,3f15.7)") l2u(cha1(i)), apt(:,i,imat)
   end do
 end do

 deallocate(apt,dmat,nabc)

 return
end subroutine aptcalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! atomic DDD charges (DDD) calculation (NR Hamiltonian only)
!
! Input:
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   irelc          : scalar or 2-component calculation
!   natom          : number of atoms
!   iham           : NR or X2C Hamiltonians
!   nchag          : net charge
!   iza,za         : atomic nuclear charges
!   xyz            : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   npcar          : number of primitive Cartesian functions
!   npsph          : number of primitive spherical functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   pv,pw          : transformed density matrices for V and W
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dddcalc(iout,mpatt,irelc,natom,iham,nchag,iza,za,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,smat,  &
 pv,pw,scr,lscr)
 use Constant, only : maxlq, Zero, au2deb
 implicit real(kind=8) (a-h,o-z)
 dimension     :: iza(natom), za(natom), xyz(3,natom), mxla(natom), nshlla(0:maxlq,natom), gtoexp(*), c2scoef(*),  &
                  smat(npsph*npsph), pv(npsph*npsph,*), pw(npsph*npsph,*), scr(lscr)
 character*3   :: Elm
 character*1   :: l2u
 allocatable   :: dip(:), ddd(:), dmat(:), nabc(:,:)

 ! call writemat(npsph,npsph,1,pv,30)
 ntt  = npcar*(npcar+1)/2
 nss  = npcar*npcar
 npss = npsph*npsph

 ! single-value density matrices
 if(irelc == 1) then
   lden = 0
 else if(irelc == 2) then
   ipv  = 1
   lden = npss
   ! Only the real part of Pv should be taken
   call SglMat(npsph,npsph,pv,scr(ipv))
   ! Split Pw into Pwo, Pwx, Pwy, and Pwz
   if(iham > 1) then
     write(iout,"(/,' Error! DDD has been programmed for the NR Hamiltonian only.',/)")
     call estop(1)
     ! ipwo = ipv  + npss
     ! ipwx = ipwo + npss
     ! ipwy = ipwx + npss
     ! ipwz = ipwy + npss
     ! lden = lden + npss*4
     ! call SGLCOMPW(npsph,pw(1,1),pw(1,5),scr(ipwo),scr(ipwx),scr(ipwy),scr(ipwz))
   end if
 end if

 isc1 = lden +   1
 isc2 = isc1 + nss
 isc3 = isc2 + nss
 isc4 = isc3 + nss

 allocate(dip(3), ddd(natom),dmat(ntt),nabc(2,natom))

 ! contribution from nuclear charges
 ddd = za

 ! S
 ! call writemat(npsph,npsph,1,smat,30)

 nabc(1,1) = 1
 do imat = 1, natom
   if(imat > 1) nabc(1,imat) = nabc(1,imat-1) + nabc(2,imat-1)
   nabc(2,imat) = 0
   do i = 0, mxla(imat)
     nabc(2,imat) = nabc(2,imat) + (i+i+1)*nshlla(i,imat)
   end do
 end do

 ! S in npcar --> dmat. dip is used for dumb
 call pdst_drv(iout,0,0,.false.,mpatt,natom,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,dmat,dip,dip,scr(isc1),lscr-lden)
 ! ddV/dFdX --> scr1
 call ddd_drv(natom,nshell,mxla,nshlla,gtoexp,ntt,scr(isc1))
 ! scr1 * dmat
 call scamul(ntt,scr(isc1),dmat,dmat)
 call lt2sqr(0,npcar,dmat,scr(isc1))
 ! call writemat(npcar,npcar,1,scr(isc1),30)
 call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
 if(irelc == 1) then
   call STDDD(natom,npsph,ddd,nabc,pv,smat,scr(isc2))
 else if(irelc == 2) then
   call STDDD(natom,npsph,ddd,nabc,scr(ipv),smat,scr(isc2))
 end if

 if(iham > 1) then
   write(iout,"(/,' Error! DDD has been programmed for the NR Hamiltonian only.',/)")
   call estop(1)
 end if

 ! shift
 scr(isc1) = Zero
 do imat = 1, natom
   scr(isc1) = scr(isc1) + ddd(imat)
 end do
 scr(isc1) = (scr(isc1) - dble(nchag)) / dble(natom)
 do imat = 1, natom
   ddd(imat) = ddd(imat) - scr(isc1)
   if(abs(ddd(imat)) < 1.0d-8) ddd(imat) = Zero
 end do

 write(iout,*)
 do imat = 1, natom
   call ElemZA(1,Elm,iza(imat),Elm)
   write(iout,"(i4,3x,a3,f15.4)") imat, Elm, ddd(imat)
 end do

 dip = Zero
 do i = 1, 3
   do imat = 1, natom
     dip(i) = dip(i) + xyz(i,imat) * ddd(imat)
   end do
 end do

 write(iout,"(/,1x,'electric dipole moment from DDD charges')")
 write(iout,"(/,4x,'|',2(f12.5,','),f12.5,'|  =',2x,f12.5,'   (a.u.)')") dip, sqrt(dotx(3,dip,dip))
 write(iout,"(/,4x,'|',2(f12.5,','),f12.5,'|  =',2x,f12.5,'   (Debye)')") dip*au2deb, sqrt(dotx(3,dip,dip))*au2deb

 deallocate(dip,ddd,dmat,nabc)

 return
end subroutine dddcalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Electric field gradient (EFG) calculation (local version; not suggested!)
!
! Input:
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   irelc          : scalar or 2-component calculation
!   natom          : number of atoms
!   iham           : NR or X2C Hamiltonians
!   minza          : min. atomic nuclear charge for EFG calculation
!   loca           : 1CA approximation (see the manual)
!   iza,za         : atomic nuclear charges
!   rms            : nuclear RMS charge radii
!   xyz            : Cartesian coordinates
!   npsph          : number of primitive spherical functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   ialbs          : starting position of each L in primitive spherical basis functions
!   iptexp         : starting positions of GTO exponents for each L
!   gtoexp         : GTO exponents
!   pv,pw          : transformed density matrices for V and W
!   npchg          : number of point charge
!   pchg, pxyz     : point charges and their coordinates
!   qdat           : user's NQM value qdat(2,i) if qdat(1,i) = 1.0
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine efgcalclo(iout,mpatt,irelc,natom,iham,minza,loca,iza,za,rms,xyz,npsph,mxla,nshlla,ialbs,iptexp,gtoexp,c2scoef,  &
 pv,pw, npchg,pchg,pxyz,qdat, scr,lscr)
 use Constant, only : maxlq, au2wvn, clight_si, au2ang, Zero, f_4cc
 implicit real*8(A-H,O-Z)
 dimension     :: iza(natom), za(natom), rms(natom), xyz(3,natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom),  &
                  iptexp(0:maxlq,natom), gtoexp(*), c2scoef(*),  &
                  pv(npsph*npsph,*), pw(npsph*npsph,*), pchg(*), pxyz(3,*), qdat(2,natom), scr(lscr)
 allocatable   :: dvm(:,:), efglt(:,:), efgsq(:,:), efgeig(:), d1(:,:), d2(:), d5(:), xyz0(:), iziso(:)
 character*3   :: Elm

 if(iham == 1) then
   if(loca > 1) write(iout,"(/,' The Atomic block approximation is turned on in P.')")
 else
   if(loca == 2) then
     write(iout,"(/,' The Atomic block approximation is turned on in Pv and Pw.')")
   else if(loca == 3) then
     write(iout,"(/,' The Atomic block approximation is turned on in P.')")
   end if
 end if

 fac1 = sqrt(1.5d0) * 1.0d5 * au2ang
 fac2 = au2wvn * clight_si * 1.0d-14 / (au2ang * au2ang)
 fac2 = fac2 * 1.0d-3   ! millibarn to barn
 npss = npsph*npsph

 allocate(efglt(6,3), efgsq(3,3), efgeig(3), d1(3,natom+npchg), d2(natom+npchg), d5(natom+npchg), xyz0(3), iziso(3))
 xyz0 = Zero

 ! check minza
 call ChkMinZa(iout,iham,natom,iza,minza,6,ncalc)
 if(ncalc == 0) return

 ! single-value density matrices
 if(irelc == 1) then
   lden = 0
 else if(irelc == 2) then
   ipv  = 1
   lden = npss
   ! Only the real part of Pv should be taken
   call SglMat(npsph,npsph,pv,scr(ipv))
   ! Split Pw into Pwo, Pwx, Pwy, and Pwz
   if(iham > 1) then
     lden = npss*5
     ipwo = ipv  + npss
     ipwx = ipwo + npss
     ipwy = ipwx + npss
     ipwz = ipwy + npss
     call SGLCOMPW(npsph,pw(1,1),pw(1,5),scr(ipwo),scr(ipwx),scr(ipwy),scr(ipwz))
   end if
 end if

 write(iout,"(1x,106('-'),/,'   No.   Atom   ZA     RMS (fm)',10x,'EFG eigen. (a.u.)',23x,'EFG tensor (a.u.)',/,1x,106('-'))")
 isc1 = lden +   1

 do iatom = 1, natom
   if(minza < 0) then
     if(iza(iatom) /= abs(minza)) cycle
   else
     if(iza(iatom) < minza) cycle
   end if

   nshell0 = 0
   I0 = ialbs(0,iatom)
   Ic = 0
   Ip = 0
   do i = 0, mxla(iatom)
     ngto = nshlla(i,iatom)
     nshell0  = nshell0  + ngto
     Ic = Ic + (i+1)*(i+2)*ngto/2
     Ip = Ip + (i+i+1)*ngto
   end do
   ntt = Ic*(Ic+1)/2
   nss = Ic*Ic
   isc2 = isc1 + nss
   isc3 = isc2 + nss
   isc4 = isc3 + nss

   allocate(dvm(ntt,6))
   efglt = Zero
   call pnai_drv(iout,0,0,0,-1,mpatt,1,za(iatom),xyz0,mxla(iatom),nshell0,Ic,mxla(iatom),nshlla(0,iatom),  &
     gtoexp(iptexp(0,iatom)),rms(iatom),dvm,scr(isc1),lscr-lden)
   do ixy = 1, 6
     call lt2sqr(0,Ic,dvm(1,ixy),scr(isc1))
     ! dV/dQ --> sc2
     call car2sph(0,1,mxla(iatom),nshlla(0,iatom),Ic,Ip,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
     if(irelc == 1) then
       call TakDBlk(npsph,I0,Ip,pv,scr(isc3))
     else if(irelc == 2) then
       call TakDBlk(npsph,I0,Ip,scr(ipv),scr(isc3))
     end if
     efglt(ixy,1) = traceSS(Ip,scr(isc3),scr(isc2))
   end do

   if(iham > 1) then
     if(irelc == 1) then
       call pnai_drv(iout,irelc,0,0,-1,mpatt,1,za(iatom),xyz0,mxla(iatom),nshell0,Ic,mxla(iatom),nshlla(0,iatom),  &
         gtoexp(iptexp(0,iatom)),rms(iatom),dvm,scr(isc1),lscr-lden)
       do ixy = 1, 6
         call lt2sqr(0,Ic,dvm(1,ixy),scr(isc1))
         ! dW/dQ --> sc2
         call car2sph(0,1,mxla(iatom),nshlla(0,iatom),Ic,Ip,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
         call TakDBlk(npsph,I0,Ip,pw,scr(isc3))
         efglt(ixy,1) = efglt(ixy,1) + traceSS(Ip,scr(isc3),scr(isc2)) / f_4cc
       end do
     else if(irelc == 2) then
       jpwso = ipwo
       do ipvp = 1, 4
         call pnai_drv(iout,irelc,ipvp,0,-1,mpatt,1,za(iatom),xyz0,mxla(iatom),nshell0,Ic,mxla(iatom),nshlla(0,iatom),  &
           gtoexp(iptexp(0,iatom)),rms(iatom),dvm,scr(isc1),lscr-lden)
         do ixy = 1, 6
           call lt2sqr((ipvp+1)/3,Ic,dvm(1,ixy),scr(isc1))
           ! dW/dQ --> sc2
           call car2sph(0,1,mxla(iatom),nshlla(0,iatom),Ic,Ip,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
           call TakDBlk(npsph,I0,Ip,scr(jpwso),scr(isc3))
           if(ipvp == 1) then
             efglt(ixy,1) = efglt(ixy,1) + traceSS(Ip,scr(isc3),scr(isc2)) / f_4cc
           else
             ! * -1 due to antisymmetry
             efglt(ixy,1) = efglt(ixy,1) - traceSS(Ip,scr(isc3),scr(isc2)) / f_4cc
           end if
         end do
         jpwso = jpwso + npss
       end do
     end if
   end if

   call EFGNuc(iatom,efglt(1,2), natom,za,xyz, npchg,pchg,pxyz, d1,d2,d5)

   efglt(:,3) = efglt(:,1) + efglt(:,2)
   call LT2Sqr(0,3,efglt(1,3),efgsq)
   efgeig = Zero
   scr(isc1:isc1+8) = Zero
   call dsyev('N','L',3,efgsq,3,efgeig,scr(isc1),9,INFO)
   if(INFO /= 0) then
     write(*,"(/,' Diagnolization failed in efgcalclo.')")
     call estop(1)
   end if
   ! reorder 3 eigenvalues |v1| <= |v2| <= |v3|
   call reord3abs(efgeig)
   ! load EFG tensor again since efgsq has been destroyed
   call LT2Sqr(0,3,efglt(1,3),efgsq)

   call ElemZA(1,Elm,iza(iatom),Elm)
   write(iout,"(i6,4x,a3,i5,f13.5,9x,'Vaa =',f13.5,7x,3f14.5)") iatom, Elm, iza(iatom), rms(iatom)*fac1, efgeig(1), efgsq(:,1)
   write(iout,"(40x,'Vbb =',f13.5,7x,3f14.5)") efgeig(2), efgsq(:,2)
   write(iout,"(40x,'Vcc =',f13.5,7x,3f14.5)") efgeig(3), efgsq(:,3)
   eta = 0.0d0
   if(abs(efgeig(3)) > 1.0d-8) eta = abs( (efgeig(1) - efgeig(2)) / efgeig(3) )
   write(iout,"(/,40x,'eta =',f13.5)") eta

   ! correction for open-shell atom & linear molecule
   ! d1 and scr are used for scratch
   if(npchg == 0) call efgfix(iout,natom,xyz,eta,efgeig,efgsq, d1,scr(isc1))

   if(nint(qdat(1,iatom)) == 0) then
     call NQMlib(iza(iatom),niso,iziso,scr(isc1))
   else
     niso = 1
     iziso(1) = 0
     scr(isc1) = qdat(2,iatom)
   end if
   if(niso > 0) then
     scr(isc1+niso) = fac2*efgeig(3)*scr(isc1)
     write(iout,"(/,40x,'eQq =',f13.5,' MHz  with  Q(',a3,'-',i3,') =',f9.2,' millibarn')") scr(isc1+niso),  &
       Elm,iziso(1),scr(isc1)
     ! dEq for 57Fe
     if(iza(iatom) == 26 .and. iziso(1) == 57) then
       scr(isc1+niso+1) = 0.5d0 * scr(isc1+niso) * sqrt(1.0d0 + eta*eta/3.0d0)
       write(iout,"(40x,'dEq =',f13.5,' MHz  or',f13.5,' mm/s for (57)Fe')") scr(isc1+niso+1), scr(isc1+niso+1)/11.6248d0 
     end if
     do i = 2, niso  ! niso <= 3
       scr(isc1+niso) = fac2*efgeig(3)*scr(isc1+i-1)
       write(iout,"(45x,f13.5,14x,a3,'-'i3,3x,f9.2)") scr(isc1+niso), Elm, iziso(i), scr(isc1+i-1)
     end do
   else
     write(iout,"(/,40x,'eQq :      N.A.')")
   end if

   call LT2Sqr(0,3,efglt(1,2),efgsq)
   write(iout,"(/,40x,'Nuclear contributions    ',3f14.5,2(/,65x,3f14.5))") efgsq
   call LT2Sqr(0,3,efglt(1,1),efgsq)
   write(iout,"(/,40x,'Electronic contributions ',3f14.5,2(/,65x,3f14.5))") efgsq
   write(iout,*)

   deallocate(dvm)

 end do
 write(iout,"(1x,106('-'))")

 deallocate(efglt, efgsq, efgeig, d1, d2, d5, xyz0, iziso)

 return
end subroutine efgcalclo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Electric field gradient (EFG) calculation
!
! Input:
!   mpatt          : pattern mode of basis functions. 0: Gaussian, -1: Molden
!   irelc          : scalar or 2-component calculation
!   natom          : number of atoms
!   iham           : NR or X2C Hamiltonians
!   minza          : min. atomic nuclear charge for EFG calculation
!   ipop           : population if =1. Scalar calculation only.
!   iza,za         : atomic nuclear charges
!   rms            : nuclear RMS charge radii
!   xyz            : Cartesian coordinates
!   maxtyp         : max_L of the system
!   nshell         : number of primitive shells
!   npcar          : number of primitive Cartesian functions
!   npsph          : number of primitive spherical functions
!   mxla           : max_L for each atom
!   nshlla         : numbers of primitive shells for each atom
!   gtoexp         : GTO exponents
!   pv,pw          : transformed density matrices for V and W
!   npchg          : number of point charge
!   pchg, pxyz     : point charges and their coordinates
!   qdat           : user's NQM values
!   nmo,occ,ene,ispn,symm
!                  : information of MOs
!   cr,cu          : transformed MOs
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine efgcalc(iout,mpatt,irelc,natom,iham,minza,ipop,iza,za,rms,xyz,maxtyp,nshell,npcar,npsph,mxla,nshlla,gtoexp,c2scoef,  &
 pv,pw, npchg,pchg,pxyz,qdat, nmo,occ,ene,ispn,symm,cr,cu, scr,lscr)
 use Constant, only : t33nam, maxlq, au2wvn, clight_si, au2ang, kev2mhz, Zero, One, f_4cc, occtol
 implicit real*8(A-H,O-Z)
 parameter(maxiso=3)
 character*10  :: symm
 dimension     :: iza(natom), za(natom), rms(natom), xyz(3,natom), mxla(natom), nshlla(0:maxlq,natom), gtoexp(*), c2scoef(*),  &
                  pv(npsph*npsph,*), pw(npsph*npsph,*), pchg(*), pxyz(3,*), qdat(2,natom), &
                  occ(nmo), ene(nmo), ispn(nmo), symm(nmo), cr(npsph*irelc*nmo,irelc), cu(npsph*irelc*nmo,irelc), scr(lscr)
 allocatable   :: dvm(:,:), efglt(:,:), efgsq(:,:), efgeig(:), d1(:,:), d2(:), d5(:)
 allocatable   :: I0(:),IW(:),popefg(:,:,:)
 allocatable   :: iziso(:), is2(:), iex(:), egm(:), qvl(:)
 character*3   :: Elm, Cst(2)

 fac1 = sqrt(1.5d0) * 1.0d5 * au2ang
!  NOTE. to eQq in MHz: CB * eQ (in barn) * Vqq(a.u.), where
!    CB = 2Ry * c * 1.0e-28 / a0^2
!       = 219474.63137 * 2.99792458e10 * 1.0e-14 / (0.5291772083^2) = 234.96478
!  See p. 406 in Computer Physics Communications 74, 399 (1993).
 fac2 = au2wvn * clight_si * 1.0d-14 / (au2ang * au2ang)
 fac2 = fac2 * 1.0d-3   ! * millibarn to barn
 fac3 = clight_si * 1.0d-10 / kev2mhz
 ntt  = npcar*(npcar+1)/2
 nss  = npcar*npcar
 npss = npsph*npsph

 Cst = (/'gr.','ex.'/)

 ! single-value density matrices
 if(irelc == 1) then
   ipop0 = ipop
   lden = 0
 else if(irelc == 2) then
   ipop0 = 0
   ipv  = 1
   lden = npss
   ! Only the real part of Pv should be taken
   call SglMat(npsph,npsph,pv,scr(ipv))
   ! Split Pw into Pwo, Pwx, Pwy, and Pwz
   if(iham > 1) then
     ipwo = ipv  + npss
     ipwx = ipwo + npss
     ipwy = ipwx + npss
     ipwz = ipwy + npss
     lden = lden + npss*4
     call SGLCOMPW(npsph,pw(1,1),pw(1,5),scr(ipwo),scr(ipwx),scr(ipwy),scr(ipwz))
   end if
 end if

 isc1 = lden +   1
 isc2 = isc1 + nss
 isc3 = isc2 + nss
 isc4 = isc3 + nss

 ! check minza
 call ChkMinZa(iout,iham,natom,iza,minza,6,ncalc)
 if(ncalc == 0) return

 ! starting position and width of primitive sperical functions for each atom
 if(ipop0 == 1) then
   allocate(I0(natom), IW(natom), popefg(natom,nmo,6))

   do iatom = 1, natom
     if(iatom == 1) then
       I0(iatom) = 1
     else
       I0(iatom) = I0(iatom-1) + Iw(iatom-1)
     end if
     Iw(iatom) = 0
     do i = 0, mxla(iatom)
       Iw(iatom) = Iw(iatom) + (i+i+1)*nshlla(i,iatom)
     end do
   end do
   matom = min(4,natom)

   write(iout,"(1x,132('-'),/,'   No.   Atom   ZA     RMS (fm)',10x,'EFG eigen. (a.u.)',23x,'EFG tensor (a.u.)',/,1x,132('-'))")
 else
   write(iout,"(1x,106('-'),/,'   No.   Atom   ZA     RMS (fm)',10x,'EFG eigen. (a.u.)',23x,'EFG tensor (a.u.)',/,1x,106('-'))")
 end if

 allocate(dvm(ntt,6), efglt(6,3), efgsq(3,3), efgeig(3), d1(3,natom+npchg), d2(natom+npchg), d5(natom+npchg),  &
   iziso(maxiso), is2(maxiso), iex(maxiso), egm(maxiso), qvl(maxiso))

 do iatom = 1, natom
   if(minza < 0) then
     if(iza(iatom) /= abs(minza)) cycle
   else
     if(iza(iatom) < minza) cycle
   end if

   efglt = Zero
   if(ipop0 == 1) popefg = Zero

   call pnai_drv(iout,0,0,0,-iatom,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,dvm,scr(isc1),lscr-lden)
   do ixy = 1, 6
     call lt2sqr(0,npcar,dvm(1,ixy),scr(isc1))
     ! dV/dQ --> sc2
     call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
     if(irelc == 1) then
       efglt(ixy,1) = traceSS(npsph,pv,scr(isc2))
       ! population
       if(ipop0 == 1) call cdpop(npsph,nmo,natom,One,I0,Iw,scr(isc2),occ,cr,popefg(1,1,ixy),scr(isc3))
       !if(ipop0 == 1 .and. ixy == 6) call cdpopx(npsph,nmo,One,scr(isc2),occ,cr,scr(isc3),scr(isc4))
     else if(irelc == 2) then
       efglt(ixy,1) = traceSS(npsph,scr(ipv),scr(isc2))
     end if
   end do

   if(iham > 1) then
     if(irelc == 1) then
       call pnai_drv(iout,irelc,0,0,-iatom,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,dvm,scr(isc1),lscr-lden)
       do ixy = 1, 6
         call lt2sqr(0,npcar,dvm(1,ixy),scr(isc1))
         ! dW/dQ --> sc2
         call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
         efglt(ixy,1) = efglt(ixy,1) + traceSS(npsph,pw,scr(isc2)) / f_4cc
         ! population
         if(ipop0 == 1) call cdpop(npsph,nmo,natom,(One/f_4cc),I0,Iw,scr(isc2),occ,cu,popefg(1,1,ixy),scr(isc3))
         !if(ipop0 == 1 .and. ixy == 6) call cdpopx(npsph,nmo,(One/f_4cc),scr(isc2),occ,cu,scr(isc3),scr(isc4))
       end do
     else if(irelc == 2) then
       jpwso = ipwo
       do ipvp = 1, 4
         call pnai_drv(iout,irelc,ipvp,0,-iatom,mpatt,natom,za,xyz,maxtyp,nshell,npcar,mxla,nshlla,gtoexp,rms,dvm,  &
           scr(isc1),lscr-lden)
         do ixy = 1, 6
           call lt2sqr((ipvp+1)/3,npcar,dvm(1,ixy),scr(isc1))
           ! dW/dQ --> sc2
           call car2sph(0,natom,mxla,nshlla,npcar,npsph,scr(isc1),scr(isc2),c2scoef,scr(isc3),scr(isc4))
           if(ipvp == 1) then
             efglt(ixy,1) = efglt(ixy,1) + traceSS(npsph,scr(jpwso),scr(isc2)) / f_4cc
           else
             ! * -1 due to antisymmetry
             efglt(ixy,1) = efglt(ixy,1) - traceSS(npsph,scr(jpwso),scr(isc2)) / f_4cc
           end if
         end do
         jpwso = jpwso + npss
       end do
     end if

   end if

   call EFGNuc(iatom,efglt(1,2), natom,za,xyz, npchg,pchg,pxyz, d1,d2,d5)

   efglt(:,3) = efglt(:,1) + efglt(:,2)
   call LT2Sqr(0,3,efglt(1,3),efgsq)
   efgeig = Zero
   scr(isc1:isc1+8) = Zero
   call dsyev('N','L',3,efgsq,3,efgeig,scr(isc1),9,INFO)
   if(INFO /= 0) then
     write(*,"(/,' Diagnolization failed in efgcalc.')")
     call estop(1)
   end if
   ! reorder 3 eigenvalues |v1| <= |v2| <= |v3|
   call reord3abs(efgeig)
   ! load EFG tensor again since efgsq has been destroyed
   call LT2Sqr(0,3,efglt(1,3),efgsq)

   call ElemZA(1,Elm,iza(iatom),Elm)
   write(iout,"(i6,4x,a3,i5,f13.5,9x,'Vaa =',f13.5,7x,3f14.5)") iatom, Elm, iza(iatom), rms(iatom)*fac1, efgeig(1), efgsq(:,1)
   write(iout,"(40x,'Vbb =',f13.5,7x,3f14.5)") efgeig(2), efgsq(:,2)
   write(iout,"(40x,'Vcc =',f13.5,7x,3f14.5)") efgeig(3), efgsq(:,3)
   eta = 0.0d0
   if(abs(efgeig(3)) > 1.0d-8) eta = abs( (efgeig(1) - efgeig(2)) / efgeig(3) )
   write(iout,"(/,40x,'eta =',f13.5)") eta

   ! correction for open-shell atom & linear molecule
   ! d1 and scr are used for scratch
   if(npchg == 0) call efgfix(iout,natom,xyz,eta,efgeig,efgsq, d1,scr(isc1))

   call LT2Sqr(0,3,efglt(1,2),efgsq)
   write(iout,"(/,40x,'Nuclear contributions    ',3f14.5,2(/,65x,3f14.5))") efgsq
   call LT2Sqr(0,3,efglt(1,1),efgsq)
   write(iout,"(/,40x,'Electronic contributions ',3f14.5,2(/,65x,3f14.5))") efgsq
   write(iout,*)

   write(iout,"(26x,'<< Nuclear Quadrupole Coupling Constant >>')")
   if(nint(qdat(1,iatom)) == 0) then
     call NQMlib(iza(iatom),niso,iziso,qvl)
   else
     niso = 1
     iziso(1) = 0
     qvl(1) = qdat(2,iatom)
   end if
   if(niso > 0) then
     do i = 1, niso  ! niso <= maxiso
       scr(isc1) = fac2 * efgeig(3) * qvl(i)
       write(iout,"(29x,a3,'-',i3,':  eQq =',f13.5,' MHz  with  Q =',f9.2,' millibarn')") Elm, iziso(i), scr(isc1), qvl(i)
     end do
   else
     write(iout,"(29x,'N.A.')")
   end if

   write(iout,"(/,26x,'<< Mossbauer Nuclear Quadrupole Splitting >>')")
   call Mosslib(iza(iatom),niso,iziso,is2,iex,egm,qvl)
   iso0 = 0
   if(niso > 0) then
     do i = 1, niso  ! niso <= maxiso
       iso0 = iso0 + 1
       scr(isc1) = fac2 * dRval_dEq(is2(i)) * efgeig(3) * qvl(i)
       ! I = 3/2 only
       if(is2(i) == 3) scr(isc1) = scr(isc1) * sqrt(1.0d0 + eta*eta/3.0d0)
       if(mod(is2(i),2) == 0) then
         write(iout,"(29x,a3,'-',i3,':  dEq =',f13.5,' MHz  with  Q =',f9.2,' millibarn & I =',i3,2x,'(',a3,')')")  &
           Elm, iziso(i), scr(isc1), qvl(i), is2(i)/2, Cst(iex(i)+1)
       else
         write(iout,"(29x,a3,'-',i3,':  dEq =',f13.5,' MHz  with  Q =',f9.2,' millibarn & I =',i2,'/2 (',a3,')')")  &
           Elm, iziso(i), scr(isc1), qvl(i), is2(i), Cst(iex(i)+1)
       end if
       write(iout,"(43x,'=',f13.5,' mm/s with Eg =',f9.2,' KeV')") fac3*scr(isc1)/egm(i), egm(i)
     end do
   end if
   if(iso0 == 0) then
     write(iout,"(29x,'N.A.')")
   else
     write(iout,"(/,29x,'NOTE. The term sqrt(1 + eta * eta / 3) is included in dEq only for I = 3/2.')")
   end if

   if(ipop0 == 1) then
     write(iout,"(/,11x,'#MO      Irrep  S',10x,'Ene',8x,'Occ',12x,'dV/dQ sum',10x,'dV/dQ population on atoms')")
     do imo = 1, nmo
       if(abs(occ(imo)) < occtol) cycle
       scr(isc1) = Zero
       do ixy = 1, 6
         scr(isc1+ixy) = sum(popefg(:,imo,ixy))
         scr(isc1) = scr(isc1) + abs(scr(isc1+ixy))
       end do
       if(scr(isc1) < 1.0d-5) cycle

       write(iout,"(8x,i6,1x,a10,i3,f13.4,f11.6)") imo, adjustr(symm(imo)), ispn(imo), ene(imo), occ(imo)
       do ixy = 1, 6
         write(iout,"(54x,a2,f17.4,4f15.4)") t33nam(ixy), scr(isc1+ixy), popefg(1:matom,imo,ixy)
         if(natom > 4) write(iout,"(73x,4f15.4)") popefg(matom+1:natom,imo,ixy)
       end do
     end do
     write(iout,*)
   end if

 end do

 deallocate(dvm, efglt, efgsq, efgeig, d1, d2, d5, iziso, is2, iex, egm, qvl)

 if(ipop0 == 1) then
   write(iout,"(1x,132('-'))")
   deallocate(I0, IW, popefg)
 else
   write(iout,"(1x,106('-'))")
 end if

 return
end subroutine efgcalc

!subroutine cdpopx(npsph,nmo,fac,dnai,occ,cmoL,cmoR,scr)
! use Constant, only : occtol, One, Zero
! implicit real(kind=8) (a-h,o-z)
! dimension         :: dnai(npsph,npsph), occ(nmo), cmoL(npsph,nmo), cmoR(npsph,nmo), scr(npsph,nmo)
!
! call MatMult(1,npsph,npsph,nmo,One,Zero,dnai,cmoL,cmoR)
!
! scr = Zero
! do imo = 1, nmo
!   if(abs(occ(imo)) < occtol) cycle
!
!   do ibs = 1, npsph
!     scr(ibs,imo) = fac * occ(imo) * cmoL(ibs,imo) * cmoR(ibs,imo)
!   end do
! end do
!
! write(*,"(/,'npsph, nmo',2i6)") npsph, nmo
! do ibs = 1, npsph
!   write(*,"(10f12.6)") scr(ibs,:)
! end do
!
! return
!end subroutine cdpopx

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Fix EFG values for open-shell atom & linear molecule.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine efgfix(iout,natom,xyz,eta,efgeig,efgsq, cxyz,eigscr)
 use Constant, only : Zero, Half, One
 implicit real*8(A-H,O-Z)
 dimension xyz(3,natom), efgeig(3), efgsq(3,3), cxyz(3,natom), eigscr(3,7)
 logical :: fixed, linear

 fixed = .false.
 linear= .false.

 if(natom == 1) then
   efgcc = Zero
   iaxis = 3
   if(abs(efgcc-efgeig(3)) > 1.0d-5 .or. eta > 1.0d-5) fixed = .true.
!else if(natom == 2) then
!  d = distance(xyz(1,1),xyz(1,2))
!  do i = 1, 3
!    x = abs(xyz(i,1) - xyz(i,2))
!    x = abs(x - d)
!    if(x < 1.0d-6) then    ! molecular axis
!      iaxis = i
!      efgcc = efgsq(iaxis,iaxis)
!      if(abs(efgcc-efgeig(3)) > 1.0d-5) fixed = .true.
!      exit
!    end if
!  end do
 else
   ! calculate mass-unweighted principal axes and moments of inertia
   ! eigenvalues --> eigscr(:,1), eigenvectors --> eigscr(:,2:4)
   call RotCons(iout,natom,xyz,cxyz,eigscr)
   ! write(iout,"(/, ' Mass-unweighted principal axes and moments of inertia in atomic units:',
   !   /,37x,'1',19x,'2',19x,'3',/, &
   !   5x,'Eigenvalues --  ',4x,3f20.6)") eigscr(1:3,1)
   ! write(iout,"(11x,'X',13x,3f20.6)")(eigscr(1,i),i=2,4)
   ! write(iout,"(11x,'Y',13x,3f20.6)")(eigscr(2,i),i=2,4)
   ! write(iout,"(11x,'Z',13x,3f20.6)")(eigscr(3,i),i=2,4)

   ! for linear molecule, an eigenvalue is zero
   do i = 1, 3
     if(abs(eigscr(i,1)) < 1.0d-5) then
       iaxis = i
       linear= .true.
       exit
     end if
   end do
   if(linear) then
     ! the molecular axis vector is eigscr(:,iaxis+1)
     eigscr(:,1) = eigscr(:,iaxis+1)
     ! Check whether it is in the x-, y, or z-direction
     do i = 1, 3
       if(abs(abs(eigscr(i,1)) - One) < 1.0d-5) then    ! molecular axis
         iaxis = i
         efgcc = efgsq(iaxis,iaxis)
         if(abs(efgcc-efgeig(3)) > 1.0d-5 .or. eta > 1.0d-5) fixed = .true.
         exit
       end if
     end do
   end if
 end if

 if(fixed) then
   efgeig(1) = -efgcc * Half
   efgeig(2) = -efgcc * Half
   efgeig(3) =  efgcc
   efgsq = Zero
   efgsq(1,1) = efgeig(1)
   efgsq(2,2) = efgeig(1)
   efgsq(3,3) = efgeig(1)
   efgsq(iaxis,iaxis) = efgeig(3)
   write(iout,"(/,39x,'Corrected EFG results due to broken symmetry of wavefunction:')")
   write(iout,"(40x,'Vaa =',f13.5,7x,3f14.5)") efgeig(1), efgsq(:,1)
   write(iout,"(40x,'Vbb =',f13.5,7x,3f14.5)") efgeig(2), efgsq(:,2)
   write(iout,"(40x,'Vcc =',f13.5,7x,3f14.5)") efgeig(3), efgsq(:,3)
   eta = 0.0d0
   if(abs(efgeig(3)) > 1.0d-8) eta = abs( (efgeig(1) - efgeig(2)) / efgeig(3) )
   write(iout,"(/,40x,'eta =',f13.5)") eta
 end if

end subroutine efgfix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Nuclear contribution to EFG tensor for the iAtm-th atom. See Eq. (5b) in J. Chem. Phys. 137, 054113 (2012).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine EFGNuc(iAtm,efglt, natom,za,xyz, npchg,pchg,pxyz, d1,d2,d5)
 use Constant, only : Three
 implicit real*8(A-H,O-Z)
 dimension efglt(6), za(*), xyz(3,*), pchg(*), pxyz(3,*), d1(3,*), d2(*), d5(*)

 do jAtm = 1, natom+npchg
   if(jAtm <= natom) then
     d1(1:3,jAtm) = xyz(1:3,iAtm) - xyz(1:3,jAtm)
   else
     ! point charges
     d1(1:3,jAtm) = xyz(1:3,iAtm) - pxyz(1:3,jAtm-natom)
   endif
   d2(jAtm) = d1(1,jAtm)**2 + d1(2,jAtm)**2 + d1(3,jAtm)**2
   d5(jAtm) = sqrt(d2(jAtm)) * d2(jAtm) * d2(jAtm)
 enddo

 do jAtm = 1, natom+npchg
   if(iAtm /= jAtm) then
     if(jAtm <= natom) then
       cha = za(jAtm)
     else
       cha = pchg(jAtm-natom)
     end if
     ! (3X^2-R^2)/R^5
     efglt(1) = efglt(1) + cha * (Three * d1(1,jAtm)**2 - d2(jAtm)) / d5(jAtm)
     ! 3XY/R^5
     efglt(2) = efglt(2) + cha *  Three * d1(1,jAtm) * d1(2,jAtm)   / d5(jAtm)
     ! (3Y^2-R^2)/R^5
     efglt(3) = efglt(3) + cha * (Three * d1(2,jAtm)**2 - d2(jAtm)) / d5(jAtm)
     ! 3XZ/R^5
     efglt(4) = efglt(4) + cha *  Three * d1(1,jAtm) * d1(3,jAtm)   / d5(jAtm)
     ! 3YZ/R^5
     efglt(5) = efglt(5) + cha *  Three * d1(2,jAtm) * d1(3,jAtm)   / d5(jAtm)
     ! (3Z^2-R^2)/R^5
     efglt(6) = efglt(6) + cha * (Three * d1(3,jAtm)**2 - d2(jAtm)) / d5(jAtm)
   endif
 enddo

 return
end

!--- END
