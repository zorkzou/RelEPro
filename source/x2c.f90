!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! X2C driver
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x2c_drv(iout,irelc,iham,natom,nuqatm,mapatm,iza,mxla,nshlla,iptexp,gtoexp,ialbs,  npsph,smat,tmat,vmat,wmat,xmat,  &
 rmat,scr)
 implicit real(kind=8) (a-h,o-z)
 dimension mapatm(*),iza(*),mxla(*),nshlla(*),iptexp(*),gtoexp(*),ialbs(*)
 dimension smat(*),tmat(*),vmat(*),wmat(*),xmat(*),rmat(*),scr(*)

 select case(iham)
   case(1)
     if(irelc == 1) then
       write(iout,"('   This is a non-relativistic calculation.')")
     else
       write(iout,"('   This is a non-relativistic two-component calculation.')")
     end if
   case(2)
     if(irelc == 1) then
       write(iout,"('   This is a sf-X2C calculation.')")
       call x1c(iout,npsph,smat,tmat,vmat,wmat,xmat,rmat,scr)
     else
       write(iout,"('   This is a X2C calculation.')")
       call x2c(iout,npsph,smat,tmat,vmat,wmat,xmat,rmat,scr)
       !write(*,"(/,'X matrix')")
       !call writemat(npsph+npsph,npsph+npsph,2,xmat,10)
       !write(*,"(/,'R matrix')")
       !call writemat(npsph+npsph,npsph+npsph,2,rmat,10)
     end if
   case(3,4)
     if(irelc == 1) then
       if(iham == 3) then
         write(iout,"('   This is a sf-X2C/DLXR calculation.')")
       else
         write(iout,"('   This is a sf-X2C/AXR calculation.')")
       end if
       ! unique atoms
       call uniatm(natom,nuqatm,mapatm,iza,mxla,nshlla,iptexp,gtoexp)
       call x1c_dlxr(iout, natom,nuqatm,mapatm,mxla,nshlla,ialbs,  npsph,smat,tmat,vmat,wmat,xmat,rmat,scr)
     else
       if(iham == 3) then
         write(iout,"('   This is a X2C/DLXR calculation.')")
       else
         write(iout,"('   This is a X2C/AXR calculation.')")
       end if
       ! unique atoms
       call uniatm(natom,nuqatm,mapatm,iza,mxla,nshlla,iptexp,gtoexp)
       call x2c_dlxr(iout, natom,nuqatm,mapatm,mxla,nshlla,ialbs,  npsph,smat,tmat,vmat,wmat,xmat,rmat,scr)
       !write(*,"(/,'X matrix')")
       !call writemat(npsph+npsph,npsph+npsph,2,xmat,10)
       !write(*,"(/,'R matrix')")
       !call writemat(npsph+npsph,npsph+npsph,2,rmat,10)
     end if
   case(5,6)
     if(irelc == 1) then
       if(iham == 5) then
         write(iout,"('   This is a sf-X2C/DLU calculation.')")
       else
         write(iout,"('   This is a sf-X2C/AU calculation.')")
       end if
       ! unique atoms
       call uniatm(natom,nuqatm,mapatm,iza,mxla,nshlla,iptexp,gtoexp)
       call x1c_dlu(iout, natom,nuqatm,mapatm,mxla,nshlla,ialbs,  npsph,smat,tmat,vmat,wmat,xmat,rmat,scr)
     else
       if(iham == 5) then
         write(iout,"('   This is a X2C/DLU calculation.')")
       else
         write(iout,"('   This is a X2C/AU calculation.')")
       end if
       ! unique atoms
       call uniatm(natom,nuqatm,mapatm,iza,mxla,nshlla,iptexp,gtoexp)
       call x2c_dlu(iout, natom,nuqatm,mapatm,mxla,nshlla,ialbs,  npsph,smat,tmat,vmat,wmat,xmat,rmat,scr)
       !write(*,"(/,'X matrix')")
       !call writemat(npsph+npsph,npsph+npsph,2,xmat,npsph+npsph)
       !write(*,"(/,'R matrix')")
       !call writemat(npsph+npsph,npsph+npsph,2,rmat,5)
     end if
   case default
     write(iout,"('   This is an unknown Hamiltonian.')")
     call estop(1)
 end select

 return
end subroutine x2c_drv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Driver of X1C.
!
! Input:
! np                : number of primitive spherical functions
! smat,tmat,vmat,wmat
!                   : S, T, V, W (pVp scaled by 1/(4cc))
!
! Output:
! xmat,rmat         : X and R
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x1c(iout,np,smat,tmat,vmat,wmat,xmat,rmat,scr)
 implicit real(kind=8) (a-h,o-z)
 dimension     :: smat(np*np), tmat(np*np), vmat(np*np), wmat(np*np), xmat(np*np), rmat(np*np), scr(np*np,4)
 allocatable   :: dirac(:,:), tmp(:,:)

 ! X matrix
 allocate(dirac(np*np,4), tmp(np*np,3))
 ! smat will be destroyed in x1c_atml, so copy smat to tmp(:,3)
 tmp(:,3) = smat
 call x1c_atml(iout,0,0,np,tmp(1,3),tmat,vmat,wmat,xmat,dirac,tmp,scr)

 ! R matrix
 call rcalc(np,smat,tmat,xmat,rmat,scr,tmp)
 deallocate(dirac, tmp)

 return
end subroutine x1c

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Driver of X1C/DLXR and X1C/AXR.
!
! Input:
! natom             : number of atoms
! nuqatm            : number of unique atoms
! mapatm            : map relationship of unique atoms. see sub. uniatm
! mxla              : max_L for each atom
! nshlla            : numbers of primitive shells for each atom
! ialbs             : starting position of each L in primitive sperical basis functions
! np                : number of primitive spherical functions
! smat,tmat,vmat,wmat
!                   : S, T, V, W (pVp scaled by 1/(4cc))
!
! Output:
! xmat,rmat         : X(DLX or AX) and R
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x1c_dlxr(iout, natom,nuqatm,mapatm,mxla,nshlla,ialbs, np,smat,tmat,vmat,wmat,xmat,rmat,scr)
 use Constant, only : maxlq, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension     :: mapatm(natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom), smat(np*np), tmat(np*np),  &
                  vmat(np*np), wmat(np*np), xmat(np*np), rmat(np*np), scr(np*np,4)
 allocatable   :: listeq(:), sla(:), tla(:), vla(:), wla(:), dirac(:,:), tmp(:,:)


 ! X matrix
 xmat = Zero

 npla = 10
 do iatm = 1, natom
   do lq = 0, mxla(iatm)
     npla = max(npla, (lq+lq+1)*nshlla(lq,iatm))
   end do
 end do
 nssla = npla*npla

 allocate(listeq(natom), sla(nssla), tla(nssla), vla(nssla), wla(nssla), dirac(nssla,4), tmp(nssla,2))

 do iua = 1, nuqatm
   call GetEqAtm(natom,iua,mapatm,neqa,listeq)
   iatm = listeq(1)
   do lq = 0, mxla(iatm)
     I0 = ialbs(lq,iatm)
     Iw = (lq+lq+1)*nshlla(lq,iatm)
     if(Iw < 1) cycle

     call TakDBlk(np,I0,Iw,smat,sla)
     call TakDBlk(np,I0,Iw,tmat,tla)
     call TakDBlk(np,I0,Iw,vmat,vla)
     call TakDBlk(np,I0,Iw,wmat,wla)

     ! atomic X will be saved at the very beginning of rmat
     call x1c_atml(iout,iatm,lq,Iw,sla,tla,vla,wla,rmat,dirac,tmp,scr)

     ! construct molecular X
     do jatm=1,neqa
       J0 = ialbs(lq,listeq(jatm))
       Jw = (lq+lq+1)*nshlla(lq,listeq(jatm))
       if(Iw /= Jw) then
         write(iout,"(/,' Error in sub. x1c_dlxr: Iw /= Jw')")
         call estop(1)
       end if

       call PutDBlk(np,J0,Jw,xmat,rmat)
     end do

   end do
 end do

 deallocate(listeq, sla, tla, vla, wla, dirac, tmp)

 ! R matrix
 allocate(tmp(np*np,3))
 call rcalc(np,smat,tmat,xmat,rmat,scr,tmp)
 deallocate(tmp)

 return
end subroutine x1c_dlxr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Driver of X1C/DLU and X1C/AU.
!
! Input:
! natom             : number of atoms
! nuqatm            : number of unique atoms
! mapatm            : map relationship of unique atoms. see sub. uniatm
! mxla              : max_L for each atom
! nshlla            : numbers of primitive shells for each atom
! ialbs             : starting position of each L in primitive sperical basis functions
! np                : number of primitive spherical functions
! smat,tmat,vmat,wmat
!                   : S, T, V, W (pVp scaled by 1/(4cc))
!
! Output:
! xmat,rmat         : X(DLX or AX) and R(DLR or AR)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x1c_dlu(iout, natom,nuqatm,mapatm,mxla,nshlla,ialbs, np,smat,tmat,vmat,wmat,xmat,rmat,scr)
 use Constant, only : maxlq, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension     :: mapatm(natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom), smat(np*np), tmat(np*np),  &
                  vmat(np*np), wmat(np*np), xmat(np*np), rmat(np*np), scr(np*np,4)
 allocatable   :: listeq(:), sla(:), tla(:), vla(:), wla(:), dirac(:,:), tmp(:,:)

 ! X matrix
 xmat = Zero

 npla = 10
 do iatm = 1, natom
   do lq = 0, mxla(iatm)
     npla = max(npla, (lq+lq+1)*nshlla(lq,iatm))
   end do
 end do
 nssla = npla*npla

 allocate(listeq(natom), sla(nssla), tla(nssla), vla(nssla), wla(nssla), dirac(nssla,4), tmp(nssla,2))

 do iua = 1, nuqatm
   call GetEqAtm(natom,iua,mapatm,neqa,listeq)
   iatm = listeq(1)
   do lq = 0, mxla(iatm)
     I0 = ialbs(lq,iatm)
     Iw = (lq+lq+1)*nshlla(lq,iatm)
     if(Iw < 1) cycle

     call TakDBlk(np,I0,Iw,smat,sla)
     call TakDBlk(np,I0,Iw,tmat,tla)
     call TakDBlk(np,I0,Iw,vmat,vla)
     call TakDBlk(np,I0,Iw,wmat,wla)

     ! atomic X will be saved at the very beginning of rmat
     call x1c_atml(iout,iatm,lq,Iw,sla,tla,vla,wla,rmat,dirac,tmp,scr)

     ! construct molecular X
     do jatm=1,neqa
       J0 = ialbs(lq,listeq(jatm))
       Jw = (lq+lq+1)*nshlla(lq,listeq(jatm))
       if(Iw /= Jw) then
         write(iout,"(/,' Error in sub.x1c_dlu: Iw /= Jw')")
         call estop(1)
       end if

       call PutDBlk(np,J0,Jw,xmat,rmat)
     end do

   end do
 end do

 ! R matrix
 rmat = Zero

 do iua = 1, nuqatm
   call GetEqAtm(natom,iua,mapatm,neqa,listeq)
   iatm = listeq(1)
   do lq = 0, mxla(iatm)
     I0 = ialbs(lq,iatm)
     Iw = (lq+lq+1)*nshlla(lq,iatm)
     if(Iw < 1) cycle

     call TakDBlk(np,I0,Iw,smat,sla)
     call TakDBlk(np,I0,Iw,tmat,tla)
     call TakDBlk(np,I0,Iw,xmat,vla)

     ! atomic R will be saved in wla
     call rcalc(Iw,sla,tla,vla,wla,scr,dirac)

     ! construct molecular R
     do jatm=1,neqa
       J0 = ialbs(lq,listeq(jatm))
       Jw = (lq+lq+1)*nshlla(lq,listeq(jatm))
       if(Iw /= Jw) then
         write(iout,"(/,' Error in sub.x1c_dlu: Iw /= Jw')")
         call estop(1)
       end if

       call PutDBlk(np,J0,Jw,rmat,wla)
     end do

   end do
 end do

 deallocate(listeq, sla, tla, vla, wla, dirac, tmp)

 return
end subroutine x1c_dlu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! R calculation (spin-free)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rcalc(np,smat,tmat,xmat,rmat,scr,tmp)
 use Constant, only : f_2cc, Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension     :: smat(np*np), tmat(np*np), xmat(np*np), rmat(np*np), scr(np*np,4), tmp(np*np,3), dum(1)

 call dgemm('n','n',np,np,np,One,tmat,np,xmat,np,Zero,scr(1,1),np)
 call dgemm('t','n',np,np,np,One,xmat,np,scr(1,1),np,Zero,rmat,np)
 rmat = smat + (One/f_2cc) * rmat

 ! S^(1/2) --> tmp(1), S^(-1/2) --> tmp(2)
 call sqivmt(np,smat,0,dum,1,tmp(1,1),1,tmp(1,2),scr(1,1),scr(1,2),scr(1,3))
 call dgemm('n','n',np,np,np,One,tmp(1,2),np,rmat,np,Zero,scr(1,1),np)
 call dgemm('n','n',np,np,np,One,scr(1,1),np,tmp(1,2),np,Zero,rmat,np)

 ! K^(-1/2) --> tmp(3)
 call sqivmt(np,rmat,0,dum,0,dum,1,tmp(1,3),scr(1,1),scr(1,2),scr(1,3))

 call dgemm('n','n',np,np,np,One,tmp(1,2),np,tmp(1,3),np,Zero,scr(1,1),np)
 call dgemm('n','n',np,np,np,One,scr(1,1),np,tmp(1,1),np,Zero,rmat,np)

 return
end subroutine rcalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Atomic X1C for L = lqa of atom iatm (iatm > 0), or molecular X1C (iatm = 0). Only the X matrix will be calculated.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x1c_atml(iout,iatm,lqa,np,sm,tm,vm,wm,xm,dirac,tmp,scr)
 use Constant, only : s_2cc, Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: sm(np*np),tm(np*np),vm(np*np),wm(np*np),xm(np*np),dirac(np*np,4),tmp(np*np,2),scr(np*np,4),dum(1)

 ! W-T --> W
 wm = wm - tm

 !<<< Construct Dirac matrix Equation and diagonalize

 ! S^(-1/2) --> tmp(1)
 call sqivmt(np,sm,0,dum,0,dum,1,tmp(1,1),scr(1,1),scr(1,2),scr(1,3))

 ! T^(1/2) --> sm, T^(-1/2) --> tmp(2)
 call sqivmt(np,tm,0,dum,1,sm,1,tmp(1,2),scr(1,1),scr(1,2),scr(1,3))

 ! make an overlap-weighted Dirac hamiltonian
 call mkdiracol(np,vm,wm,tmp(1,1),sm,tmp(1,2),dirac,scr(1,1),scr(1,2))

 ! eigenvalues and eigenvectors are saved in xm and dirac, respectively
 LWork=max(np*np*4,np*6)
 call DSYEV('V','L',np*2,dirac,np*2,xm,scr,LWork,INFO)
 if(INFO /= 0) then
   write(iout,"(/,' Error in sub. x1c_atml: Diagnolization failed.')")
   call estop(1)
 end if

 ! check positron eigenvalues
 call eigcheck(iout,np,xm,INFO)
 if (INFO == 1) then
   if(iatm > 0) then
     write(iout,"(/,' IVC problem encountered in sub. x1c_atml: L = ',i2,' of atom ',i4)") lqa, iatm
   else
     write(iout,"(/,' IVC problem encountered in sub. x1c_atml.')")
   end if
   call estop(1)
 end if

 ! overlap-unweighted eigenvectors A and B in parts (1,2) and (2,2), respectively
 !call TakeAB(np,dirac,scr(1,1),scr(1,2),scr(1,3),scr(1,4))
 call TakeAB(np,dirac,scr(1,3),scr(1,4))
 call dgemm('n','n',np,np,np,One,tmp(1,1),np,scr(1,3),np,Zero,scr(1,1),np)
 call dgemm('n','n',np,np,np,One,tmp(1,2),np,scr(1,4),np,Zero,scr(1,2),np)
 scr(:,2) = scr(:,2) * s_2cc

 !>>>

 ! calculate X
 call dgemm('n','t',np,np,np,One,scr(1,1),np,scr(1,1),np,Zero,tmp(1,1),np)
 call sqivmt(np,tmp(1,1),1,tmp(1,2),0,dum,0,dum,scr(1,3),scr(1,4),dirac)
 call dgemm('n','t',np,np,np,One,scr(1,2),np,scr(1,1),np,Zero,tmp(1,1),np)
 call dgemm('n','n',np,np,np,One,tmp(1,1),np,tmp(1,2),np,Zero,xm,np)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Construct one-electron overlap-weighted Hamiltonian of spin-free Dirac equation
!
!  Input
!   fc  : sqrt(2)*c
!   V   : V
!   TW  : W-T
!   Sm  : S^(-1/2)
!   Tp  : T^(1/2)
!   Tm  : T^(-1/2)
!
!  Output
!   H   : overlap-weighted Dirac Hamiltonian
!
!  Scratch:
!   Sc1, Sc2
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MkDiracOL(N,V,TW,Sm,Tp,Tm,H,Sc1,Sc2)
 use Constant, only : s_2cc, f_2cc, Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: V(N,N),TW(N,N),Sm(N,N),Tp(N,N),Tm(N,N),H(N,2,N,2),Sc1(N,N),Sc2(N,N)

 ! part (1,1): Sm * V * Sm
 call DGEMM('N','N',N,N,N,One,Sm,N,V,N,Zero,Sc2,N)
 call DGEMM('N','N',N,N,N,One,Sc2,N,Sm,N,Zero,Sc1,N)
 do i=1,N
   H(:,1,i,1) = Sc1(:,i)
 end do

 ! part (2,1): s_2cc * Tp * Sm
 call DGEMM('N','N',N,N,N,One,Tp,N,Sm,N,Zero,Sc1,N)
 Sc1 = Sc1 * s_2cc
 do i=1,N
   H(:,2,i,1) = Sc1(:,i)
 end do

 ! part (1,2)
 Sc1 = transpose(Sc1)
 do i=1,N
   H(:,1,i,2) = Sc1(:,i)
 end do

 ! part (2,2): f_2cc * Tm * (W - T) * Tm
 call DGEMM('N','N',N,N,N,One,Tm,N,TW,N,Zero,Sc2,N)
 call DGEMM('N','N',N,N,N,One,Sc2,N,Tm,N,Zero,Sc1,N)
 do i=1,N
   H(:,2,i,2) = Sc1(:,i) * f_2cc
 end do

 Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! check positron eigenvalues: E(i) must be smaller than -2mcc.
! eigenvalues should have been reordered.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine eigcheck(iout,n,e,ierr)
 use constant, only : f_2cc
 implicit real(kind=8) (a-h,o-z)
 dimension :: e(n)

 ierr = 0
 if(f_2cc > 1.0d8) return
 if(e(n) <= -f_2cc) return

 ierr = 1
 do i=n,1,-1
   if(e(i) <= -f_2cc) exit
   write(iout,"(2x,'---> ',i6,f30.8)")i,e(i)
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! take eigenvectors of the spin-free Dirac equation
!   | A C | --> | A- A+ |
!   | B D |     | B- B+ |
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TakeAB(n,v,c,d) 
 implicit real(kind=8) (a-h,o-z)
 dimension v(n,2,n,2),c(n,n),d(n,n)
 ! dimension a(n,n),b(n,n)

 do i=1,n
   ! a(:,i)=v(:,1,i,1)
   ! b(:,i)=v(:,2,i,1)
   c(:,i)=v(:,1,i,2)
   d(:,i)=v(:,2,i,2)
 enddo

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! take complex eigenvectors of the 4c-Dirac equation
!   | A C | --> | A- A+ |
!   | B D |     | B- B+ |
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TakeABcplx(n,v,cr,ci,dr,di)
 implicit real(kind=8) (a-h,o-z)
 dimension v(2,n,2,n,2),cr(n,n),ci(n,n),dr(n,n),di(n,n)

 do i=1,n
   cr(:,i)=v(1,:,1,i,2)
   ci(:,i)=v(2,:,1,i,2)
   dr(:,i)=v(1,:,2,i,2)
   di(:,i)=v(2,:,2,i,2)
 enddo

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Construct one-electron overlap-weighted Hamiltonian of 4c-Dirac equation
!
!  Input
!   fc : sqrt(2)*c
!   T, V, Wo, Wx, Wy, Wz
!   Z  : S^(-1/2)
!   P  : T^(1/2)
!   Q  : T^(-1/2)
!
!  Output
!   H   : overlap-weighted Dirac Hamiltonian
!
!  Scratch:
!   Sc1, Sc2, QQ
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Mk4cDiracOL(N,T,V,Wo,Wx,Wy,Wz,Z,P,Q,H,Sc1,Sc2,QQ)
 use Constant, only : s_2cc, f_2cc, Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: T(N,N),V(N,N),Wo(N,N),Wx(N,N),Wy(N,N),Wz(N,N),Z(N,N),P(N,N),Q(N,N),H(2,4*N,4*N),Sc1(N,N*4),Sc2(N*2,N*2),QQ(N,N*4)

 H = Zero
 N2= N+N

 ! part (1,1): Z * V * Z
 call DGEMM('N','N',N,N,N,One,Z,N,V,N,Zero,Sc2,N)
 call DGEMM('N','N',N,N,N,One,Sc2,N,Z,N,Zero,Sc1,N)
 do i = 1, N
   i2 = i + i
   i1 = i2- 1
   do j = 1, N
     j2 = j + j
     j1 = j2- 1
     H(1,j1,i1) = Sc1(j,i)
     H(1,j2,i2) = Sc1(j,i)
   end do
 end do

 ! part (2,1): s_2cc * P * Z
 call DGEMM('N','N',N,N,N,One,P,N,Z,N,Zero,Sc1,N)
 call AScale(N*N,s_2cc,Sc1,Sc1)
 do i = 1, N
   i2 = i + i
   i1 = i2- 1
   do j = 1, N
     j2 = j + j + N2
     j1 = j2- 1
     H(1,j1,i1) = Sc1(j,i)
     H(1,j2,i2) = Sc1(j,i)
   end do
 end do

 ! part (1,2)
 call TranspSq(N,Sc1)
 do i = 1,N
   i2 = i + i + N2
   i1 = i2- 1
   do j = 1, N
     j2 = j + j
     j1 = j2- 1
     H(1,j1,i1) = Sc1(j,i)
     H(1,j2,i2) = Sc1(j,i)
   end do
 end do

 ! part (2,2): f_2cc * Q * (W - T) * Q
 call DblMat(N,N,Q,QQ)
 ! real
 call ASMD(2,N*N,Wo,T,Sc1)
 call WoaddWy(N,Sc1,Wy,Sc2)
 call DGEMM('N','N',N2,N2,N2,One,QQ,N2,Sc2,N2,Zero,Sc1,N2)
 call DGEMM('N','N',N2,N2,N2,One,Sc1,N2,QQ,N2,Zero,Sc2,N2)
 do i = 1, N2
   i1 = i + N2
   do j = 1, N2
     j1 = j + N2
     H(1,j1,i1)= Sc2(j,i) * f_2cc
   end do
 end do
 ! imag
 call WxaddWz(N,Wx,Wz,Sc2)
 call DGEMM('N','N',N2,N2,N2,One,QQ,N2,Sc2,N2,Zero,Sc1,N2)
 call DGEMM('N','N',N2,N2,N2,One,Sc1,N2,QQ,N2,Zero,Sc2,N2)
 do i = 1, N2
   i1 = i + N2
   do j = 1, N2
     j1 = j + N2
     H(2,j1,i1)= Sc2(j,i) * f_2cc
   end do
 end do

 Return
End subroutine Mk4cDiracOL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Construct the real part of W(o,y) matrix
!                 i1      i2
!           j1  | Wo      Wy|
!           j2  |-Wy      Wo|
!
! NOTE. Minus signs of Wx, Wy, and Wz in Eq. (26) of J. Chem. Phys. 142, 214106 (2015) are missing. See Eq. (14.17) in
! Reiher & Wolf, Relativistic Quantum Chemistry, Ed.2, 2015,
! or Eq. (3.9) in
! Zhang, Theor. Chem. Acc. 133, 1588 (2014).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WoaddWy(N,Wo,Wy,WR)
 implicit real(kind=8) (a-h,o-z)
 dimension :: Wo(N,N),Wy(N,N),WR(2,N,2,N)

 do i = 1, N
   do j = 1, N
     WR(1,j,1,i) = Wo(j,i)
     WR(2,j,1,i) =-Wy(j,i)
     WR(1,j,2,i) = Wy(j,i)
     WR(2,j,2,i) = Wo(j,i)
   end do
 end do

 Return
End subroutine WoaddWy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Construct the imaginary part of W matrix
!                 i1      i2
!           j1  | Wz      Wx|
!           j2  | Wx     -Wz|
!
! NOTE. Minus signs of Wx, Wy, and Wz in Eq. (26) of J. Chem. Phys. 142, 214106 (2015) are missing. See Eq. (14.17) in
! Reiher & Wolf, Relativistic Quantum Chemistry, Ed.2, 2015,
! or Eq. (3.9) in
! Zhang, Theor. Chem. Acc. 133, 1588 (2014).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WxaddWz(N,Wx,Wz,WI)
 implicit real(kind=8) (a-h,o-z)
 dimension :: Wx(N,N),Wz(N,N),WI(2,N,2,N)

 do i = 1, N
   do j = 1, N
     WI(1,j,1,i) = Wz(j,i)
     WI(2,j,1,i) = Wx(j,i)
     WI(1,j,2,i) = Wx(j,i)
     WI(2,j,2,i) =-Wz(j,i)
   end do
 end do

 Return
End subroutine WxaddWz

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Atomic X2C for L = lqa of atom iatm (iatm > 0), or molecular X2C (iatm = 0). Only the X matrix will be calculated.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x2c_atml(iout,iatm,lqa,np,sm,tm,vm,wo,wx,wy,wz,xm,dirac,tmp,scr)
 use Constant, only : s_2cc, Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: sm(np*np),tm(np*np),vm(np*np),wo(np*np),wx(np*np),wy(np*np),wz(np*np),xm(np*np*4,2),dirac(np*np,32),  &
              tmp(np*np,6),scr(np*np,32),dum(1)

 nd = np + np
 nq = nd + nd

 !<<< Construct Dirac matrix Equation and diagonalize

 ! let Z = S^(-1/2) --> tmp(1)
 call sqivmt(np,sm,0,dum,0,dum,1,tmp(1,1),scr(1,1),scr(1,2),scr(1,3))

 ! let P = T^(1/2) --> sm and Q = T^(-1/2) --> tmp(2)
 call sqivmt(np,tm,0,dum,1,sm,1,tmp(1,2),scr(1,1),scr(1,2),scr(1,3))

 ! make an overlap-weighted Dirac hamiltonian
 call mk4cdiracol(np,tm,vm,wo,wx,wy,wz,tmp(1,1),sm,tmp(1,2),dirac,tmp(1,3),scr(1,1),scr(1,5))

 ! eigenvalues and eigenvectors are saved in xm and dirac, respectively
 LWork=nq*nq
 call ZHEEV('V','L',nq,dirac,nq,xm,scr,LWork,tmp(1,3),INFO)
 if(INFO /= 0) then
   write(iout,"(/,' Error in sub. x2c_atml: Diagnolization failed.')")
   call estop(1)
 end if

 ! check positron eigenvalues
 call eigcheck(iout,nd,xm,INFO)
 if (INFO == 1) then
   if(iatm > 0) then
     write(iout,"(/,' IVC problem encountered in sub. x2c_atml: L = ',i2,' of atom ',i4)") lqa, iatm
   else
     write(iout,"(/,' IVC problem encountered in sub. x2c_atml.')")
   end if
   call estop(1)
 end if

 ! overlap-unweighted eigenvectors A and B in parts (1,2) and (2,2), respectively
 call TakeABcplx(nd,dirac,scr(1,17),scr(1,21),scr(1,25),scr(1,29))
 ! A+: real & imag
 call DblMat(np,np,tmp(1,1),tmp(1,3))
 call MatMult(1,nd,nd,nd,One,Zero,tmp(1,3),scr(1,17),scr(1, 1))
 call MatMult(1,nd,nd,nd,One,Zero,tmp(1,3),scr(1,21),scr(1, 5))
 ! B+: real & imag
 call AScale(np*np,s_2cc,tmp(1,2),sm)
 call DblMat(np,np,sm,tmp(1,3))
 call MatMult(1,nd,nd,nd,One,Zero,tmp(1,3),scr(1,25),scr(1, 9))
 call MatMult(1,nd,nd,nd,One,Zero,tmp(1,3),scr(1,29),scr(1,13))

 !>>>

 ! calculate X = B * C
 call ReImComb(nd,nd,scr(1,1),scr(1,17))               ! complex A+
 call InvMatC(nd,scr(1,17),scr(1,25),tmp(1,3),INFO)    ! complex C = inv(A+)
 If(INFO /= 0) then
   write(iout,"(/,' Inversion failed in sub. x2c_atml.')")
   call estop(1)
 end if
 call ReImSplit(nd,nd,scr(1,25),scr(1,1))              ! real & imag C
 ! Xr = Br * Cr - Bi * Ci
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,13),scr(1, 5),xm(1,1))
 call MatMult(1,nd,nd,nd,One,-One,scr(1, 9),scr(1, 1),xm(1,1))
 ! Xi = Bi * Cr + Br * Ci
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,13),scr(1, 1),xm(1,2))
 call MatMult(1,nd,nd,nd,One, One,scr(1, 9),scr(1, 5),xm(1,2))

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Driver of X2C.
!
! Input:
! np                : number of primitive spherical functions
! smat,tmat,vmat,wmat
!                   : S, T, V, W (pVp scaled by 1/(4cc))
!
! Output:
! xmat,rmat         : X and R
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x2c(iout,np,smat,tmat,vmat,wmat,xmat,rmat,scr)
 implicit real(kind=8) (a-h,o-z)
 dimension     :: smat(np*np), tmat(np*np), vmat(np*np), wmat(np*np,4), xmat(np*np*8), rmat(np*np,8), scr(np*np,32)
 allocatable   :: dirac(:,:)

 ! smat will be destroyed in x2c_atml, so copy smat to rmat(:,1)
 rmat(:,1) = smat

 ! X matrix
 allocate(dirac(np*np,32))
 call x2c_atml(iout,0,0,np,rmat(1,1),tmat,vmat,wmat(1,1),wmat(1,2),wmat(1,3),wmat(1,4),xmat,dirac,rmat(1,2),scr)
 deallocate(dirac)

 ! R matrix
 call rcalc2(np,smat,tmat,xmat,rmat,scr)

 return
end subroutine x2c

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! R calculation (two-component)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rcalc2(np,smat,tmat,xmat,rmat,scr)
 use Constant, only : f_2cc, Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension     :: smat(np*np), tmat(np*np), xmat(np*np*4,2), rmat(np*np*4,2), scr(np*np,32), dum(1)

 nd = np + np

 ! S~ = S + X^H * T * X / (2cc) --> rmat
 call AScale(np*np,(One/f_2cc),tmat,scr(1,1))
 call DblMat(np,np,scr(1,1),scr(1,2))
 ! real(S~) = S + (Xr^T * T * Xr + Xi^T * T * Xi) / (2cc)
 call DblMat(np,np,smat,rmat(1,1))
 call MatMult(2,nd,nd,nd,One,Zero,xmat(1,1),scr(1,2),scr(1,6))
 call MatMult(1,nd,nd,nd,One, One,scr(1,6),xmat(1,1),rmat(1,1))
 call MatMult(2,nd,nd,nd,One,Zero,xmat(1,2),scr(1,2),scr(1,10))
 call MatMult(1,nd,nd,nd,One, One,scr(1,10),xmat(1,2),rmat(1,1))
 ! imag(S~) = (Xr^T * T * Xi - Xi^T * T * Xr) / (2cc)
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,6),xmat(1,2),scr(1,10))
 call TrAddSub(1,nd,scr(1,10),rmat(1,2))

 ! S^(1/2) --> scr(1), S^(-1/2) --> scr(2)
 call sqivmt(np,smat,0,dum,1,scr(1,1),1,scr(1,2),scr(1,3),scr(1,4),scr(1,5))
 ! doubled S^(-1/2) --> scr(3)
 call DblMat(np,np,scr(1,2),scr(1,3))
 ! KK --> rmat
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,3),rmat(1,1),scr(1,7))
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,7),scr(1,3),rmat(1,1))
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,3),rmat(1,2),scr(1,7))
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,7),scr(1,3),rmat(1,2))

 ! K --> rmat
 call ReImComb(nd,nd,rmat,scr(1,7))
 call SqrtMCh(nd,1,scr(1,7),scr(1,15),rmat,scr(1,31),scr(1,23),INFO)
 If(INFO /= 0) then
   write(*,"(/,' SqrtMCh failed in sub. rcalc2.')")
   call estop(1)
 end if
 call ReImSplit(nd,nd,scr(1,15),rmat)

 ! R matrix
 ! doubled S^(1/2) --> scr(7)
 call DblMat(np,np,scr(1,1),scr(1,7))
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,3),rmat(1,1),scr(1,11))
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,11),scr(1,7),rmat(1,1))
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,3),rmat(1,2),scr(1,11))
 call MatMult(1,nd,nd,nd,One,Zero,scr(1,11),scr(1,7),rmat(1,2))

 return
end subroutine rcalc2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Driver of X2C/DLXR and X2C/AXR.
!
! Input:
! natom             : number of atoms
! nuqatm            : number of unique atoms
! mapatm            : map relationship of unique atoms. see sub. uniatm
! mxla              : max_L for each atom
! nshlla            : numbers of primitive shells for each atom
! ialbs             : starting position of each L in primitive sperical basis functions
! np                : number of primitive spherical functions
! smat,tmat,vmat,wmat
!                   : S, T, V, W (pVp scaled by 1/(4cc))
!
! Output:
! xmat,rmat         : X(DLX or AX) and R
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x2c_dlxr(iout, natom,nuqatm,mapatm,mxla,nshlla,ialbs, np,smat,tmat,vmat,wmat,xmat,rmat,scr)
 use Constant, only : maxlq, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension     :: mapatm(natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom), smat(np*np), tmat(np*np),  &
                  vmat(np*np), wmat(np*np,4), xmat(np*np*4,2), rmat(np*np*8), scr(np*np,32)
 allocatable   :: listeq(:), sla(:), tla(:), vla(:), wla(:,:), dirac(:,:), tmp(:,:)

 nd = np + np

 ! X matrix
 xmat = Zero

 npla = 10
 do iatm = 1, natom
   do lq = 0, mxla(iatm)
     npla = max(npla, (lq+lq+1)*nshlla(lq,iatm))
   end do
 end do
 nssla = npla*npla

 allocate(listeq(natom), sla(nssla), tla(nssla), vla(nssla), wla(nssla,4), dirac(nssla,32), tmp(nssla,6))

 do iua = 1, nuqatm
   call GetEqAtm(natom,iua,mapatm,neqa,listeq)
   iatm = listeq(1)
   do lq = 0, mxla(iatm)
     I0 = ialbs(lq,iatm)
     Iw = (lq+lq+1)*nshlla(lq,iatm)
     if(Iw < 1) cycle

     call TakDBlk(np,I0,Iw,smat,sla)
     call TakDBlk(np,I0,Iw,tmat,tla)
     call TakDBlk(np,I0,Iw,vmat,vla)
     call TakDBlk(np,I0,Iw,wmat(1,1),wla(1,1))
     call TakDBlk(np,I0,Iw,wmat(1,2),wla(1,2))
     call TakDBlk(np,I0,Iw,wmat(1,3),wla(1,3))
     call TakDBlk(np,I0,Iw,wmat(1,4),wla(1,4))

     ! atomic X will be saved at the very beginning of rmat
     call x2c_atml(iout,iatm,lq,Iw,sla,tla,vla,wla(1,1),wla(1,2),wla(1,3),wla(1,4),rmat,dirac,tmp,scr)

     ! construct molecular X
     do jatm=1,neqa
       J0 = ialbs(lq,listeq(jatm))
       J0 = J0+J0-1
       Jw = (lq+lq+1)*nshlla(lq,listeq(jatm))
       Jw = Jw+Jw
       if(Iw+Iw /= Jw) then
         write(iout,"(/,' Error in sub. x2c_dlxr: 2Iw /= Jw')")
         call estop(1)
       end if

       call PutDBlk(nd,J0,Jw,xmat(1,1),rmat)
       call PutDBlk(nd,J0,Jw,xmat(1,2),rmat(1+Jw*Jw))
     end do

   end do
 end do

 deallocate(listeq, sla, tla, vla, wla, dirac, tmp)

 ! R matrix
 call rcalc2(np,smat,tmat,xmat,rmat,scr)

 return
end subroutine x2c_dlxr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Driver of X2C/DLU and X2C/AU.
!
! Input:
! natom             : number of atoms
! nuqatm            : number of unique atoms
! mapatm            : map relationship of unique atoms. see sub. uniatm
! mxla              : max_L for each atom
! nshlla            : numbers of primitive shells for each atom
! ialbs             : starting position of each L in primitive sperical basis functions
! np                : number of primitive spherical functions
! smat,tmat,vmat,wmat
!                   : S, T, V, W (pVp scaled by 1/(4cc))
!
! Output:
! xmat,rmat         : X(DLX or AX) and R(DLR or AR)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine x2c_dlu(iout, natom,nuqatm,mapatm,mxla,nshlla,ialbs, np,smat,tmat,vmat,wmat,xmat,rmat,scr)
 use Constant, only : maxlq, Zero
 implicit real(kind=8) (a-h,o-z)
 dimension     :: mapatm(natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom), smat(np*np), tmat(np*np),  &
                  vmat(np*np), wmat(np*np,4), xmat(np*np*4,2), rmat(np*np*8), scr(np*np,32)
 allocatable   :: listeq(:), sla(:), tla(:), vla(:), wla(:,:), dirac(:), tmp(:)

 nd = np + np

 ! X matrix
 xmat = Zero

 npla = 10
 do iatm = 1, natom
   do lq = 0, mxla(iatm)
     npla = max(npla, (lq+lq+1)*nshlla(lq,iatm))
   end do
 end do
 nssla = npla*npla

 allocate(listeq(natom), sla(nssla), tla(nssla), vla(nssla), wla(nssla,4), dirac(nssla*32), tmp(nssla*8))

 do iua = 1, nuqatm
   call GetEqAtm(natom,iua,mapatm,neqa,listeq)
   iatm = listeq(1)
   do lq = 0, mxla(iatm)
     I0 = ialbs(lq,iatm)
     Iw = (lq+lq+1)*nshlla(lq,iatm)
     if(Iw < 1) cycle

     call TakDBlk(np,I0,Iw,smat,sla)
     call TakDBlk(np,I0,Iw,tmat,tla)
     call TakDBlk(np,I0,Iw,vmat,vla)
     call TakDBlk(np,I0,Iw,wmat(1,1),wla(1,1))
     call TakDBlk(np,I0,Iw,wmat(1,2),wla(1,2))
     call TakDBlk(np,I0,Iw,wmat(1,3),wla(1,3))
     call TakDBlk(np,I0,Iw,wmat(1,4),wla(1,4))

     ! atomic X will be saved at the very beginning of rmat
     call x2c_atml(iout,iatm,lq,Iw,sla,tla,vla,wla(1,1),wla(1,2),wla(1,3),wla(1,4),rmat,dirac,tmp,scr)

     ! construct molecular X
     do jatm=1,neqa
       J0 = ialbs(lq,listeq(jatm))
       J0 = J0+J0-1
       Jw = (lq+lq+1)*nshlla(lq,listeq(jatm))
       Jw = Jw+Jw
       if(Iw+Iw /= Jw) then
         write(iout,"(/,' Error in sub. x2c_dlu: 2Iw /= Jw')")
         call estop(1)
       end if

       call PutDBlk(nd,J0,Jw,xmat(1,1),rmat)
       call PutDBlk(nd,J0,Jw,xmat(1,2),rmat(1+Jw*Jw))
     end do

   end do
 end do

 ! R matrix
 rmat = Zero

 do iua = 1, nuqatm
   call GetEqAtm(natom,iua,mapatm,neqa,listeq)
   iatm = listeq(1)
   do lq = 0, mxla(iatm)
     I0 = ialbs(lq,iatm)
     Iw = (lq+lq+1)*nshlla(lq,iatm)
     if(Iw < 1) cycle

     call TakDBlk(np,I0,Iw,smat,sla)
     call TakDBlk(np,I0,Iw,tmat,tla)
     call TakDBlk(nd,(I0+I0-1),(Iw+Iw),xmat(1,1),tmp)
     call TakDBlk(nd,(I0+I0-1),(Iw+Iw),xmat(1,2),tmp(1+4*Iw*Iw))

     ! atomic R will be saved in dirac
     call rcalc2(Iw,sla,tla,tmp,dirac,scr)

     ! construct molecular R
     do jatm=1,neqa
       J0 = ialbs(lq,listeq(jatm))
       J0 = J0+J0-1
       Jw = (lq+lq+1)*nshlla(lq,listeq(jatm))
       Jw = Jw+Jw
       if(Iw+Iw /= Jw) then
         write(iout,"(/,' Error in sub. x2c_dlu: 2Iw /= Jw')")
         call estop(1)
       end if

       call PutDBlk(nd,J0,Jw,rmat,dirac)
       call PutDBlk(nd,J0,Jw,rmat(1+nd*nd),dirac(1+Jw*Jw))
     end do

   end do
 end do

 deallocate(listeq, sla, tla, vla, wla, dirac, tmp)

 return
end subroutine x2c_dlu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Density transfermation of X2C (Pv & Pw only)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine X2CDenTran(irelc,np,xmat,rmat,pdens,pdenv,pdenw,scr)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: xmat(np*np,*), rmat(np*np,*), pdens(np*np,*), pdenv(np*np,*), pdenw(np*np,*), scr(np*np,*)

 nd = np + np

 ! Pv = R P R'; the <Pxv> term has been ignored
 if(irelc == 1) then
   call MatMult(1,np,np,np,One,Zero,rmat,pdens,scr)
   call MatMult(3,np,np,np,One,Zero,scr,rmat,pdenv)
 else if(irelc == 2) then
   ! Re(Pv) = Rr * Pr * Rr^T + Ri * Pr * Ri^T + TrAdd(Rr * Pi * Ri^T)
   call MatMult(1,nd,nd,nd,One,Zero,rmat(1,1),pdens(1,5),scr(1,1))
   call MatMult(3,nd,nd,nd,One,Zero,scr(1,1),rmat(1,5),scr(1,5))
   call TrAddSub(0,nd,scr(1,5),pdenv(1,1))
   call MatMult(1,nd,nd,nd,One,Zero,rmat(1,1),pdens(1,1),scr(1,1))
   call MatMult(3,nd,nd,nd,One,One ,scr(1,1),rmat(1,1),pdenv(1,1))
   call MatMult(1,nd,nd,nd,One,Zero,rmat(1,5),pdens(1,1),scr(1,1))
   call MatMult(3,nd,nd,nd,One,One ,scr(1,1),rmat(1,5),pdenv(1,1))
   ! Im(Pv) = Rr * Pi * Rr^T + Ri * Pi * Ri^T + TrSub(Ri * Pr * Rr^T)
   call MatMult(1,nd,nd,nd,One,Zero,rmat(1,5),pdens(1,1),scr(1,1))
   call MatMult(3,nd,nd,nd,One,Zero,scr(1,1),rmat(1,1),scr(1,5))
   call TrAddSub(1,nd,scr(1,5),pdenv(1,5))
   call MatMult(1,nd,nd,nd,One,Zero,rmat(1,1),pdens(1,5),scr(1,1))
   call MatMult(3,nd,nd,nd,One,One ,scr(1,1),rmat(1,1),pdenv(1,5))
   call MatMult(1,nd,nd,nd,One,Zero,rmat(1,5),pdens(1,5),scr(1,1))
   call MatMult(3,nd,nd,nd,One,One ,scr(1,1),rmat(1,5),pdenv(1,5))
 end if

 ! Pw = X Pv X'; the <Pxw> term has been ignored
 if(irelc == 1) then
   call MatMult(1,np,np,np,One,Zero,xmat,pdenv,scr)
   call MatMult(3,np,np,np,One,Zero,scr,xmat,pdenw)
 else if(irelc == 2) then
   ! Re(Pw)
   call MatMult(1,nd,nd,nd,One,Zero,xmat(1,1),pdenv(1,5),scr(1,1))
   call MatMult(3,nd,nd,nd,One,Zero,scr(1,1),xmat(1,5),scr(1,5))
   call TrAddSub(0,nd,scr(1,5),pdenw(1,1))
   call MatMult(1,nd,nd,nd,One,Zero,xmat(1,1),pdenv(1,1),scr(1,1))
   call MatMult(3,nd,nd,nd,One,One ,scr(1,1),xmat(1,1),pdenw(1,1))
   call MatMult(1,nd,nd,nd,One,Zero,xmat(1,5),pdenv(1,1),scr(1,1))
   call MatMult(3,nd,nd,nd,One,One ,scr(1,1),xmat(1,5),pdenw(1,1))
   ! Im(Pw)
   call MatMult(1,nd,nd,nd,One,Zero,xmat(1,5),pdenv(1,1),scr(1,1))
   call MatMult(3,nd,nd,nd,One,Zero,scr(1,1),xmat(1,1),scr(1,5))
   call TrAddSub(1,nd,scr(1,5),pdenw(1,5))
   call MatMult(1,nd,nd,nd,One,Zero,xmat(1,1),pdenv(1,5),scr(1,1))
   call MatMult(3,nd,nd,nd,One,One ,scr(1,1),xmat(1,1),pdenw(1,5))
   call MatMult(1,nd,nd,nd,One,Zero,xmat(1,5),pdenv(1,5),scr(1,1))
   call MatMult(3,nd,nd,nd,One,One ,scr(1,1),xmat(1,5),pdenw(1,5))
 end if

 return
end subroutine X2CDenTran

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MO transformation for X2C : CR = R * C, CU = U * C = X * CR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine X2CMOTran(irelc,np,nmo,xmat,rmat,cm,cr,cu)
  use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension :: xmat(np*np,*), rmat(np*np,*), cm(irelc*np*nmo,irelc), cr(irelc*np*nmo,irelc), cu(irelc*np*nmo,irelc)

 nd = np + np

 if(irelc == 1) then
   call MatMult(1,np,np,nmo,One,Zero,rmat,cm,cr)
   call MatMult(1,np,np,nmo,One,Zero,xmat,cr,cu)
 else if(irelc == 2) then
   ! Re(CR) = Rr * Cr - Ri * Ci
   call MatMult(1,nd,nd,nmo,One,Zero,rmat(1,5),cm(1,2),cr(1,1))
   call MatMult(1,nd,nd,nmo,One,-One,rmat(1,1),cm(1,1),cr(1,1))
   ! Im(CR) = Ri * Cr + Rr * Ci
   call MatMult(1,nd,nd,nmo,One,Zero,rmat(1,5),cm(1,1),cr(1,2))
   call MatMult(1,nd,nd,nmo,One, One,rmat(1,1),cm(1,2),cr(1,2))
   ! Re(CU) = Xr * Cr - Xi * Ci
   call MatMult(1,nd,nd,nmo,One,Zero,xmat(1,5),cr(1,2),cu(1,1))
   call MatMult(1,nd,nd,nmo,One,-One,xmat(1,1),cr(1,1),cu(1,1))
   ! Im(CU) = Xi * Cr + Xr * Ci
   call MatMult(1,nd,nd,nmo,One,Zero,xmat(1,5),cr(1,1),cu(1,2))
   call MatMult(1,nd,nd,nmo,One, One,xmat(1,1),cr(1,2),cu(1,2))
 end if

 return
end subroutine X2CMOTran

!--- END

