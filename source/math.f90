!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Density matrix calculation: P = C O C^T
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine calc_den(ncbas,nmo,cmo,occ,pden)
 use Constant, only : Zero, occtol
 implicit real(kind=8) (a-h,o-z)
 dimension         :: cmo(ncbas,nmo), occ(nmo) ,pden(ncbas,ncbas)

 pden = Zero

 do i=1,nmo
   if(abs(occ(i)) < occtol) cycle
   do ibs=1,ncbas
     do jbs=1,ibs
       pden(jbs,ibs) = pden(jbs,ibs) + occ(i)*cmo(ibs,i)*cmo(jbs,i)
     end do
   end do
 end do

 do ibs=1,ncbas
   do jbs=1,ibs-1
     pden(ibs,jbs) = pden(jbs,ibs)
   end do
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! B(*) = c * A(*)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine AScale(N,c,A,B)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(*),B(*)

 Do I = 1,N
   B(I) = c*A(I)
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Mode = 1: C = beta * C + alpha * A * B
!        2: C = beta * C + alpha * A^T * B
!        3: C = beta * C + alpha * A * B^T
!        4: C = beta * C + alpha * A^T * B^T
! where C(MxN), op(A)(MxL), and op(B)(LxN)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine MatMult(Mode,M,L,N,alpha,beta,A,B,C)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: A(*),B(*),C(*)

 LDC=M
 if(Mode == 1)then
   LDA=M
   LDB=L
   call dgemm('N','N',M,N,L,alpha,A,LDA,B,LDB,beta,C,LDC)
 else if(Mode == 2)then
   LDA=L
   LDB=L
   call dgemm('T','N',M,N,L,alpha,A,LDA,B,LDB,beta,C,LDC)
 else if(Mode == 3)then
   LDA=M
   LDB=N
   call dgemm('N','T',M,N,L,alpha,A,LDA,B,LDB,beta,C,LDC)
 else if(Mode == 4)then
   LDA=L
   LDB=N
   call dgemm('T','T',M,N,L,alpha,A,LDA,B,LDB,beta,C,LDC)
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Trace of the product of two SYMMETRIC matrices A and B.
! If one is SYMMETRIC and the other is ANTISYMMETRIC, the product is zero.
! The product should be multiplied by -1 if A and B are ANTISYMMETRIC.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double precision function traceSS(n,A,B)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension A(n,n), B(n,n)

 su = Zero
 do i = 1, n
   su = su + dotx(n, A(1,i), B(1,i))
 end do

 traceSS = su

 return
end function traceSS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! (mode = 0) B = A + A^T, where B is symmetric
! (mode = 1) B = A - A^T, where B is antisymmetric
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TrAddSub(mode,N,A,B)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension         :: A(N,N), B(N,N)

 if(mode == 0) then

   do i=1,N
     do j=1,i-1
       B(j,i)=A(j,i)+A(i,j)
       B(i,j)=B(j,i)
     end do
     B(i,i)=A(i,i)+A(i,i)
   end do

 else if(mode == 1) then

   do i=1,N
     do j=1,i-1
       B(j,i)=A(j,i)-A(i,j)
       B(i,j)=-B(j,i)
     end do
     ! B(i,i)=A(i,i)-A(i,i)
     B(i,i)=Zero
   end do

 end if

 Return
End subroutine TrAddSub

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! (mode = 1) C = A + B
! (mode = 2) C = A - B
! (mode = 3) C[i] = A[i] * B[i]
! (mode = 4) C[i] = A[i] / B[i]
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ASMD(mode,N,A,B,C)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension         :: A(N), B(N), C(N)

 if(mode == 1) then

   C = A + B

 else if(mode == 2) then

   C = A - B

 else if(mode == 3) then

   do i=1,N
     C(i) = A(i) * B(i)
   end do

 else if(mode == 4) then

   do i=1,N
     C(i) = A(i) / B(i)
   end do

 end if

 Return
End subroutine ASMD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! dot product of A and B
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double precision function dotx(n,A,B)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: A(n), B(n)

 dotx = dot_product(A,B)

 return
end function dotx

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Distance between points A and B
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double precision function distance(A,B)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: A(3), B(3)

 distance = sqrt((A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2)

 return
end function distance

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! sum(abs())
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function abssum(n,a)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension         :: a(n)

 x = Zero
 do i = 1, n
   x = x + abs(a(i))
 end do
 abssum = x

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! isum()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isum(n,ia)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: ia(n)

 ix = 0
 do i = 1, n
   ix = ix + ia(i)
 end do
 isum = ix

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! c = a . b
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine scamul(n,a,b,c)
 Implicit Real*8(A-H,O-Z)
 Dimension a(n), b(n), c(n)

 do i = 1, n
   c(i) = b(i) * a(i)
 end do

 return
end subroutine scamul

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! A(N,N)^T
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine TranspSq(N,A)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: A(N,N)

 A = transpose(A)

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! ATr(N,M) = A(M,N)^T
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine Transp(M,N,A,ATr)
 implicit real(kind=8) (a-h,o-z)
 dimension         :: A(M,N),ATr(N,M)

 do i=1,N
   do j=1,M
     ATr(i,j) = A(j,i)
   end do
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! (I)!, I=0~N
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FacCal(N,facdat)
 use Constant, only : One
 implicit real(kind=8) (a-h,o-z)
 dimension :: facdat(N+1)

 facdat(1) = One
 Do I = 1, N
   facdat(I+1) = facdat(I) * dble(I)
 end do

 Return
End Subroutine FacCal

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tlm1ff = (2 (N-1) -1)!!
! So if you calculate (2*J-1)!!, you have to use Tlm1ff(J+1)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tlm1ff(N)
 use Constant, only : One
 implicit real(kind=8) (a-h,o-z)

 Tlm1ff = One
 Do m = 2,N-1
   Tlm1ff = Tlm1ff * dble(m+m-1)
 end do

 return
end function Tlm1ff

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the combination factors of i!/[j!(i-j)!], where i=0,...,N; j=0,...,i
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CombiFac(N,cij,fac)
  implicit real(kind=8) (a-h,o-z)
  dimension cij(N+1,N+1),fac(*)

  do i=0,N
    do j=0,i
      cij(j+1,i+1)=fac(i+1)/(fac(j+1)*fac(i-j+1))
    end do
  end do

  Return
End subroutine CombiFac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! save the (I0,I0)-(I0+Iw-1,I0+Iw-1) diagonal block of matrix A into B.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TakDBlk(N,I0,IW,A,B)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N,N), B(IW,IW)

 j=0
 do i=I0, I0+IW-1
   j=j+1
   B(:,j) = A(I0:I0+IW-1,i)
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Put B at the (I0,I0)-(I0+Iw-1,I0+Iw-1) diagonal block of matrix A.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PutDBlk(N,I0,IW,A,B)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: A(N,N), B(IW,IW)

 j=0
 do i=I0, I0+IW-1
   j=j+1
   A(I0:I0+IW-1,i) = B(:,j)
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! lower triangular matrix --> symmetric (mode = 0) or antisymmetric (mode = 1) square matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine lt2sqr(mode,n,t,s)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: t(*),s(n,n)

 k=0
 if(mode == 0) then
   do i=1,N
     do j=1,i-1
       k=k+1
       S(j,i) = T(k)
       S(i,j) = T(k)
     end do
     k=k+1
     S(i,i) = T(k)
   end do
 else if(mode == 1) then
   do i=1,N
     do j=1,i-1
       k=k+1
       S(j,i) =-T(k)
       S(i,j) = T(k)
     end do
     k=k+1
     S(i,i) = T(k)
   end do
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculate a small value Eps = zero in the accuracy of double precision.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double precision function EpsFun(M)
 use Constant, only : Zero, Tenth, One
 implicit real(kind=8) (a-h,o-z)

 ! Calculate N which makes 10^(-N) = zero in the accuracy of double precision.
 eps=One
 N=0
 do while(.true.)
   eps=eps*Tenth
   N=N+1
   if(eps > Zero) cycle
   exit
 end do

 EpsFun=Tenth**(min(N/2,M))

 return
end function EpsFun

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! C = B * op(A), where A is a diagonal matrix with the diagonal elements saved in array D
! Iop:
!  <1 = A^{p}
!   1 = A^{-1}
!   2 = A^{1/2}
!   3 = A^{-1/2}
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DAxMT(N,Iop,p,D,B,C)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: D(N), B(N,N), C(N,N)

 if(Iop == 1) then
   do i=1,N
     if(abs(D(i)) <= Zero) then
       write(*,"(/,' Error in sub. DAxMT: D(A) = 0.')")
       call estop(1)
     end if
     C(:,i) = B(:,i) / D(i)
   end do
 else if(Iop == 2) then
   do i=1,N
     if(D(i) < Zero) then
       write(*,"(/,' Error in sub. DAxMT: D(A) < 0.')")
       call estop(1)
     end if
     C(:,i) = B(:,i) * sqrt(D(i))
   end do
 else if(Iop == 3) then
   do i=1,N
     if(D(i) <= Zero) then
       write(*,"(/,' Error in sub. DAxMT: D(A) <= 0.')")
       call estop(1)
     end if
     C(:,i) = B(:,i) / sqrt(D(i))
   end do
 else if(Iop < 1) then
   do i=1,N
     C(:,i) = B(:,i) * D(i)**p
   end do
 else
   write(*,"(/,' Error in sub. DAxMT: Unknown operation.')")
   call estop(1)
 end if

 Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculate
!   X1 = S^{-1} (if Idx1 = 1),
!   X2 = S^{1/2} (if Idx2 = 1), and
!   X3 = S^{-1/2} (if Idx3 = 1) for a symmetric matrix S.
!
! The eigenvectors and eigenvalues of S are saved in E and A, respectively.
! For scratch Scr(N,M), M = max(2*N,4).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SqIvMt(N,S,Idx1,X1,Idx2,X2,Idx3,X3,E,A,Scr)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: S(N,N),X1(N,N),X2(N,N),X3(N,N),E(N,N),A(N),Scr(N,*)

 E = S
 LWork=(N+N)*max(N,2)
 call DSYEV('V','L', N, E, N, A, Scr, LWork, Info)
 if(INFO /= 0) then
   write(*,"(/,' Error in sub. SqIvMt: Diagnolization failed.')")
   call estop(1)
 end if

 ! X1 = S^{-1}
 if(Idx1 == 1) then
   call DAxMT(N,1,One,A,E,Scr)
   call DGEMM('N','T',N,N,N,One,Scr,N,E,N,Zero,X1,N)
 end if

 ! X2 = S^{1/2}
 if(Idx2 == 1) then
   call DAxMT(N,2,One,A,E,Scr)
   call DGEMM('N','T',N,N,N,One,Scr,N,E,N,Zero,X2,N)
 end if

 ! X3 = S^{-1/2}
 if(Idx3 == 1) then
   call DAxMT(N,3,One,A,E,Scr)
   call DGEMM('N','T',N,N,N,One,Scr,N,E,N,Zero,X3,N)
 end if

 Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! real and imaginary parts of a complex array CA are saved separately into CB
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReImSplit(m,n,ca,cb)
 implicit real(kind=8) (a-h,o-z)
 dimension  :: ca(2,m,n), cb(m,n,2)

 do i=1,n
   do j=1,m
     cb(j,i,1) = ca(1,j,i)
   end do
 end do

 do i=1,n
   do j=1,m
     cb(j,i,2) = ca(2,j,i)
   end do
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! lower triangular complex matrix --> hermitian square complex matrix
!
! if splt is .true., the real and imaginary parts are saved separately in S2, and otherwise in S1.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine lt2sqrC(splt,n,t,s1,s2)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: t(2,*),s1(2,n,n),s2(n,n,2)
 logical      :: splt

 k=0
 if(splt) then
   do i=1,N
     do j=1,i-1
       k=k+1
       S2(j,i,1)= T(1,k)
       S2(i,j,1)= T(1,k)
       S2(j,i,2)= T(2,k)  ! the u.t. imag. part is anti-symmetric
       S2(i,j,2)=-T(2,k)
     end do
     k=k+1
     S2(i,i,1)=T(1,k)
     S2(i,i,2)=0.0d0      ! T(2,k) should be zero due to anti-symmetry
   end do
 else
   do i=1,N
     do j=1,i-1
       k=k+1
       S1(1,j,i)= T(1,k)
       S1(2,j,i)= T(2,k)
       S1(1,i,j)= T(1,k)
       S1(2,i,j)=-T(2,k)
     end do
     k=k+1
     S1(1,i,i)=T(1,k)
     S1(2,i,i)=0.0d0      ! T(2,k) should be zero due to anti-symmetry
   end do
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Single the dimension of matrix A.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SglMat(M,N,A,B)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension :: A(2,M,2,N),B(M,N)

 do i = 1, N
   do j = 1, M
     B(j,i) = A(1,j,1,i) + A(2,j,2,i)
   end do
 end do

 Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Double the dimension of matrix A
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DblMat(M,N,A,B)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension :: A(M,N),B(2,M,2,N)

 B = Zero

 do i = 1, N
   do j = 1, M
     B(1,j,1,i) = A(j,i)
     B(2,j,2,i) = A(j,i)
   end do
 end do

 Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 2-component density matrix calculation: P = C O C^H = (A*A^T + B*B^T) + i(B*A^T - A*B^T), where C = A + iB.
!
! Cccupation numbers of GHF/GKS spinor can be either 0 or 1, and only the first NOCC occupied spinors are calculated
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine calc_den2c(nbas,nmo,cmo2,nocc,den)
 use Constant, only : Zero, One
 implicit double precision (a-h,o-z)
 dimension cmo2(nbas,nmo,2), den(nbas,nbas,2)

 call MatMult(3,nbas,nocc,nbas,One,Zero,cmo2(1,1,1),cmo2(1,1,1),den(1,1,1))
 call MatMult(3,nbas,nocc,nbas,One, One,cmo2(1,1,2),cmo2(1,1,2),den(1,1,1))

 call MatMult(3,nbas,nocc,nbas,One,Zero,cmo2(1,1,1),cmo2(1,1,2),den(1,1,2))
 call MatMult(3,nbas,nocc,nbas,One,-One,cmo2(1,1,2),cmo2(1,1,1),den(1,1,2))

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Combine the real and imaginary matrices into a complex matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReImComb(m,n,ca,cb)
 implicit real(kind=8) (a-h,o-z)
 dimension  :: ca(m,n,2), cb(2,m,n)

 do i=1,n
   do j=1,m
     cb(1,j,i) = ca(j,i,1)
   end do
 end do

 do i=1,n
   do j=1,m
     cb(2,j,i) = ca(j,i,2)
   end do
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! B = inv (A), where A and B are complex square matrices.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine InvMatC(N,A,B,IW,INFO)
 implicit real(kind=8) (a-h,o-z)
 dimension :: A(2*N*N),B(2*N*N),IW(N)

 call UnitMC(N,B)
 call ZGESV( N, N, A, N, IW, B, N, INFO)

 Return

 contains

 Subroutine UnitMC(N,A)
  use Constant, only : Zero, One
  implicit real(kind=8) (a-h,o-z)
  ! Complex unit matrix A(N,N)
  dimension :: A(2,N,N)

  A = Zero
  Do i=1,N
    A(1,i,i)=One
  end Do

  Return
 End Subroutine UnitMC

End Subroutine InvMatC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Complex version of sub. DAxMT.
! C = B * op(A), where A is a diagonal real matrix with the diagonal elements saved in array D, and B and C are complex.
! Iop:
!  <1 = A^{p}
!   1 = A^{-1}
!   2 = A^{1/2}
!   3 = A^{-1/2}
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DAxMTC(N,Iop,p,D,B,C,INFO)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 Dimension :: D(N), B(2,N,N), C(2,N,N)

 INFO = 0
 N2 = N + N

 if(Iop == 1) then
   do i=1,N
     if(abs(D(i)) <= Zero) then
       write(*,"(/,' Error in sub. DAxMTC: D(A) = 0.')")
       INFO = -1
       exit
     end if
     call AScale(N2,(One/D(i)),B(1,1,i),C(1,1,i))
   end do
 else if(Iop == 2) then
   do i=1,N
     if(D(i) < Zero) then
       write(*,"(/,' Error in sub. DAxMTC: D(A) < 0.')")
       INFO = -2
       exit
     end if
     call AScale(N2,sqrt(D(i)),B(1,1,i),C(1,1,i))
   end do
 else if(Iop == 3) then
   do i=1,N
     if(D(i) <= Zero) then
       write(*,"(/,' Error in sub. DAxMTC: D(A) <= 0.')")
       INFO = -3
       exit
     end if
     call AScale(N2,(One/sqrt(D(i))),B(1,1,i),C(1,1,i))
   end do
 else if(Iop < 1) then
   do i=1,N
     call AScale(N2,(D(i)**p),B(1,1,i),C(1,1,i))
   end do
 else
   write(*,"(/,' Error in sub. DAxMTC: Unknown operation.')")
   INFO = Iop
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Complex version of sub. MatMult.
! C = beta * C + alpha * op(A) * op(B), where C(MxN), op(A)(MxL), and op(B)(LxN).
! ModA & ModB can be 1 (op = Normal), 2 (op = Transpose), or 3 (op = Transpose & Conjugation).
! Note that alpha and beta are real instead of complex values.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine MatMultC(ModA,ModB,M,L,N,alpha,beta,A,B,C)
 use Constant, only : Zero
 implicit real(kind=8) (a-h,o-z)
 dimension         :: A(2,*),B(2,*),C(2,*)
 Character*1       :: Cha(3)
 Save Cha
 Data Cha/'N','T','C'/
 complex(kind=8)   :: aa, bb

 aa=Cmplx(alpha,Zero)
 bb=Cmplx(beta,Zero)

 LDC=M
 if(ModA == 1 .and. ModB == 1)then
   LDA=M
   LDB=L
   call zgemm(Cha(ModA),Cha(ModB),M,N,L,aa,A,LDA,B,LDB,bb,C,LDC)
 else if(ModA /= 1 .and. ModB == 1)then
   LDA=L
   LDB=L
   call zgemm(Cha(ModA),Cha(ModB),M,N,L,aa,A,LDA,B,LDB,bb,C,LDC)
 else if(ModA == 1 .and. ModB /= 1)then
   LDA=M
   LDB=N
   call zgemm(Cha(ModA),Cha(ModB),M,N,L,aa,A,LDA,B,LDB,bb,C,LDC)
 else if(ModA /= 1 .and. ModB /= 1)then
   LDA=L
   LDB=N
   call zgemm(Cha(ModA),Cha(ModB),M,N,L,aa,A,LDA,B,LDB,bb,C,LDC)
 end if

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Complex version of X = S^(1/2) if Indx = 0 or X = S^(-1/2) if Indx = 1.
! S must be hermitian, othewise the eigenvalues are complex.
! The eigenvectors and eigenvalues of S are saved in E and A, respectively.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SqrtMCh(N,Indx,S,X,E,A,Scr,INFO)
 use Constant, only : Zero, One
 implicit real(kind=8) (a-h,o-z)
 dimension  ::  S(2,N,N),X(2,N,N),E(2,N,N),A(N),Scr(2,N,N)

 LWork=N*N
 E = S
 call ZHEEV('V','L', N, E, N, A, Scr, LWork, X, Info)
 If (INFO /= 0) return

 ! eigenvalues saved in A are p.d.
 call DAxMTC(N,Indx+2,1.0d0,A,E,Scr,Info)
 If (INFO /= 0) return

 call MatMultC(1,3,N,N,N,One,Zero,Scr,E,X)

 Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! compress doubled complex density matrices Pr + i * Pi into single-basis real ones Po, Px, Py, and Pz.
!               i1  i2           i1  i2 
!   j1        | Wo  Wy|        | Wz  Wx|
!   j2        |-Wy  W0|  + I * | Wx -Wz|
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SGLCOMPW(N,Pr,Pi,Po,Px,Py,Pz)
 implicit real(kind=8) (a-h,o-z)
 dimension     :: Pr(2,N,2,N),Pi(2,N,2,N),Po(N,N),Px(N,N),Py(N,N),Pz(N,N)

 do I=1,N
   do J=1,N
     Po(J,I) = Pr(1,J,1,I) + Pr(2,J,2,I)
     Py(J,I) = Pr(1,J,2,I) - Pr(2,J,1,I)
   end do
 end do

 do I=1,N
   do J=1,N
     Px(J,I) = Pi(2,J,1,I) + Pi(1,J,2,I)
     Pz(J,I) = Pi(1,J,1,I) - Pi(2,J,2,I)
   end do
 end do

 return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Modified STForc for APT.
! Trace(P * dO/dl), where O = S or D, and dX collects the first-order derivatives of all the atoms.
!
! In block 1 in D, the taken values for    For ifx=idx, blocks 3 and 2
! Atom-I must be multiplied by -1.         in -S are also needed.
!
!    |-------|---|---------|               |-------|---|---------|
!    |       |///|         |               |       |   |         |
!    |   0   |/1/|    0    |               |   0   | 0 |    0    |
!    |       |///|         |               |       |   |         |
!    |-------|---|---------|               |-------|---|---------|
!    |///1///| I |////2////|               |   0   |/3/|////2////|
!    |-------|---|---------|               |-------|---|---------|
!    |       |///|         |               |       |///|         |
!    |       |///|         |               |       |///|         |
!    |   0   |/2/|    0    |               |   0   |/2/|    0    |
!    |       |///|         |               |       |///|         |
!    |-------|---|---------|               |-------|---|---------|
!
! NABC(1,:) : starting position of primitive functions for each atom
! NABC(2,:) : number of primitive functions for each atom
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine STAPT(NAtm,ifx,idx,N,apt,NABC,PX,SX,DX)
 Implicit Real*8(A-H,O-Z)
 Dimension apt(3,3,NAtm), NABC(2,NAtm), PX(N,N), SX(N,N), DX(N,N)

 Do I = 1, NAtm
   IS = NABC(1,I)
   IW = NABC(2,I)
   ISW= IS + IW
   Do j = IS, ISW-1
     ! left & upper blocks
     If(IS > 1) then
       Lenth = IS-1
       X = dotx(Lenth,PX(1,j),DX(1,j))
       apt(ifx,idx,I) = apt(ifx,idx,I) - X - X
     end if
     ! atomic 1-center blocks
     if(ifx == idx) then
       X = dotx(IW,PX(IS,j),SX(IS,j))
       apt(ifx,idx,I) = apt(ifx,idx,I) - X
     end if
     ! right & lower blocks
     If(ISW <= N) then
       Lenth = N-ISW+1
       X = dotx(Lenth,PX(ISW,j),DX(ISW,j))
       apt(ifx,idx,I) = apt(ifx,idx,I) + X + X
       if(ifx == idx) then
         X = dotx(Lenth,PX(ISW,j),SX(ISW,j))
         apt(ifx,idx,I) = apt(ifx,idx,I) - X - X
       end if
     end if
   end do
 end do

 Return
End subroutine STAPT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Modified STForc for DDD.
!
! NABC(1,:) : starting position of primitive functions for each atom
! NABC(2,:) : number of primitive functions for each atom
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine STDDD(NAtm,N,ddd,NABC,PX,SX,DX)
 Implicit Real*8(A-H,O-Z)
 Dimension ddd(NAtm), NABC(2,NAtm), PX(N,N), SX(N,N), DX(N,N)

 Do I = 1, NAtm
   IS = NABC(1,I)
   IW = NABC(2,I)
   ISW= IS + IW
   Do j = IS, ISW-1
     ! left & upper blocks
     If(IS > 1) then
       Lenth = IS-1
       X = dotx(Lenth,PX(1,j),DX(1,j))
       ddd(I) = ddd(I) - X - X
     end if
     ! atomic 1-center blocks
     if(ifx == idx) then
       X = dotx(IW,PX(IS,j),SX(IS,j))
       ddd(I) = ddd(I) - X
     end if
     ! right & lower blocks
     If(ISW <= N) then
       Lenth = N-ISW+1
       X = dotx(Lenth,PX(ISW,j),DX(ISW,j))
       ddd(I) = ddd(I) + X + X
       if(ifx == idx) then
         X = dotx(Lenth,PX(ISW,j),SX(ISW,j))
         ddd(I) = ddd(I) - X - X
       end if
     end if
   end do
 end do

 Return
End subroutine STDDD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! reorder three values to get |c| <= |b| <= |a|
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine reord3abs(v3)
  implicit real*8(A-H,O-Z)
  dimension v3(3)

  ! a = max(abs(v3))
  j = 1
  a = v3(j)
  do i = 2, 3
    if(abs(v3(i)) > abs(a)) then
      j = i
      a = v3(j)
    end if
  end do

  ! exchange v3(j) and v3(3)
  if(j /= 3) then
    v3(j) = v3(3)
    v3(3) = a
  end if

  ! order v3(1) and v3(2)
  if(abs(v3(1)) > abs(v3(2))) then
    a = v3(1)
    v3(1) = v3(2)
    v3(2) = a
  end if

  return
end subroutine reord3abs

!--- END
