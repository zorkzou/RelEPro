!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! population
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine population(iout,natom,np,irelc,maxtyp,iza,za,mxla,ialbs,nshlla,smat,pmat,scr)
 use Constant, only : Zero,One
 implicit real(kind=8) (a-h,o-z)
 dimension :: iza(*), za(*), mxla(*), ialbs(*), nshlla(*), smat(np*np), pmat(np*np), scr(np*np,4)

 write(iout,"(//,23x,'.',23('-'),'.',/,23x,'|  P O P U L A T I O N  |',/,23x,'`',23('-'),'`')")

 ! S^1/2 --> Scr1
 call sqivmt(np,smat,0,dum,1,scr(1,1),0,dum,scr(1,2),scr(1,3),scr(1,4))

 ! P or Rel(P) --> Scr2
 if(irelc == 1) then
   scr(:,2) = pmat
 else if(irelc == 2) then
   call SglMat(np,np,pmat,scr(1,2))
 end if

 ! Lowdin charge
 call dgemm('n','n',np,np,np,One,scr(1,2),np,scr(1,1),np,Zero,scr(1,3),np)
 call acharge(iout,1,natom,np,maxtyp,iza,za,mxla,ialbs,nshlla,scr(1,1),scr(1,3),scr(1,4),scr(np+1,4))

 ! Mulliken charge
 call acharge(iout,2,natom,np,maxtyp,iza,za,mxla,ialbs,nshlla,smat,scr(1,2),scr(1,4),scr(np+1,4))

 ! Mulliken's BO and generalized Mayer's BO
 ! Re(D) = Re(P * S) --> Scr3
 call dgemm('n','n',np,np,np,One,scr(1,2),np,smat,np,Zero,scr(1,3),np)
 ! Im(D) = Im(P * S) --> Scr1
 if(irelc == 2) call dgemm('n','n',np,np,np,One,pmat(np*np*4+1),np,smat,np,Zero,scr(1,1),np)
 call mubo(iout,natom,np,irelc,iza,mxla,ialbs,nshlla,scr(1,2),smat,scr(1,4))
 call gmbo(iout,natom,np,irelc,iza,mxla,ialbs,nshlla,scr(1,3),scr(1,1),scr(1,4))

 return

 contains

 ! Lowdin (itype = 1) or Mulliken (itype = 2) charge
 subroutine acharge(iout,itype,natom,np,mxl,iza,za,mxla,ialbs,nshlla,sm,pm,tmp,chl)
  use Constant, only : maxlq, Zero
  implicit real(kind=8) (a-h,o-z)
  dimension :: iza(natom), za(natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom), &
    sm(np,np), pm(np,np), tmp(np), chl(-1:mxl,natom)
  character*3   :: Elm

  ! diag(sm * pm) --> tmp
  ! here sm must be a square matrix, thus tmp(i) = dot(sm(:,i),pm(:,i))
  do i = 1, np
    tmp(i) = dotx(np,sm(1,i),pm(1,i))
  end do

  chl = Zero
  do i = 1, natom
    chl(-1,i) = za(i)
    do lq = 0, mxla(i)
      j1 = ialbs(lq,i)
      j2 = j1 + (lq+lq+1)*nshlla(lq,i) - 1
      do j = j1, j2
        chl(lq,i) = chl(lq,i) + tmp(j)
      end do
      chl(-1,i) = chl(-1,i) - chl(lq,i)
    end do
  end do

  if(itype == 1) then
    write(iout,"(/,' Lowdin charge')")
  else if(itype == 2) then
    write(iout,"(/,' Mulliken charge')")
  end if
  write(iout,"(1x,106('-'))")
  write(iout,"(3x,'No.   Atom   ZA        QA    QA for L =',i5,10(5x,i5))") (i, i=0,mxl)
  write(iout,"(1x,106('-'))")
  do i = 1, natom
    call ElemZA(1,Elm,iza(i))
    write(iout,"(i5,5x,a3,i5,f10.3,9x,10f10.4)") i, Elm, iza(i), chl(-1,i), chl(0:mxla(i),i)
  end do
  write(iout,"(1x,106('-'))")

  return
 end subroutine acharge

 ! Mulliken's BO
 subroutine mubo(iout,natom,np,irelc,iza,mxla,ialbs,nshlla,pm,sm,bo)
  use Constant, only : maxlq, Zero, lqnam
  implicit real(kind=8) (a-h,o-z)
  parameter (tolbo=5.0d-3)
  dimension :: iza(natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom), &
    pm(np,np), sm(np,np), bo(0:maxlq,natom,0:maxlq,natom)
  character*3   :: Elmi, Elmj

  if(natom < 2) return

  write(iout,"(/,' Mulliken''s Bond Order')")
  write(iout,"(1x,106('-'))")

  bo = Zero
  do i = 1, natom-1
    i1 = ialbs(0,i) - 1
    do iq = 0, mxla(i)
      i2 = (iq+iq+1)*nshlla(iq,i)
      do i3 = 1, i2
        i1 = i1 + 1

        do j = i+1, natom
          j1 = ialbs(0,j) - 1
          do jq = 0, mxla(j)
            j2 = (jq+jq+1)*nshlla(jq,j)
            do j3 = 1, j2
              j1 = j1 + 1
              bo(jq,j,iq,i) = bo(jq,j,iq,i) + 2.0d0 * pm(j1,i1) * sm(j1,i1)
            end do
          end do
        end do

      end do
    end do
  end do

  do i = 1, natom-1
    call ElemZA(1,Elmi,iza(i))
    do j = i+1, natom
      call ElemZA(1,Elmj,iza(j))
      tmp = Zero
      do iq = 0, mxla(i)
        do jq = 0, mxla(j)
          tmp = tmp + bo(jq,j,iq,i)
        end do
      end do
      if(abs(tmp) > tolbo) then    ! Mulliken BO can be negative
        ! if(tmp > tolbo) then
        write(iout,"(2(1x,i4,1x,a3),f10.3)",advance='no') i, Elmi, j, Elmj, tmp
        i1 = 0
        do iq = 0, mxla(i)
          do jq = 0, mxla(j)
            if(abs(bo(jq,j,iq,i)) > tolbo) then
              i1 = i1 + 1
              if(i1 == 1) then
                write(iout,"(6x,'Compositions (%)',2x)",advance='no')
              else if(mod(i1-1,5) == 0) then
                write(iout,"(/,52x)",advance='no')
              end if
              write(iout,"(2x,a1,'-',a1,':',f5.1)",advance='no') lqnam(iq), lqnam(jq), 1.0d2*bo(jq,j,iq,i)/tmp
            end if
          end do
        end do
        write(iout,*)
      end if
    end do
  end do

  write(iout,"(1x,106('-'))")

  return
 end subroutine mubo

 ! generalized MBO
 subroutine gmbo(iout,natom,np,irelc,iza,mxla,ialbs,nshlla,dr,di,bo)
  use Constant, only : maxlq, Zero, lqnam
  implicit real(kind=8) (a-h,o-z)
  parameter (tolbo=5.0d-3)
  dimension :: iza(natom), mxla(natom), nshlla(0:maxlq,natom), ialbs(0:maxlq,natom), &
    dr(np,np), di(np,np), bo(0:maxlq,natom,0:maxlq,natom)
  character*3   :: Elmi, Elmj

  if(natom < 2) return

  write(iout,"(/,' (Generalized) Mayer''s Bond Order')")
  write(iout,"(1x,106('-'))")

  ! Q = D .* D^T
  ! Dr will be destroyed!
  do i = 1, np
    do j = 1, i-1
      dr(j,i) = dr(i,j) * dr(j,i)
      if(irelc == 2) dr(j,i) = dr(j,i) + abs(di(i,j) * di(j,i))
      dr(i,j) = dr(j,i)
    end do
  end do

  bo = Zero
  do i = 1, natom-1
    i1 = ialbs(0,i) - 1
    do iq = 0, mxla(i)
      i2 = (iq+iq+1)*nshlla(iq,i)
      do i3 = 1, i2
        i1 = i1 + 1

        do j = i+1, natom
          j1 = ialbs(0,j) - 1
          do jq = 0, mxla(j)
            j2 = (jq+jq+1)*nshlla(jq,j)
            do j3 = 1, j2
              j1 = j1 + 1
              bo(jq,j,iq,i) = bo(jq,j,iq,i) + dr(j1,i1)
            end do
          end do
        end do

      end do
    end do
  end do

  do i = 1, natom-1
    call ElemZA(1,Elmi,iza(i))
    do j = i+1, natom
      call ElemZA(1,Elmj,iza(j))
      tmp = Zero
      do iq = 0, mxla(i)
        do jq = 0, mxla(j)
          tmp = tmp + bo(jq,j,iq,i)
        end do
      end do
      if(abs(tmp) > tolbo) then
        write(iout,"(2(1x,i4,1x,a3),f10.3)",advance='no') i, Elmi, j, Elmj, tmp
        i1 = 0
        do iq = 0, mxla(i)
          do jq = 0, mxla(j)
            if(bo(jq,j,iq,i) > tolbo) then
              i1 = i1 + 1
              if(i1 == 1) then
                write(iout,"(6x,'Compositions (%)',2x)",advance='no')
              else if(mod(i1-1,5) == 0) then
                write(iout,"(/,52x)",advance='no')
              end if
              write(iout,"(2x,a1,'-',a1,':',f5.1)",advance='no') lqnam(iq), lqnam(jq), 1.0d2*bo(jq,j,iq,i)/tmp
            end if
          end do
        end do
        write(iout,*)
      end if
    end do
  end do

  write(iout,"(1x,106('-'))")

  return
 end subroutine gmbo

end

!--- END
