MODULE Constant

  character*5            :: Ver="1.2.0"
  character*12           :: Date="May 14, 2024"

  integer,parameter      :: maxlq=5   ! max LQ number for MOLDEN (5=H)
  ! integer,parameter      :: maxlq_G=6 ! max LQ number for Gaussian FCHK (6=I)
  integer,parameter      :: maxlq_G=maxlq  ! not used!
  integer,parameter      :: maxcn=50  ! max contraction

  character*1,parameter  :: uqnam(0:maxlq)=(/'S','P','D','F','G','H'/)
  character*1,parameter  :: lqnam(0:maxlq)=(/'s','p','d','f','g','h'/)
  character*2,parameter  :: t33nam(6)=(/'XX','XY','YY','XZ','YZ','ZZ'/)
  character*3,parameter  :: pronam(6)=(/'DIP','APT','DDD','ED ','CD ','EFG'/)

  real(kind=8),parameter :: exptol=1.0d-8
  real(kind=8),parameter :: occtol=1.0d-10
  real(kind=8),parameter :: pi=acos(-1.0d0)
  real(kind=8),parameter :: dpi=pi+pi
  real(kind=8),parameter :: fpi=dpi+dpi
  real(kind=8),parameter :: Zero=0.0d0
  real(kind=8),parameter :: Tenth=0.1d0
  real(kind=8),parameter :: Quarter=0.25d0
  real(kind=8),parameter :: Half=0.5d0
  real(kind=8),parameter :: One=1.0d0
  real(kind=8),parameter :: Two=2.0d0
  real(kind=8),parameter :: Three=3.0d0
  real(kind=8),parameter :: au2wvn=219474.63137d0
  real(kind=8),parameter :: au2ang=0.52917720859d0
  real(kind=8),parameter :: clight=137.035999074d0
  real(kind=8),parameter :: clight_si=2.99792458d10
  real(kind=8),parameter :: au2deb=2.541746231d0
  real(kind=8),parameter :: f_4cc=4.0d0*clight*clight
  real(kind=8),parameter :: f_2cc=2.0d0*clight*clight
  real(kind=8),parameter :: s_2cc=sqrt(f_2cc)

  contains

  real(kind=8) function atom_nucrad(IZA)
  ! Nuclear charge radii.
  !
  ! IZA < 110: the nuclear charge radii are taken from Visscher-Dyall (in a.u.)
  !   Visscher and Dyall, At. Data and Nucl. Data Tables 67, 207 (1997).
  !
  ! IZA>= 110: r0 = 0.57 + 0.836 * A^1/3 (in fm), where the isotope mass number A is determined by IZA
  !   according to the relationship A(Z) = 0.004467 * Z^2 + 2.163 * Z - 1.168. See
  !   Appendix A in D. Andrae, Phys. Rep. 336, 414 (2000).
  !   and
  !   D. Andrae, Nuclear charge density distributions in quantum chemistry, in Relativistic Electronic
  !   Structure Theory, Part 1: Fundamentals, P. Schwerdtfeger Ed., Theoretical and Computational
  !   Chemistry, Vol. 11, Elsevier, 2002.

    implicit none
    integer, intent(in)  :: IZA
    real(kind=8) :: r0, A

    r0 = Zero
    select case(IZA)
      case(001)                 ! H
        r0 = 2.6569547399d-5
      case(002)                 ! He
        r0 = 3.5849373401d-5
      case(003)                 ! Li
        r0 = 4.0992133976d-5
      case(004)                 ! Be
        r0 = 4.3632829651d-5
      case(005)                 ! B
        r0 = 4.5906118608d-5
      case(006)                 ! C
        r0 = 4.6940079496d-5
      case(007)                 ! N
        r0 = 4.8847128967d-5
      case(008)                 ! O
        r0 = 5.0580178957d-5
      case(009)                 ! F
        r0 = 5.2927138943d-5
      case(010)                 ! Ne
        r0 = 5.3654104231d-5
      case(011)                 ! Na
        r0 = 5.5699159416d-5
      case(012)                 ! Mg
        r0 = 5.6341070732d-5
      case(013)                 ! Al
        r0 = 5.8165765928d-5
      case(014)                 ! Si
        r0 = 5.8743802504d-5
      case(015)                 ! P
        r0 = 6.0399312923d-5
      case(016)                 ! S
        r0 = 6.0927308666d-5
      case(017)                 ! Cl
        r0 = 6.2448101115d-5
      case(018)                 ! Ar
        r0 = 6.4800211825d-5
      case(019)                 ! K
        r0 = 6.4346167051d-5
      case(020)                 ! Ca
        r0 = 6.4800211825d-5
      case(021)                 ! Sc
        r0 = 6.6963627201d-5
      case(022)                 ! Ti
        r0 = 6.8185577480d-5
      case(023)                 ! V
        r0 = 6.9357616830d-5
      case(024)                 ! Cr
        r0 = 6.9738057221d-5
      case(025)                 ! Mn
        r0 = 7.0850896638d-5
      case(026)                 ! Fe
        r0 = 7.1212829817d-5
      case(027)                 ! Co
        r0 = 7.2273420879d-5
      case(028)                 ! Ni
        r0 = 7.1923970253d-5
      case(029)                 ! Cu
        r0 = 7.3633018675d-5
      case(030)                 ! Zn
        r0 = 7.3963875193d-5
      case(031)                 ! Ga
        r0 = 7.5568424848d-5
      case(032)                 ! Ge
        r0 = 7.7097216161d-5
      case(033)                 ! As
        r0 = 7.7394645153d-5
      case(034)                 ! Se
        r0 = 7.8843427408d-5
      case(035)                 ! Br
        r0 = 7.8558604038d-5
      case(036)                 ! Kr
        r0 = 7.9959560033d-5
      case(037)                 ! Rb
        r0 = 8.0233033713d-5
      case(038)                 ! Sr
        r0 = 8.1040799081d-5
      case(039)                 ! Y
        r0 = 8.1305968993d-5
      case(040)                 ! Zr
        r0 = 8.1569159980d-5
      case(041)                 ! Nb
        r0 = 8.2347219923d-5
      case(042)                 ! Mo
        r0 = 8.3607614434d-5
      case(043)                 ! Tc
        r0 = 8.3607614434d-5
      case(044)                 ! Ru
        r0 = 8.4585397905d-5
      case(045)                 ! Rh
        r0 = 8.4825835954d-5
      case(046)                 ! Pd
        r0 = 8.5537941156d-5
      case(047)                 ! Ag
        r0 = 8.5772320442d-5
      case(048)                 ! Cd
        r0 = 8.7373430179d-5
      case(049)                 ! In
        r0 = 8.7596760865d-5
      case(050)                 ! Sn
        r0 = 8.8694413774d-5
      case(051)                 ! Sb
        r0 = 8.8910267995d-5
      case(052)                 ! Te
        r0 = 9.0801452955d-5
      case(053)                 ! I
        r0 = 9.0181040290d-5
      case(054)                 ! Xe
        r0 = 9.1209776425d-5
      case(055)                 ! Cs
        r0 = 9.1412392742d-5
      case(056)                 ! Ba
        r0 = 9.2410525664d-5
      case(057)                 ! La
        r0 = 9.2607247118d-5
      case(058)                 ! Ce
        r0 = 9.2803027311d-5
      case(059)                 ! Pr
        r0 = 9.2997877424d-5
      case(060)                 ! Nd
        r0 = 9.3576955934d-5
      case(061)                 ! Pm
        r0 = 9.3768193375d-5
      case(062)                 ! Sm
        r0 = 9.5082839751d-5
      case(063)                 ! Eu
        r0 = 9.5267329183d-5
      case(064)                 ! Gd
        r0 = 9.6177915369d-5
      case(065)                 ! Tb
        r0 = 9.6357719009d-5
      case(066)                 ! Dy
        r0 = 9.6892647152d-5
      case(067)                 ! Ho
        r0 = 9.6892647152d-5
      case(068)                 ! Er
        r0 = 9.7943009317d-5
      case(069)                 ! Tm
        r0 = 9.8115626740d-5
      case(070)                 ! Yb
        r0 = 9.8968651305d-5
      case(071)                 ! Lu
        r0 = 9.9137288835d-5
      case(072)                 ! Hf
        r0 = 9.9970978172d-5
      case(073)                 ! Ta
        r0 = 1.0013585755d-4
      case(074)                 ! W
        r0 = 1.0062688070d-4
      case(075)                 ! Re
        r0 = 1.0111259523d-4
      case(076)                 ! Os
        r0 = 1.0191070333d-4
      case(077)                 ! Ir
        r0 = 1.0206865731d-4
      case(078)                 ! Pt
        r0 = 1.0238293593d-4
      case(079)                 ! Au
        r0 = 1.0269507292d-4
      case(080)                 ! Hg
        r0 = 1.0346628039d-4
      case(081)                 ! Tl
        r0 = 1.0392291259d-4
      case(082)                 ! Pb
        r0 = 1.0437511130d-4
      case(083)                 ! Bi
        r0 = 1.0452487744d-4
      case(084)                 ! Po
        r0 = 1.0452487744d-4
      case(085)                 ! At
        r0 = 1.0467416660d-4
      case(086)                 ! Rn
        r0 = 1.0642976299d-4
      case(087)                 ! Fr
        r0 = 1.0657317899d-4
      case(088)                 ! Ra
        r0 = 1.0700087100d-4
      case(089)                 ! Ac
        r0 = 1.0714259349d-4
      case(090)                 ! Th
        r0 = 1.0784503195d-4
      case(091)                 ! Pa
        r0 = 1.0770535752d-4
      case(092)                 ! U
        r0 = 1.0867476102d-4
      case(093)                 ! Np
        r0 = 1.0853744903d-4
      case(094)                 ! Pu
        r0 = 1.0949065967d-4
      case(095)                 ! Am
        r0 = 1.0935561268d-4
      case(096)                 ! Cm
        r0 = 1.0989359973d-4
      case(097)                 ! Bk
        r0 = 1.0989359973d-4
      case(098)                 ! Cf
        r0 = 1.1042580946d-4
      case(099)                 ! Es
        r0 = 1.1055797721d-4
      case(100)                 ! Fm
        r0 = 1.1121362374d-4
      case(101)                 ! Md
        r0 = 1.1134373034d-4
      case(102)                 ! No
        r0 = 1.1147350119d-4
      case(103)                 ! Lr
        r0 = 1.1186082063d-4
      case(104)                 ! Rf
        r0 = 1.1173204420d-4
      case(105)                 ! Db
        r0 = 1.1186082063d-4
      case(106)                 ! Sg
        r0 = 1.1198926979d-4
      case(107)                 ! Bh
        r0 = 1.1186082063d-4
      case(108)                 ! Hs
        r0 = 1.1224519460d-4
      case(109)                 ! Mt
        r0 = 1.1237267433d-4
      case(110:)  ! Ds -
        A = 4.467d-3 * dble(IZA*IZA) + 2.163d0 * dble(IZA) -1.168d0
        A = dble(NINT(A))
        r0 = 0.57d0 + 0.836d0*A**(1.0d0/3.0d0)
        ! fm --> a.u.
        r0 = r0 * 1.0d-5 / au2ang
    end select

    atom_nucrad = r0

    return
  end function atom_nucrad

END MODULE Constant

!--- END

