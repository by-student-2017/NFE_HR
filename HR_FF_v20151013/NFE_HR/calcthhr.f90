!----------------------------------------------------------------------
! HISTORY: version
! new avga.f for As3Ca (No.2)                           2013.6.6
! PROJECTED (DECOMPOSED) STATE-DENSITY PROGRAMME
!                                04.3.12, 09.8.30, 11.2.22, 3.15
! 2013/06/15 F90 type :M. Inukai
! 2013/07/19 modified :M. Inukai
! 2013/07/29 Monkhost-Pack type data is OK : M. Inukai
!----------------------------------------------------------------------
! common sentence for Fortran 90
! module constants
! module common
! module allocate_data
!----------------------------------------------------------------------
! Compile option for gfortran
! gfortran -o calcthhr -Wall -pedantic -std=f95 -fbounds-check -O
! -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace
! -g calcthhr.f90 -free -m64
!----------------------------------------------------------------------
module constants
  implicit none
  real(8),parameter:: PI=3.1415926535897932d0
  real(8),parameter:: EV=13.605698
end module constants
!----------------------------------------------------------------------
module common
  implicit none
  real(8):: V
  real(8):: EB,EF,DE,GEE
  integer(4):: NE,NEMAX
  integer(4):: NPX,NPY,NPZ,NPT,MMIN,MAXII,NTK
  integer(4):: ML,MML
  integer(4):: IATO,NATO
  integer(4):: READ_LMAX
  real(8):: COK, fac
  !
contains
  subroutine read_input
    open(13, file="f13")
      read(13,'(4x,i6)') NPX
      read(13,'(4x,i6)') NPY
      read(13,'(4x,i6)') NPZ
      read(13,'(4x,i10)') NPT
    close(13)
    !
    open(14, file="f14")
      read(14,'(8x,f10.5)') COK
      read(14,'(8x,f10.5)') fac
    close(14)
    !
    open(16, file="f16")
      read(16,'(4x,i6)') MMIN
      read(16,'(4x,i6)') MAXII
      read(16,'(4x,i6)') NTK
    close(16)
    !
    open(70, file="f70")
      read(70,'(5x,i5)') IATO
      read(70,'(5x,i5)') NATO
    close(70)
    !
    open(75, file="f75")
      read(75,'(5x,i2)') READ_LMAX
    close(75)
    !
    write(6,'(A)') '===============READ INPUT DATA================'
    write(6,'(A,i10)') 'Number of x axis mesh=',NPX
    write(6,'(A,i10)') 'Number of y axis mesh=',NPY
    write(6,'(A,i10)') 'Number of z axis mesh=',NPZ
    write(6,'(A,i10)') 'Total mesh in tetrahedoron method=', NPT
    write(6,'(A,f10.5)') 'Vabc/Va=', COK
    write(6,'(A,f10.5)') 'V/Vcub =', fac
    write(6,'(A,i10)') 'check data =', MMIN
    write(6,'(A,i10)') 'Number of total K^2 and C^2 data=', MAXII
    write(6,'(A,i10)') 'Number of irreducible total k point=',NTK
    write(6,'(A,i5)') 'INDEPENDENT ATOMS IN UNITCELL=',IATO
    write(6,'(A,i5)') 'TOTAL ATOMS IN UNITCELL=',NATO
    write(6,'(A,i2)') 'LMAX=', READ_LMAX
    return
  end subroutine read_input
end module common
!----------------------------------------------------------------------
module allocate_data
  implicit none
  ! NPT
  integer(4),allocatable:: N1(:),N2(:),N3(:),N4(:),N5(:),N6(:),N7(:),IEX(:)
  real(8),allocatable:: IWW(:)
  ! NPXD, NPYD, NPZD
  integer(4),allocatable:: IC(:,:,:)
  ! MAXII
  integer(4),allocatable:: IAB(:)
  real(8),allocatable:: ES(:)
  ! MAXII,MMIN
  real(8),allocatable:: AKKK(:,:,:),AFFF(:,:,:) ! K^2, C^2
  real(8),allocatable:: AMA1(:,:) ! energy
  real(8),allocatable:: AMA2(:,:),BKKK(:,:),BFFF(:,:) ! energy, K^2, C^2
  real(8),allocatable:: PR(:),WV(:)
  real(8),allocatable:: AMA(:,:),AKKB(:,:), AKK(:,:) ! energy, K^2, C^2
  ! NPT
  real(8),allocatable:: E(:),FK(:),FKB(:),FKC(:) ! energy, K^2, C^2, (K^2)^2
  ! NE
  real(8),allocatable:: EM(:) ! energy
  real(8),allocatable:: CNT(:),SNT(:),FNT(:),FNTB(:),FNTC(:),FNTD(:)
  !
end module allocate_data
!----------------------------------------------------------------------
subroutine initial_data1(NE,SNT,CNT,FNT,FNTB,FNTC,FNTD)
  implicit none
  integer(4):: i
  integer(4),intent(in):: NE
  real(8),intent(inout):: SNT(NE),CNT(NE),FNT(NE),FNTB(NE),FNTC(NE),FNTD(NE)
  do i=1, NE
    SNT(i) = 0.d0
    CNT(i) = 0.d0
    FNT(i) = 0.d0
    FNTB(i)= 0.d0
    FNTC(i)= 0.d0
    FNTD(i)= 0.d0
  end do
  return
end subroutine initial_data1
!----------------------------------------------------------------------
subroutine initial_data2(NPXD,NPYD,NPZD,IC,NPT,N1,N2,N3,N4,N5,N6,N7,IEX,IWW)
  implicit none
  integer(4):: i,j,k,m
  integer(4),intent(in):: NPXD,NPYD,NPZD,NPT
  integer(4),intent(inout):: IC(NPXD,NPYD,NPZD) &
    ,N1(NPT),N2(NPT),N3(NPT),N4(NPT),N5(NPT),N6(NPT),N7(NPT),IEX(NPT)
  real(8),intent(inout)::IWW(NPT)
  do i=1, NPXD
    do j=1, NPYD
      do k=1, NPZD
        IC(i,j,k)=0
      end do
    end do
  end do
  do m=1, NPT
    N1(m)=0
    N2(m)=0
    N3(m)=0
    N4(m)=0
    N5(m)=0
    N6(m)=0
    N7(m)=0
    IEX(m)=0
    IWW(m)=0
  end do
  return
end subroutine initial_data2
!----------------------------------------------------------------------
subroutine initial_data3(MAXII,IAB,ES)
  implicit none
  integer(4)::i
  integer(4)::MAXII,IAB(MAXII)
  real(8):: ES(MAXII)
  do i=1, MAXII
    IAB(i)=0
    ES(i) =0.d0
  end do
  return
end subroutine initial_data3
!----------------------------------------------------------------------

subroutine TGEN(NPXD,NPYD,NPZD,IC,N1,N2,N3,N4,N5,N6,N7)
  use common
  implicit none
!---------------------------------------------------------------
!   *   DEFINE TETRAHEDRA ASSOCIATED WITH EACH K-POINT         *
!   *   THE K-SPACE MESH IMPLIED SHOULD BE CONSISTENT WITH THAT*
!   *   GENERATED BY SUBROUTINE <MESH> IN <STR>                *
!   *   SIMPLE CUBIC                                           *
!---------------------------------------------------------------
  integer(4),intent(in):: NPXD,NPYD,NPZD
  integer(4),intent(inout):: IC(NPXD,NPYD,NPZD)
  integer(4),intent(inout):: N1(NPT),N2(NPT),N3(NPT),N4(NPT),N5(NPT),N6(NPT),N7(NPT)
  !
  integer(4):: X,Y,Z
  integer(4):: NPXM,NPYM,NPZM
  integer(4):: NP
  real(8):: DKX,DKY,DKZ
  !
  NPX=NPXD
  NPY=NPYD
  NPZ=NPZD
  !
  NPXM=NPX-1
  NPYM=NPY-1
  NPZM=NPZ-1
  !
  DKX=1.d0/(float(NPXM))
  DKY=1.d0/(float(NPYM))
  DKZ=1.d0/(float(NPZM))
  V=  DKX*DKY*DKZ / 3.0
  write(6,*) ' start!! PNT=', 'V=',V
!
! be careful! h,k,l start at 0, but fortran array start at 1.
! Please, remember X:1 = h:0, Y:1 = k:0 , Z:1 = l:0
!
  NP=0
  ! write(6,*) NPX, NPY, NPZ
  do X=1,NPX
    do Y=1,NPY
      do Z=1,NPZ
        NP=NP+1
        IC(X,Y,Z)=NP
      end do
    end do
  end do
  !
  ! write(6,*) ' total k-point',NP
  !
  NP=0
  DO X=1,NPX-1
    DO Y=1,NPY-1
      DO Z=1,NPZ-1
        NP=IC(X,Y,Z)
        N1(NP)=IC(X+1,Y  ,Z  )
        N2(NP)=IC(X+1,Y+1,Z  )
        N3(NP)=IC(X+1,Y+1,Z+1)
        N4(NP)=IC(X+1,Y  ,Z+1)
        N5(NP)=IC(X  ,Y  ,Z+1)
        N6(NP)=IC(X  ,Y+1,Z+1)
        N7(NP)=IC(X  ,Y+1,Z  )
      end do
    end do
  end do
  return
end subroutine TGEN
!----------------------------------------------------------------------
subroutine FDENS(NPXD,NPYD,NPZD,IC,NPT,N1,N2,N3,N4,N5,N6,N7)
  implicit none
  integer(4),intent(in):: NPXD, NPYD, NPZD, NPT
  integer(4),intent(inout):: IC(NPXD,NPYD,NPZD)
  integer(4),intent(inout):: N1(NPT),N2(NPT),N3(NPT),N4(NPT),N5(NPT),N6(NPT),N7(NPT)
  integer(4):: I,J,K
  integer(4):: NPX,NPY,NPZ
  integer(4):: LP
  integer(4):: I0,I1,I2,I3,I4,I5,I6,I7
  !
  NPX=NPXD
  NPY=NPYD
  NPZ=NPZD
  LP=0
  DO I=1,NPX-1
    DO J=1,NPY-1
      DO K=1,NPZ-1
        ! be careful! Tihs method is calculate tetrahedron including I0+1 range.
        ! NPX-1 is OK.
        ! NP=NP+1
        LP=IC(I,J,K)
        I0=LP
        I1=N1(LP)
        I2=N2(LP)
        I3=N3(LP)
        I4=N4(LP)
        I5=N5(LP)
        I6=N6(LP)
        I7=N7(LP)
        !
        CALL PDNS(I0,I1,I2,I3)
        CALL PDNS(I0,I1,I3,I4)
        CALL PDNS(I0,I3,I4,I5)
        CALL PDNS(I0,I3,I5,I6)
        CALL PDNS(I0,I3,I6,I7)
        CALL PDNS(I0,I2,I3,I7)
        !
      end do
    end do
  end do
  return
end subroutine FDENS
!----------------------------------------------------------------------

subroutine PDNS(I1,I2,I3,I4)
  use constants
  use common
  use allocate_data
  implicit none
!------------------------------------------------------------------
!   *   CALCULATE STATE DENSITY AND NUMBER OF STATES FOR A SINGLE *
!   *   TETRAHEDRON                                               *
!   *   CN : STATE DENSITY                                        *
!   *   SN : NUMBER OF STATES                                     *
!------------------------------------------------------------------
  !
  integer(4),intent(in):: I1,I2,I3,I4
  !
  real(8):: S(4),F(4),FB(4),FC(4),G(4)
  real(8):: DELTA,AM,D1,D2,D4
  real(8):: V1,DEE,FRAC
  real(8):: CN,SN,SSS,CN1,CN2,DF,DF1,DF2
  real(8):: FN,FNB,FNC,FND
  real(8):: EE
  integer(4):: LE,NE1,NE2
  !
  G(1)=1.0
  G(2)=1.0
  G(3)=1.0
  G(4)=1.0
  !
  S(1)=E(I1)
  S(2)=E(I2)
  S(3)=E(I3)
  S(4)=E(I4)
  F(1)=FK(I1)
  F(2)=FK(I2)
  F(3)=FK(I3)
  F(4)=FK(I4)
  FB(1)=FKB(I1)
  FB(2)=FKB(I2)
  FB(3)=FKB(I3)
  FB(4)=FKB(I4)
  FC(1)=FKC(I1)
  FC(2)=FKC(I2)
  FC(3)=FKC(I3)
  FC(4)=FKC(I4)
  !
  CALL SORT(S,F,FB,FC)
  !
  if( S(1)>=(EM(NE)) ) return
  EB=EM(1)
  NE1=int ( (S(1)-EB)/DE )
  if( S(1)<EB ) NE1=-1
  NE1=NE1+2
  NE2=NEMAX
  DO 25 LE=NE1,NE2
    EE=EM(LE)
    if( EE <= S(1) ) goto 25
    if( EE >= S(4) ) goto 20
    if( EE <= S(2) ) goto 21
    if( EE >= S(3) ) goto 22
! CASE_2 E2 <E< E3
    DELTA=S(4)+S(3)-S(2)-S(1)
    AM=(S(4)*S(3)-S(2)*S(1))/DELTA
    D1=S(1)-AM
    D2=S(2)-AM
    V1=V/DELTA
    DEE=EE-AM
    FRAC=DEE*DEE/D2/D1
    CN=3.D0*V1*(1.D0-FRAC)
    SN=3.D0*V1*DEE-V1*(D1+D2+DEE*FRAC)
    if( ABS(S(2)-S(1)) <= 0.0 ) then
      SSS=3.D0*V*(EE-S(1))/(S(3)-S(1))/(S(4)-S(1))
      CN1=SSS*(S(3)-EE)/(S(3)-S(1))
      CN2=SSS*(S(4)-EE)/(S(4)-S(1))
      ! FN
      DF1=(F(2)-F(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
        +2*(F(3)-F(1))*(EE-S(1))/(S(3)-S(1)) &
        +  (F(4)-F(1))*(EE-S(1))/(S(4)-S(1))
      DF2=(F(2)-F(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
        +  (F(3)-F(1))*(EE-S(1))/(S(3)-S(1)) &
        +2*(F(4)-F(1))*(EE-S(1))/(S(4)-S(1))
      FN=F(1)*CN+ (DF1*CN1+DF2*CN2)/3.D0
      ! FNB
      DF1=(FB(2)-FB(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
        +2*(FB(3)-FB(1))*(EE-S(1))/(S(3)-S(1)) &
        +  (FB(4)-FB(1))*(EE-S(1))/(S(4)-S(1))
      DF2=(FB(2)-FB(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
        +  (FB(3)-FB(1))*(EE-S(1))/(S(3)-S(1)) &
        +2*(FB(4)-FB(1))*(EE-S(1))/(S(4)-S(1))
      FNB=FB(1)*CN+(DF1*CN1+DF2*CN2)/3.D0
      ! FNC
      DF1=(F(2)*FB(2)-F(1)*FB(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
        +2*(F(3)*FB(3)-F(1)*FB(1))*(EE-S(1))/(S(3)-S(1)) &
        +  (F(4)*FB(4)-F(1)*FB(1))*(EE-S(1))/(S(4)-S(1))
      DF2=(F(2)*FB(2)-F(1)*FB(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
        +  (F(3)*FB(3)-F(1)*FB(1))*(EE-S(1))/(S(3)-S(1)) &
        +2*(F(4)*FB(4)-F(1)*FB(1))*(EE-S(1))/(S(4)-S(1))
      FNC=F(1)*FB(1)*CN+(DF1*CN1+DF2*CN2)/3.D0
      ! FND
      DF1=(G(2)*FC(2)-G(1)*FC(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
        +2*(G(3)*FC(3)-G(1)*FC(1))*(EE-S(1))/(S(3)-S(1)) &
        +  (G(4)*FC(4)-G(1)*FC(1))*(EE-S(1))/(S(4)-S(1))
      DF2=(G(2)*FC(2)-G(1)*FC(1))*( (S(3)-EE)/(S(3)-S(1))+ (S(4)-EE)/(S(4)-S(1)) ) &
        +  (G(3)*FC(3)-G(1)*FC(1))*(EE-S(1))/(S(3)-S(1)) &
        +2*(G(4)*FC(4)-G(1)*FC(1))*(EE-S(1))/(S(4)-S(1))
      FND=G(1)*FC(1)*CN+(DF1*CN1+DF2*CN2)/3.D0
      !
      goto 23
    end if
    if( ABS(S(3)-S(2)) <= 0.0) then
      ! FN
      DF=F(2)-F(1)+F(3)-F(1)+(F(4)-F(1))*(EE-S(1))/(S(4)-S(1))
      FN=F(1)*CN+DF*CN/3.D0
      ! FNB
      DF=FB(2)-FB(1)+FB(3)-FB(1)+(FB(4)-FB(1))*(EE-S(1))/(S(4)-S(1))
      FNB=FB(1)*CN+DF*CN/3.D0
      ! FNC
      DF=F(2)*FB(2)-F(1)*FB(1)+F(3)*FB(3)-F(1)*FB(1)+ &
      (F(4)*FB(4)-F(1)*FB(1))*(EE-S(1))/(S(4)-S(1))
      FNC=F(1)*FB(1)*CN+DF*CN/3.D0
      ! FND
      DF=G(2)*FC(2)-G(1)*FC(1)+G(3)*FC(3)-G(1)*FC(1)+ &
      (G(4)*FC(4)-G(1)*FC(1))*(EE-S(1))/(S(4)-S(1))
      FND=G(1)*FC(1)*CN+DF*CN/3.D0
      !
      goto 23
    end if
! CN1=small tri, CN2=big tri
!     CN2=3.D0*V1*(-DDE*DDE/D1+2.D0*DEE-D1) /(S(2)-S(1))
    CN2=3.D0*V *(EE-S(1))**2/(S(2)-S(1))/(S(3)-S(1))/(S(4)-S(1))
    CN1=CN2-CN
    ! FN
    DF2=( (F(2)-F(1))/(S(2)-S(1)) +(F(3)-F(1))/(S(3)-S(1)) &
      +(F(4)-F(1))/(S(4)-S(1)) )*(EE-S(1))
    DF1=(F(2)-F(1))*( (EE-S(1))/(S(2)-S(1))+(S(3)-EE)/(S(3)-S(2)) &
      +(S(4)-EE)/(S(4)-S(2)) ) &
      +( (F(3)-F(1))/(S(3)-S(2))+(F(4)-F(1))/(S(4)-S(2)) )*(EE-S(2))
    FN=F(1)*CN +(DF2*CN2-DF1*CN1)/3.D0
    ! FNB
    DF2=( (FB(2)-FB(1))/(S(2)-S(1)) +(FB(3)-FB(1))/(S(3)-S(1)) &
      +(FB(4)-FB(1))/(S(4)-S(1)) )*(EE-S(1))
    DF1=(FB(2)-FB(1))*( (EE-S(1))/(S(2)-S(1))+(S(3)-EE)/(S(3)-S(2)) &
      +(S(4)-EE)/(S(4)-S(2)) ) &
      +( (FB(3)-FB(1))/(S(3)-S(2))+(FB(4)-FB(1))/(S(4)-S(2)) )*(EE-S(2))
    FNB=FB(1)*CN+(DF2*CN2-DF1*CN1)/3.D0
    ! FNC
    DF2=( (F(2)*FB(2)-F(1)*FB(1))/(S(2)-S(1)) +(F(3)*FB(3)-F(1)*FB(1))/(S(3)-S(1)) &
      +(F(4)*FB(4)-F(1)*FB(1))/(S(4)-S(1)) )*(EE-S(1))
    DF1=(F(2)*FB(2)-F(1)*FB(1))*( (EE-S(1))/(S(2)-S(1))+(S(3)-EE)/(S(3)-S(2)) &
      +(S(4)-EE)/(S(4)-S(2)) ) &
      +( (F(3)*FB(3)-F(1)*FB(1))/(S(3)-S(2))+(F(4)*FB(4)-F(1)*FB(1))/ &
      (S(4)-S(2)) )*(EE-S(2))
    FNC=F(1)*FB(1)*CN+(DF2*CN2-DF1*CN1)/3.D0
    ! FND
    DF2=( (G(2)*FC(2)-G(1)*FC(1))/(S(2)-S(1)) +(G(3)*FC(3)-G(1)*FC(1))/(S(3)-S(1)) &
      +(G(4)*FC(4)-G(1)*FC(1))/(S(4)-S(1)) )*(EE-S(1))
    DF1=(G(2)*FC(2)-G(1)*FC(1))*( (EE-S(1))/(S(2)-S(1))+(S(3)-EE)/(S(3)-S(2)) &
      +(S(4)-EE)/(S(4)-S(2)) ) &
      +( (G(3)*FC(3)-G(1)*FC(1))/(S(3)-S(2))+(G(4)*FC(4)-G(1)*FC(1))/ &
      (S(4)-S(2)) )*(EE-S(2))
    FND=G(1)*FC(1)*CN+(DF2*CN2-DF1*CN1)/3.D0
    !
    goto 23
! CASE_3 E1<E2< E3<E<E4
 22 if( ABS(S(4)-S(3)) <= 0.0 ) then
      CN=0.0
      SN=0.0
      FN=0.0
      FNB=0.0
      FNC=0.0
      FND=0.0
      goto 23
    end if
    DELTA=S(4)+S(3)-S(2)-S(1)
    AM=(S(4)*S(3)-S(2)*S(1))/DELTA
    D4=S(4)-AM
    V1=V/DELTA
    DEE=EE-S(4)
    FRAC=DEE*DEE/D4/(S(4)-S(3))
    CN=3.D0*V1*FRAC
    SN=V1*(DELTA+DEE*FRAC)
    ! FN
    DF=(F(2)-F(1))*(S(4)-EE)/(S(4)-S(2)) &
        +(F(3)-F(1))*(S(4)-EE)/(S(4)-S(3)) &
        +(F(4)-F(1))*( (EE-S(1))/(S(4)-S(1)) &
        +(EE-S(2))/(S(4)-S(2))+(EE-S(3))/(S(4)-S(3)) )
    FN=F(1)*CN+DF*CN/3.D0
    ! FNB
    DF=(FB(2)-FB(1))*(S(4)-EE)/(S(4)-S(2)) &
        +(FB(3)-FB(1))*(S(4)-EE)/(S(4)-S(3)) &
        +(FB(4)-FB(1))*( (EE-S(1))/(S(4)-S(1)) &
        +(EE-S(2))/(S(4)-S(2))+(EE-S(3))/(S(4)-S(3)) )
    FNB=FB(1)*CN+DF*CN/3.D0
    ! FNC
    DF=(F(2)*FB(2)-F(1)*FB(1))*(S(4)-EE)/(S(4)-S(2)) &
        +(F(3)*FB(3)-F(1)*FB(1))*(S(4)-EE)/(S(4)-S(3)) &
        +(F(4)*FB(4)-F(1)*FB(1))*( (EE-S(1))/(S(4)-S(1)) &
        +(EE-S(2))/(S(4)-S(2))+(EE-S(3))/(S(4)-S(3)) )
    FNC=F(1)*FB(1)*CN+DF*CN/3.D0
    ! FND
    DF=(G(2)*FC(2)-G(1)*FC(1))*(S(4)-EE)/(S(4)-S(2)) &
        +(G(3)*FC(3)-G(1)*FC(1))*(S(4)-EE)/(S(4)-S(3)) &
        +(G(4)*FC(4)-G(1)*FC(1))*( (EE-S(1))/(S(4)-S(1)) &
        +(EE-S(2))/(S(4)-S(2))+(EE-S(3))/(S(4)-S(3)) )
    FND=G(1)*FC(1)*CN+DF*CN/3.D0
    !
    goto 23
! CASE_1 E1<E<E2 <E3<E4
 21 if( ABS(S(2)-S(1)) <= 0.0 ) then
      CN=0.0
      SN=0.0
      FN=0.0
      FNB=0.0
      FNC=0.0
      FND=0.0
      goto 23
    endif
    DELTA=S(4)+S(3)-S(2)-S(1)
    AM=(S(4)*S(3)-S(2)*S(1))/DELTA
    D1=S(1)-AM
    V1=V/DELTA
    DEE=EE-S(1)
    FRAC=DEE*DEE/D1/(S(1)-S(2))
    CN=3.D0*V1*FRAC
    SN=V1*FRAC*DEE
    ! FN
    DF=( (F(2)-F(1))/(S(2)-S(1))+ (F(3)-F(1))/(S(3)-S(1)) &
           +(F(4)-F(1))/(S(4)-S(1)) )*DEE
    FN=F(1)*CN+DF*CN/3.D0
    ! FNB
    DF=( (FB(2)-FB(1))/(S(2)-S(1))+ (FB(3)-FB(1))/(S(3)-S(1)) &
           +(FB(4)-FB(1))/(S(4)-S(1)) )*DEE
    FNB=FB(1)*CN+DF*CN/3.D0
    ! FNC
    DF=( (F(2)*FB(2)-F(1)*FB(1))/(S(2)-S(1))+ (F(3)*FB(3)-F(1)*FB(1))/(S(3)-S(1)) &
           +(F(4)*FB(4)-F(1)*FB(1))/(S(4)-S(1)) )*DEE
    FNC=F(1)*FB(1)*CN+DF*CN/3.D0
    ! FND
    DF=( (G(2)*FC(2)-G(1)*FC(1))/(S(2)-S(1))+ (G(3)*FC(3)-G(1)*FC(1))/(S(3)-S(1)) &
           +(G(4)*FC(4)-G(1)*FC(1))/(S(4)-S(1)) )*DEE
    FND=G(1)*FC(1)*CN+DF*CN/3.D0
    !
    goto 23
!
 20 SN =V
    CN =0.d0
    FN =0.d0
    FNB=0.d0
    FNC=0.d0
    FND=0.d0
!
 23 SNT(LE) =SNT(LE) +SN
    CNT(LE) =CNT(LE) +CN
    FNT(LE) =FNT(LE) +FN
    FNTB(LE)=FNTB(LE)+FNB
    FNTC(LE)=FNTC(LE)+FNC
    FNTD(LE)=FNTD(LE)+FND
 25 CONTINUE
  return
end subroutine PDNS
!----------------------------------------------------------------------
subroutine SORT(S,F,G,O)
  implicit none
  integer(4):: I,J,NP
  real(8),intent(inout):: S(4),F(4),G(4),O(4)
  real(8):: ST, FT, GT, OT
  DO J=1,3
    NP=4-J
    DO I=1,NP
      if(S(I+1)<S(I)) then
        ST=S(I)
        FT=F(I)
        GT=G(I)
        OT=O(I)
        !
        S(I)=S(I+1)
        F(I)=F(I+1)
        G(I)=G(I+1)
        O(I)=O(I+1)
        !
        S(I+1)=ST
        F(I+1)=FT
        G(I+1)=GT
        O(I+1)=OT
      end if
    end do
  end do
  return
end subroutine SORT
!----------------------------------------------------------------------

program calcthhr
!----------------------------------------------------------------
! HISTORY: version
! new avga.f for As3Ca (No.2)                           2013.6.6
! PROJECTED (DECOMPOSED) STATE-DENSITY PROGRAMME
!                                04.3.12, 09.8.30, 11.2.22, 3.15
!----------------------------------------------------------------
! calcthhr
! module constants
! module allocate_data
!----------------------------------------------------------------
  use constants
  use common
  use allocate_data
  implicit none
  integer(4)::I,J,L,II,LL,J1,J2,I2,LE,LB,IP
  ! user setting variables
  integer(4)::NB1,LMAX
  ! C2 for L
  integer(4)::ix
  !
  ! integer(4)::NPY,NPZ,NPO,NPXD,NPYD,NPZD
  integer(4)::NPO,NPXD,NPYD,NPZD
  real(8):: BB1,BB2,BB3,BB4,BB5,BB6,BB7,BB8,BB9,BB10 &
    ,BB11,BB12,BB13,BB14,BB15,BB16,BB17,BB18,BB19,BB20 &
    ,BB21,BB22,BB23,BB24,BB25,BB26,BB27,BB28,BB29,BB30
  real(8):: CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9,CC10 &
    ,CC11,CC12,CC13,CC14,CC15,CC16,CC17,CC18,CC19,CC20 &
    ,CC21,CC22,CC23,CC24,CC25,CC26,CC27,CC28,CC29,CC30
  !
  real(8)::ES1,AAA, ENE, EOA
  !real(8)::TEN,TF1,TFF
  real(8)::AK2,AF2,AE2
  real(8)::EE
  ! output data
  real(8)::WR,CC,C2,C2B,C2C,C2D,TN
  real(8)::CF,CCC,CCC1,CFC,CFC1,CFD,SG2,SGP,SSS
  real(8)::CP,CCP,CYY,CFP,CFP1,CCP1,CXX,CPD,SP2,SPP
  !
  write(6,'(A)') '===============calcthhr v.1.00================'
  !
  write(6,'(A)') 'call read_input'
  call read_input
  !
  !bag fix
  TN = 0.0 
  !
  EB=0.0
  EF=4.0
  NE=8001
  NEMAX=NE
  DE=(EF-EB)/(NE-1)
  ! read(5,*) '2K^2(EF)='GEE
  GEE=1 ! auto
  !
  MML=NATO
  ML=NATO
  !ML=IATO
  !
  NB1=1
  LMAX = READ_LMAX
  !MMIN=1000
  !
  NPX=NPX+1
  NPY=NPY+1
  NPZ=NPZ+1
  NPXD=NPX;NPYD=NPY;NPZD=NPZ
  NPO=NTK
  !
  COK=COK**(2.0/3.0)
  !
  write(6,'(A)') 'allocate NE,SNT,CNT,FNT,FNTB,FNTC,FNTD'
  allocate(SNT(NE))
  allocate(CNT(NE))
  allocate(FNT(NE))
  allocate(FNTB(NE))
  allocate(FNTC(NE))
  allocate(FNTD(NE))
  call initial_data1(NE,SNT,CNT,FNT,FNTB,FNTC,FNTD)
  !
  write(6,'(A)') 'allocate IC,N1,N2,N3,N4,N5,N6,N7,IEX,IWW'
  allocate(IC(NPXD,NPYD,NPZD))
  allocate(N1(NPT))
  allocate(N2(NPT))
  allocate(N3(NPT))
  allocate(N4(NPT))
  allocate(N5(NPT))
  allocate(N6(NPT))
  allocate(N7(NPT))
  allocate(IEX(NPT))
  allocate(IWW(NPT))
  call initial_data2(NPXD,NPYD,NPZD,IC,NPT,N1,N2,N3,N4,N5,N6,N7,IEX,IWW)
  write(6,'(A)') '===========MAKE MESH FOR TETRAHEDRON=========='
  write(6,'(A)') 'make mesh for tetrahedron'
  write(6,'(A)') 'call TGEN'
  call TGEN(NPXD,NPYD,NPZD,IC,N1,N2,N3,N4,N5,N6,N7)
  !
  write(6,'(A)') '===========READ K POINT AND SPACEGROUP========'
  open(1,file='THHR_klist.dat')
    do I=1,NPT
      read(1,'(5x,i5,f10.5)') IEX(I),IWW(I)
    end do
  close(1)
  !
  I=1;J=0
  !
  write(6,'(A)') 'allocate MAXII,IAB,ES'
  allocate(IAB(MAXII))
  allocate(ES(MAXII))
  allocate(AMA1(NPO,MMIN))
  allocate(AKKK(30,NPO,MMIN))
  allocate(AFFF(30,NPO,MMIN))
  call initial_data3(MAXII,IAB,ES)
  write(6,'(A)') '=============READ K^2 AND C^2 DATA============'
  open(12,file='input_K2.dat')
  open(13,file='input_C2.dat')
  do II=1,MAXII
    READ(12,'(i4,31f13.8)') IAB(II),ES(II) &
     ,BB1 ,BB2 ,BB3 ,BB4 ,BB5 ,BB6 ,BB7 ,BB8 ,BB9 ,BB10 &
     ,BB11,BB12,BB13,BB14,BB15,BB16,BB17,BB18,BB19,BB20 &
     ,BB21,BB22,BB23,BB24,BB25,BB26,BB27,BB28,BB29,BB30
    READ(13,'(4x,31f13.8)')         ES1 &
     ,CC1 ,CC2 ,CC3 ,CC4 ,CC5 ,CC6 ,CC7 ,CC8 ,CC9 ,CC10 &
     ,CC11,CC12,CC13,CC14,CC15,CC16,CC17,CC18,CC19,CC20 &
     ,CC21,CC22,CC23,CC24,CC25,CC26,CC27,CC28,CC29,CC30
    if( II==1 ) then
      EB = ES(II) / EV
    !else if( iabs( IAB(II-1)-IAB(II) ) > 5 ) then
    else if( iabs( IAB(II-1)-IAB(II) ) > 1 ) then
      I=I+1
      if( (I>=2) .and. (J<MMIN) ) MMIN=J
      J=0
    end if
    !write(6,*) II,MMIN,IAB(II),ES(II)
    J=J+1
    ! L=#Cmax, I=kpoint, J=#band(related with Energy)
    AMA1(I,J)   =ES(II)
    AKKK(1,I,J) =BB1 *COK
    AKKK(2,I,J) =BB2 *COK
    AKKK(3,I,J) =BB3 *COK
    AKKK(4,I,J) =BB4 *COK
    AKKK(5,I,J) =BB5 *COK
    AKKK(6,I,J) =BB6 *COK
    AKKK(7,I,J) =BB7 *COK
    AKKK(8,I,J) =BB8 *COK
    AKKK(9,I,J) =BB9 *COK
    AKKK(10,I,J)=BB10*COK
    AKKK(11,I,J)=BB11*COK
    AKKK(12,I,J)=BB12*COK
    AKKK(13,I,J)=BB13*COK
    AKKK(14,I,J)=BB14*COK
    AKKK(15,I,J)=BB15*COK
    AKKK(16,I,J)=BB16*COK
    AKKK(17,I,J)=BB17*COK
    AKKK(18,I,J)=BB18*COK
    AKKK(19,I,J)=BB19*COK
    AKKK(20,I,J)=BB20*COK
    AKKK(21,I,J)=BB21*COK
    AKKK(22,I,J)=BB22*COK
    AKKK(23,I,J)=BB23*COK
    AKKK(24,I,J)=BB24*COK
    AKKK(25,I,J)=BB25*COK
    AKKK(26,I,J)=BB26*COK
    AKKK(27,I,J)=BB27*COK
    AKKK(28,I,J)=BB28*COK
    AKKK(29,I,J)=BB29*COK
    AKKK(30,I,J)=BB30*COK
    !
    AFFF(1,I,J) =CC1
    AFFF(2,I,J) =CC2
    AFFF(3,I,J) =CC3
    AFFF(4,I,J) =CC4
    AFFF(5,I,J) =CC5
    AFFF(6,I,J) =CC6
    AFFF(7,I,J) =CC7
    AFFF(8,I,J) =CC8
    AFFF(9,I,J) =CC9
    AFFF(10,I,J)=CC10
    AFFF(11,I,J)=CC11
    AFFF(12,I,J)=CC12
    AFFF(13,I,J)=CC13
    AFFF(14,I,J)=CC14
    AFFF(15,I,J)=CC15
    AFFF(16,I,J)=CC16
    AFFF(17,I,J)=CC17
    AFFF(18,I,J)=CC18
    AFFF(19,I,J)=CC19
    AFFF(20,I,J)=CC20
    AFFF(21,I,J)=CC21
    AFFF(22,I,J)=CC22
    AFFF(23,I,J)=CC23
    AFFF(24,I,J)=CC24
    AFFF(25,I,J)=CC25
    AFFF(26,I,J)=CC26
    AFFF(27,I,J)=CC27
    AFFF(28,I,J)=CC28
    AFFF(29,I,J)=CC29
    AFFF(30,I,J)=CC30
  end do
  close(12)
  close(13)
  !
  AAA=EB*EV/DE-int(EB*EV/DE)
  EB=EB-(AAA+1.d0)*DE/EV
  !
  write(6,'(A)') 'allocate BKKK,BFFF,AMA2'
  allocate(BKKK(NPO,MMIN))
  allocate(BFFF(NPO,MMIN))
  allocate(AMA2(NPO,MMIN))
  do I=1,NPO
    do J=1,MMIN    
      BKKK(I,J) = AKKK(1,I,J)
      AMA2(I,J) = AMA1(I,J) 
      ix = 0
      do L=1,LMAX
        if(AFFF(L,I,J)>0) then
          BFFF(I,J) = BFFF(I,J) + AFFF(L,I,J)
        end if
      end do
      BFFF(I,J) = BFFF(I,J)
    end do
  end do
  !
  do I=1,NPO
    do J1=1,MMIN
      AK2=BKKK(I,J1)
      AF2=BFFF(I,J1)
      AE2=AMA2(I,J1)
      do J2=J1+1,MMIN
        if(AMA2(I,J2) < AE2) then
          BKKK(I,J1) =BKKK(I,J2)
          BFFF(I,J1) =BFFF(I,J2)
          AMA2(I,J1) =AMA2(I,J2)
          BKKK(I,J2)=AK2
          BFFF(I,J2)=AF2
          AMA2(I,J2)=AE2
        endif
      end do
    end do
  end do
  !
  write(6,'(A)') '===============WRITE AKK2.DATA================'
  write(6,'(A)') 'allocate PR,WV'
  allocate(PR(LMAX))
  allocate(WV(LMAX))
  open(4,file='AKK2.DATA')
  open(14,file='AKK3.DATA')
  open(24,file='AKK4.DATA')
  write(4,*) 'EC2005, C2005'
  write(14,*) 'EC2020, C2020'
  write(24,*) 'EC2050, C2050'
  do I=1,NPO
    do LL=1,1
      do J=1,MMIN
        PR(LL)=BFFF(I,J)
        WV(LL)=BKKK(I,J)
        EOA=(WV(LL)**1.5)*(PI/3/MML)
        ENE=AMA2(I,J)
        !
        if( (PR(LL) >= 0.05) ) then
          WRITE(4,'(f13.8,",",f13.8)') ENE,WV(LL)
        endif
        !
        if( (PR(LL) >= 0.20) ) then
          WRITE(14,'(f13.8,",",f13.8)') ENE,WV(LL)
        endif
        !
        if( (PR(LL) >= 0.50) ) then
          WRITE(24,'(f13.8,",",f13.8)') ENE,WV(LL)
        endif
        !
      end do
    end do
  end do
  close(4)
  close(14)
  close(24)
  !
  write(6,'(A)') '=========CALC TH-HR(FDENS+DENS+SORT)=========='
  write(6,'(A)') 'allocate AMA,AKK,AKKB'
  allocate(AMA(NPT,MMIN))
  allocate(AKK(NPT,MMIN))
  allocate(AKKB(NPT,MMIN))
  do I=1,NPT
    I2=IEX(I)
    do J=1,MMIN
      AMA(I,J) =AMA2(I2,J)
      AKK(I,J) =BFFF(I2,J)
      AKKB(I,J)=BKKK(I2,J)
    end do
  end do
  !
  write(6,'(A)') 'allocate E,FK,FKB,FKC,EM'
  allocate(E(NPT))
  allocate(FK(NPT))
  allocate(FKB(NPT))
  allocate(FKC(NPT))
  allocate(EM(NE))
  do LB=NB1,MMIN
    do IP=1,NPT
      E(IP)=AMA(IP,LB)/EV
      FK(IP)=AKK(IP,LB)
      FKB(IP)=AKKB(IP,LB)
      FKC(IP)=FKB(IP)**(2.d0)
    end do
    do LE=1,NE
      EM(LE)=EB+(LE-1)*DE
    end do
    !
    write(6,*) 'Calculation FDENS(band): ',LB, '/ ',MMIN
    CALL FDENS(NPXD,NPYD,NPZD,IC,NPT,N1,N2,N3,N4,N5,N6,N7)
  end do
  !
  write(6,'(A)') '===============WRITE AKK.DATA================='
  open(10,file='AKK.DATA')
  write(10,'(A)') 'ENERGY,X2K_2,EOA,CKmax2,F_E,DOS,NE,VEC'
  do LE=1,NE
    EE=EM(LE)*EV
    if( (EE < 0.0) .and. (LE >= 2 ) ) then
      WR=-EM(LE-1)/DE
      CC=CNT(LE-1)+WR*( CNT(LE)-CNT(LE-1) )
      C2B=FNTB(LE-1)+WR*( FNTB(LE)-FNTB(LE-1) )
      if(CC /= 0 )then
        CCP=C2B/CC
      end if
      if(CCP > 0.0)then
        GEE = CCP
      end if
    end if
  end do
  do LE=1,NE
    EE=EM(LE)*EV
    if( (LE>100) .and. (EM(LE)*(EM(LE)-DE) < 0.0) ) then
      WR=-EM(LE-1)/DE
      CC=CNT(LE-1)+WR*( CNT(LE)-CNT(LE-1) )
      C2B=FNTB(LE-1)+WR*( FNTB(LE)-FNTB(LE-1) )
      CCP=C2B/CC
      if(CCP > 0.0)then
        GEE = CCP
      end if
    end if
  end do
  !
  do LE=1,NE
    EE=EM(LE)*EV
    if( (LE>100) .and. (EM(LE)*(EM(LE)-DE) < 0.0) ) then
      WR=-EM(LE-1)/DE
      CC=CNT(LE-1)+WR*( CNT(LE)-CNT(LE-1) )
      C2=FNT(LE-1)+WR*( FNT(LE)-FNT(LE-1) )
      C2B=FNTB(LE-1)+WR*( FNTB(LE)-FNTB(LE-1) )
      C2C=FNTC(LE-1)+WR*( FNTC(LE)-FNTC(LE-1) )
      C2D=FNTD(LE-1)+WR*( FNTD(LE)-FNTD(LE-1) )
      if(CC <= 0.0) then
        CCP=0.0
        CCP1=0.0
        CP=0.0
        SSS=0.0
      else
        CP=C2/CC
        CCP=C2B/CC
        CYY=CCP*CCP
        CFP=C2C/C2
        CFP1=(CFP**1.5)*(PI/3.0/ML)
        CCP1=(CCP**1.5)*(PI/3.0/ML)
        CXX=C2D/CC
        CPD=CXX**(2.d0)
        SP2=abs(CCP-CPD)
        SPP=ABS(CXX-CYY)
        SSS=sqrt(SPP)
      end if
! CCP
!     write(10,'(A)') 'ENERGY, K2, EOA, CKmax2, F, DOS, NE'
!     write(10,'(A)') 'ENERGY, KX2DB, EOA, CKMAXDB, FUNC, DOS, TNE'
!     WRITE(10,'(F14.8,6E14.7)') .0,CFP,CFP1,CP,SSS,CNT(LE)/ML/EV,SNT(LE)
      WRITE(10,'(8(f10.5,1a))') .0,",",CCP,",",CCP1,",",CP,",",SSS/GEE,",",CNT(LE)/ML/EV,",",SNT(LE),",",SNT(LE)/ML," "
    end if
    !
    CC=CNT(LE)
    C2=FNT(LE)
    C2B=FNTB(LE)
    C2C=FNTC(LE)
    C2D=FNTD(LE)
    TN=TN+DE*CC
    if(CC <= 0.0) then
      CF=0.0
      CCC=0.0
      CCC1=0.0
      CFC=0.0
      CFC1=0.0
      CFD=0.0
      SG2=0.0
      SSS=0.0
    else
      CF=C2/CC
      CCC=C2B/CC
      CFC=C2C/C2
      CFC1=(CFC**1.5)*(PI/3.0/ML)
      CCC1=(CCC**1.5)*(PI/3.0/ML)
      CXX=C2D/CC
      CYY=CCC**(2.d0)
      CFD=CXX**(2.d0)
      SG2=abs(CCC-CFD)
      SGP=ABS(CXX-CYY)
      SSS=sqrt(SGP)
! CCC
    end if
!   write(10,'(A)') 'ENERGY, K2, EOA, CKmax2, F, DOS, NE'
!   write(10,'(A)') 'ENERGY, KX2DB, EOA, CKMAXDB, FUNC, DOS, TNE'
!   WRITE(10,'(F14.8,6E14.7)') EE,CFC,CFC1,CF,SSS,CNT(LE)/ML/EV,SNT(LE)
    WRITE(10,'(8(f10.5,1a))') EE,",",CCC,",",CCC1,",",CF,",",SSS/GEE,",",CNT(LE)/ML/EV,",",SNT(LE),",",SNT(LE)/ML," "
  end do
  CLOSE(10)
  !
  write(6,'(A)') '==================DEALLOCATE=================='
  write(6,'(A)') 'deallocate arrange'
  deallocate(SNT)
  deallocate(CNT)
  deallocate(FNT)
  deallocate(FNTB)
  deallocate(FNTC)
  deallocate(FNTD)
  write(6,'(A)') 'deallocate IC,N1,N2,N3,N4,N5,N6,N7,IEX,IWW'
  deallocate(IC)
  deallocate(N1)
  deallocate(N2)
  deallocate(N3)
  deallocate(N4)
  deallocate(N5)
  deallocate(N6)
  deallocate(N7)
  deallocate(IEX)
  deallocate(IWW)
  write(6,'(A)') 'deallocate MAXII,IAB,ES'
  deallocate(IAB)
  deallocate(ES)
  deallocate(AMA1)
  deallocate(AKKK)
  deallocate(AFFF)
  write(6,'(A)') 'deallocate BKKK,BFFF,AMA2'
  deallocate(BKKK)
  deallocate(BFFF)
  deallocate(AMA2)
  write(6,'(A)') 'deallocate PR,WV'
  deallocate(PR)
  deallocate(WV)
  write(6,'(A)') 'deallocate AMA,AKK,AKKB'
  deallocate(AMA)
  deallocate(AKK)
  deallocate(AKKB)
  write(6,'(A)') 'deallocate E,FK,FKB,FKC,EM'
  deallocate(E)
  deallocate(FK)
  deallocate(FKB)
  deallocate(FKC)
  deallocate(EM)
  write(6,'(A)') '=============CALCURATE TH-HR END=============='
end program calcthhr
!---------------------------------------------------------------------
