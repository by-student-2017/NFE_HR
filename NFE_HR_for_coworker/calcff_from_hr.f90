!----------------------------------------------------------------------
! * NEW FF.f for sc-Na27Au27Ga31  2013.8.2,3,4 *
! 2013/08/08 F90 type :M. Inukai
!----------------------------------------------------------------------
! common sentence for Fortran 90
! module constants
! module common_nfe
! module allocate_data_nfe
!----------------------------------------------------------------------
! Compile option for gfortran
! gfortran -o calcff -Wall -pedantic -std=f95 -fbounds-check -O
! -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace
! -g calcff.f90 -free -m64
!----------------------------------------------------------------------
module constants
  implicit none
  real(8),parameter:: PI=3.1415926535897932d0
  real(8),parameter:: EV=13.605698
end module constants
!----------------------------------------------------------------------
module common_ff
  implicit none
  real(8):: V
  real(8):: EB,EF,DE,GEE
  integer(4):: NE,NEMAX
  integer(4):: NPX,NPY,NPZ,NPT,MMIN,MAXII,NTK
  integer(4):: MML,ML
  integer(4):: READ_LMAX ! nfe
  real(8):: COK, fac, TIMES
  character(130)::reador
  character(1)::SG
  !
contains
  subroutine read_input_ff
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
    open(75, file="f75")
      read(75,'(5x,i2)') READ_LMAX
    close(75)
    !
    open(60, file = 'f60') ! ff
      read(60, '( a130 )') reador
      write(6, '( a45 )') reador(1:45)
      read(60, '( a130 )') reador
      write(6, '( a45 )') reador(1:45)
      SG=reador(1:1)
    close (60)
       !
    TIMES = 1.0
    open(55, file = 'f55') ! ff
      read(55, '( 19x,f10.5 )') TIMES
    close (55)
    !
    write(6,'(A)') '===============READ INPUT DATA================'
    write(6,'(A,i10)') 'Number of x axis mesh=',NPX
    write(6,'(A,i10)') 'Number of y axis mesh=',NPY
    write(6,'(A,i10)') 'Number of z axis mesh=',NPZ
    write(6,'(A,i10)') 'Total mesh in tetrahedoron method=', NPT
    write(6,'(A,f10.5)') 'Vabc/Va=', COK
    write(6,'(A,f10.5)') 'V/Vcub =', fac
    write(6,'(A,i10)') 'check data=', MMIN
    write(6,'(A,i10)') 'Number of total K^2 and C^2 data=', MAXII
    write(6,'(A,i10)') 'Number of irreducible total k point=',NTK
    write(6,'(A,i2)') 'LMAX=', READ_LMAX 
    write(6,'(A,f10.5)') '(C^2)*TIMES, TIMES=', TIMES
    return
  end subroutine read_input_ff
end module common_ff
!----------------------------------------------------------------------
module allocate_data_ff
  implicit none
  ! NPT
  integer(4),allocatable:: N1(:),N2(:),N3(:),N4(:),N5(:),N6(:),N7(:),IEX(:)
  real(8),allocatable:: IWW(:)
  integer(4),allocatable:: IXX(:),IYY(:),IZZ(:)
  real(8),allocatable:: RIXX(:),RIYY(:),RIZZ(:)
  ! NPXD, NPYD, NPZD
  integer(4),allocatable:: IC(:,:,:)
  ! MAXII
  integer(4),allocatable:: IAB(:)
  real(8),allocatable:: ES(:)
  ! MAXII,MMIN
  real(8),allocatable:: AKKK(:,:,:),AFFF(:,:,:) ! K^2, C^2
  real(8),allocatable:: AMA1(:,:) ! energy
  real(8),allocatable:: AMA2(:,:),BKKK(:,:),BFFF(:,:) ! energy, K^2, C^2
  integer(4),allocatable:: NMAX(:) ! ff
  real(8),allocatable:: NBKKK(:,:),NBFFF(:,:) ! nfe
  real(8),allocatable:: PR(:),WV(:)
  real(8),allocatable:: AMA(:,:),AKKB(:,:), AKK(:,:) ! energy, K^2, C^2
  ! NPT
  real(8),allocatable:: E(:),FK(:),FKB(:),FKC(:) ! energy, K^2, C^2, (K^2)^2
  ! NE
  real(8),allocatable:: EM(:) ! energy
  real(8),allocatable:: CNT(:),SNT(:),FNT(:),FNTB(:),FNTC(:),FNTD(:)
  ! ff
  real(8),allocatable:: CKKK(:),AMA3(:),AMA4(:) ! ff_BG
end module allocate_data_ff
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

program calcff
!----------------------------------------------------------------------
  use constants
  use common_ff
  use allocate_data_ff
  implicit none
  !
  !integer(4)::I,J,L,II,JJ,LL,J1,J2,I2,LE,LB,IP,N
  integer(4)::I,J,II,JJ,LL,N,NN,M,MM
  ! user setting variables
  integer(4)::NB1,LMAX
  integer(4)::MMAX,NMMAX ! ff
  !integer(4)::FLNFE ! nfe
  !integer(4)::JJJ ! nfe
  !
  ! integer(4)::NPY,NPZ,NPO,NPXD,NPYD,NPZD
  integer(4)::NPO,NPXD,NPYD,NPZD
  !integer(4)::NNPO,N1MMIN,N2MMIN
  real(8):: BB1,BB2,BB3,BB4,BB5,BB6,BB7,BB8,BB9,BB10 &
    ,BB11,BB12,BB13,BB14,BB15,BB16,BB17,BB18,BB19,BB20 &
    ,BB21,BB22,BB23,BB24,BB25,BB26,BB27,BB28,BB29,BB30
  real(8):: CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9,CC10 &
    ,CC11,CC12,CC13,CC14,CC15,CC16,CC17,CC18,CC19,CC20 &
    ,CC21,CC22,CC23,CC24,CC25,CC26,CC27,CC28,CC29,CC30
  !
  !real(8)::ES1,AAA, ENE, EOA
  real(8)::ES1,AAA
  !real(8)::TEN,TF1,TFF,
  !real(8)::AD
  !real(8)::AK2,AF2,AE2
  !real(8)::F2,EE
  !real(8)::AKKKK,DDX,EAV,EB1
  real(8)::DDX,EAV,EB1
  real(8)::FAA,FAV
  !integer(4)::JP,NN
  !real(8)::X1,X2,Y,Y1,Y2
  real(8)::X1,Y,Y1,Y2
  integer(4)::IO,SGF,SGFO,FLC
  integer(4)::KPF(50)
  do I=1, 50
    KPF(I)=0
  end do
  !
  !   READ TEXTS IDENTIFYING STRUCTURE CONSTANTS AND CORRECTION
  !   MATRICES USED IN BAND CALCULATION
  !
  write(6,'(A)') '================calcff v.1.00================='
  !
  write(6,'(A)') 'call read_input_ff'
  call read_input_ff
  !
  EB=0.0
  EF=4.0
  NE=8001
  NEMAX=NE
  DE=(EF-EB)/(NE-1)
  ! read(5,*) '2K^2(EF)='GEE
  GEE=1
  !
  MML=1;ML=1
  !
  NB1=1
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
  write(6,'(A)') 'allocate IEX,IWW,IXX,IYY,IZZ,RIXX,RIYY,RIZZ'
  allocate(IEX(NPT))
  allocate(IWW(NPT))
  allocate(IXX(NPT))
  allocate(IYY(NPT))
  allocate(IZZ(NPT))
  allocate(RIXX(NPT))
  allocate(RIYY(NPT))
  allocate(RIZZ(NPT))
  write(6,'(A)') '===========READ K POINT AND SPACEGROUP========='
  open(1,file='THHR_klist.dat')
    do I=1,NPT
      !read(1,'(5x,i5,f10.5)') IEX(I),IWW(I)
      !read(1,'(2I5,f10.5,3I5)') NP,NUMOKP(I2),WTT(I2),IXX(I2),IYY(I2),IZZ(I2)
      read(1,'(5x,I5,f10.5,3I5)') IEX(I),IWW(I),IXX(I),IYY(I),IZZ(I)
      RIXX(I)=float(IXX(I))/float(NPX-1)
      RIYY(I)=float(IYY(I))/float(NPY-1)
      RIZZ(I)=float(IZZ(I))/float(NPZ-1)
      !write(6,*) RIXX(I), RIYY(I), RIZZ(I)
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
  write(6,'(A)') '=============READ K^2 AND C^2 DATA============='
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
    else if( iabs( IAB(II-1)-IAB(II) ) > 5 ) then
      I=I+1
      if( (I>=2) .and. (J<MMIN) ) MMIN=J
      J=0
    end if
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
  EB1=EB
  AAA=EB*EV/DE-int(EB*EV/DE)
  EB=EB-(AAA+1.d0)*DE/EV
  !
  !write(6,*) AAA,AAA*DE,EB1*EV,EB*EV
  !write(6,*) MMIN
  !
  write(6,'(A)') 'allocate BKKK,BFFF,AMA2,NMAX,CKKK,AMA3,AMA4'
  !allocate(BKKK(NPO,MMIN))
  !allocate(BFFF(NPO,MMIN))
  !allocate(AMA2(NPO,MMIN))
  !allocate(NMAX(NPO))
  allocate(BKKK(MMIN,MMIN*30))
  allocate(BFFF(MMIN,MMIN*30))
  allocate(AMA2(MMIN,MMIN*30))
  allocate(NMAX(MMIN))
  !
  LMAX=READ_LMAX
  !DDX=0.04
  DDX=0.0
  IO=-9999
  SGFO=0
  SGF=0
  !
  MMAX=9999
  allocate(CKKK(MMAX))
  allocate(AMA3(MMAX))
  allocate(AMA4(MMAX))
  do M=1,MMAX
    CKKK(M)=-1.0
    AMA3(M)=0.0
    AMA4(M)=0.0
  end do
  MMAX=1
  !
  open(4,file='AKK9.DATA')
  open(8,file='AKK10.DATA')
  open(3,file='AKK9.txt')
  write(6,*) "plot space group=",SG
  write(6,*) "No kx ky kz kpoint"
  do II=1, NPT
    !
    ! BCC case
    if ( SG == "B" ) then
      ! GAMMA (0.0, 0.0, 0.0) in case.klist, No.1 -> (0.0 0.0 0.0) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.0) .and. KPF(1)==0 ) then
        KPF(1)=1;SGF=1
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," G:(0.0,0.0,0.0)"
        !write(4,*) 'G_E G_C2 G_base'
        write(4,*) 'G_E, G_C2'
        write(8,*) 'G_E_BG, G_C2_BG'
        write(3,*) 'CG_G_E, CG_G_C2'
      end if
      ! DELTA (0.5, 0.0, 0.0) in case.klist, No.2-> (0.125 0.125 0.375) in case.outputkgen
      ! DELTA (0.5, 1.0, 1.0) in case.klist -> (0.125 0.125 0.375) in case.outputkgen
      if( (RIXX(II)==0.125) .and. (RIYY(II)==0.125) .and. (RIZZ(II)==0.375) .and. KPF(2)==0 ) then
        KPF(2)=1;SGF=2
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," DELTA:(0.5,0.0,0.0)"
        !write(4,*) 'DEL_E DEL_C2 DEL_base'
        write(4,*) 'DEL_E, DEL_C2'
        write(8,*) 'DEL_E_BG, DEL_C2_BG'
        write(3,*) 'CG_DEL_E, CG_DEL_C2'
      end if
      ! H (1.0, 0.0, 0.0) in case.klist, No.3 -> (0.5 0.5 0.5) in case.outputkgen
      if( (RIXX(II)==0.5) .and. (RIYY(II)==0.5) .and. (RIZZ(II)==0.5) .and. KPF(3)==0  ) then
        KPF(3)=1;SGF=3
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," H:(1.0,0.0,0.0)"
        !write(4,*) 'H_E H_C2 H_base'
        write(4,*) 'H_E, H_C2'
        write(8,*) 'H_E_BG, H_C2_BG'
        write(3,*) 'CG_H_E, CG_H_C2'
      end if
      ! N (0.5, 0.5, 0) in case.klist, No.4 -> (0.0 0.0 0.5) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.5) .and. KPF(4)==0  ) then
        KPF(4)=1;SGF=4
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," N:(0.5,0.5,0.0)"
        !write(4,*) 'N_E N_C2 N_base'
        write(4,*) 'N_E, N_C2'
        write(8,*) 'N_E_BG, N_C2_BG'
        write(3,*) 'CG_N_E, CG_N_C2'
      end if
      ! SIGMA (0.25, 0.25, 0.0) in case.klist, No.5 -> (0.0, 0.0, 0.25) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.25) .and. KPF(5)==0  ) then
        KPF(5)=1;SGF=5
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," SIGMA:(0.25,0.25,0.0)"
        !write(4,*) 'SIG_E SIG_C2 SIG_base'
        write(4,*) 'SIG_E, SIG_C2'
        write(8,*) 'SIG_E_BG, SIG_C2_BG'
        write(3,*) 'CG_SIG_E, CG_SIG_C2'
      end if
      ! LAMBDA (0.25, 0.25, 0.25) in case.klist, No.6 -> (0.125, 0.125, 0.125) in case.outputkgen
      if( (RIXX(II)==0.125) .and. (RIYY(II)==0.125) .and. (RIZZ(II)==0.125) .and. KPF(6)==0  ) then
        KPF(6)=1;SGF=6
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," LAMBDA:(0.25,0.25,0.25)"
        !write(4,*) 'LAM_E LAM_C2 LAM_base'
        write(4,*) 'LAM_E, LAM_C2'
        write(8,*) 'LAM_E_BG, LAM_C2_BG'
        write(3,*) 'CG_LAM_E, CG_LAM_C2'
      end if
      ! P (0.5, 0.5, 0.5) in case.klist, No.7 -> (0.25, 0.25, 0.25) in case.outputkgen
      if( (RIXX(II)==0.25) .and. (RIYY(II)==0.25) .and. (RIZZ(II)==0.25) .and. KPF(7)==0  ) then
        KPF(7)=1;SGF=7
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," P:(0.5,0.5,0.5)"
        !write(4,*) 'P_E P_C2 P_base'
        write(4,*) 'P_E, P_C2'
        write(8,*) 'P_E_BG, P_C2_BG'
        write(3,*) 'CG_P_E, CG_P_C2'
      end if
    !
    ! FCC case
    else if ( SG == "F" ) then
      ! W (1, 0.5, 0.0) in case.klist, No.1 -> (0.25, 0.5, 0.75) in case.outputkgen
      if( (RIXX(II)==0.25) .and. (RIYY(II)==0.5) .and. (RIZZ(II)==0.75) .and. KPF(1)==0 ) then
        KPF(1)=1;SGF=1
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," W:(1.0,0.5,0.0)"
        !write(4,*) 'W_E W_C2 W_base'
        write(4,*) 'W_E, W_C2'
        write(8,*) 'W_E_BG, W_C2_BG'
        write(3,*) 'CG_W_E, CG_W_C2'
      end if
      ! L (0.5, 0.5, 0.5) in case.klist, No.2 -> (0.5, 0.5, 0.5) in case.outputkgen
      ! L (0.5, 0.5,-0.5) in case.klist -> (0.5, 0.5, 0.5) in case.outputkgen
      if( (RIXX(II)==0.5) .and. (RIYY(II)==0.5) .and. (RIZZ(II)==0.5) .and. KPF(2)==0 ) then
        KPF(2)=1;SGF=2
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," L:(0.5,0.5,0.5)"
        !write(4,*) 'L_E L_C2 L_base'
        write(4,*) 'L_E, L_C2'
        write(8,*) 'L_E_BG, L_C2_BG'
        write(3,*) 'CG_L_E, CG_L_C2'
      end if
      ! LAMBDA (0.25, 0.25, 0.25) No.3
      ! LAMBDA (0.5, 0.25, 0.0) in case.klist -> (0.125, 0.25, 0.375) in case.outputkgen
      if( (RIXX(II)==0.125) .and. (RIYY(II)==0.25) .and. (RIZZ(II)==0.375) .and. KPF(3)==0  ) then
        KPF(3)=1;SGF=3
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," LAMBDA:(0.25,0.25,0.25)"
        !write(4,*) 'LAM_E LAM_C2 LAM_base'
        write(4,*) 'LAM_E, LAM_C2'
        write(8,*) 'LAM_E_BG, LAM_C2_BG'
        write(3,*) 'CG_LAM_E, CG_LAM_C2'
      end if
      ! GAMMA (0, 0, 0) in case.klist, No.4 -> (0.0, 0.0, 0.0) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.0) .and. KPF(4)==0  ) then
        KPF(4)=1;SGF=4
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," G:(0.0,0.0,0.0)"
        !write(4,*) 'G_E G_C2 G_base'
        write(4,*) 'G_E, G_C2'
        write(8,*) 'G_E_BG, G_C2_BG'
        write(3,*) 'CG_G_E, CG_G_C2'
      end if
      ! DELTA (0.5, 0, 0) in case.klist, No.5 -> (0.0, 0.25, 0.25) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.25) .and. (RIZZ(II)==0.25) .and. KPF(5)==0  ) then
        KPF(5)=1;SGF=5
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," DELTA:(0.5,0.0,0.0)"
        !write(4,*) 'DEL_E DEL_C2 DEL_base'
        write(4,*) 'DEL_E, DEL_C2'
        write(8,*) 'DEL_E_BG, DEL_C2_BG'
        write(3,*) 'CG_DEL_E, CG_DEL_C2'
      end if
      ! X (1.0, 0, 0) in case.klist, No.6 -> (0.0, 0.5, 0.5) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.5) .and. (RIZZ(II)==0.5) .and. KPF(6)==0  ) then
        KPF(6)=1;SGF=6
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," X:(1.0,0.0,0.0)"
        !write(4,*) 'X_E X_C2 X_base'
        write(4,*) 'X_E, X_C2'
        write(8,*) 'X_E_BG, X_C2_BG'
        write(3,*) 'CG_X_E, CG_X_C2'
      end if
      ! Z (1.0, 0.25, 0.0) in case.klist, No.7 -> (0.125, 0.5, 0.625) in case.outputkgen
      if( (RIXX(II)==0.125) .and. (RIYY(II)==0.5) .and. (RIZZ(II)==0.625) .and. KPF(7)==0  ) then
        KPF(7)=1;SGF=7
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," Z:(1.0,0.25,0.0)"
        !write(4,*) 'Z_E Z_C2 Z_base'
        write(4,*) 'Z_E, Z_C2'
        write(8,*) 'Z_E_BG, Z_C2_BG'
        write(3,*) 'CG_Z_E, CG_Z_C2'
      end if
      ! K (0.75, 0.75, 0.0) in case.klist, No.8 -> (0.0, 0.125, 0.875) in case.outputkgen
      ! K (1.0, 0.75, -0.75) in case.klist -> (0.0, 0.125, 0.875) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.125) .and. (RIZZ(II)==0.875) .and. KPF(8)==0  ) then
        KPF(8)=1;SGF=8
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," K:(0.75,0.75,0.0)"
        !write(4,*) 'K_E K_C2 K_base'
        write(4,*) 'K_E, K_C2'
        write(8,*) 'K_E_BG, K_C2_BG'
        write(3,*) 'CG_K_E, CG_K_C2'
      end if
    !
    ! HCP case
    else if ( SG == "H" ) then
      ! GAMMA (0.0, 0.0, 0.0) in case.klist, No.2 -> (0.0, 0.0, 0.0) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.0) .and. KPF(1)==0 ) then
        KPF(1)=1;SGF=1
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," G:(0.0,0.0,0.0)"
        !write(4,*) 'G_E G_C2 G_base'
        write(4,*) 'G_E, G_C2'
        write(8,*) 'G_E_BG, G_C2_BG'
        write(3,*) 'CG_G_E, CG_G_C2'
      end if
      ! SIGMA (0.25, 0.0, 0.0) in case.klist, No.2 -> (0.25, 0.0, 0.0) in case.outputkgen
      if( (RIXX(II)==0.25) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.0) .and. KPF(2)==0 ) then
        KPF(2)=1;SGF=2
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," SIGMA:(0.25,0.0,0.0)"
        !write(4,*) 'SIG_E SIG_C2 SIG_base'
        write(4,*) 'SIG_E, SIG_C2'
        write(8,*) 'SIG_E_BG, SIG_C2_BG'
        write(3,*) 'CG_SIG_E, CG_SIG_C2'
      end if
      ! M (0.5, 0.0, 0.0) in case.klist, No.3 -> (0.5, 0.0, 0.0) in case.outputkgen
      if( (RIXX(II)==0.5) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.0) .and. KPF(3)==0  ) then
        KPF(3)=1;SGF=3
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," M:(0.5,0.0,0.0)"
        !write(4,*) 'M_E M_C2 M_base'
        write(4,*) 'M_E, M_C2'
        write(8,*) 'M_E_BG, M_C2_BG'
        write(3,*) 'CG_M_E, CG_M_C2'
      end if
      ! K (1/3, 1/3, 0.0) in case.klist, No.4 -> (1/3, 1/3, 0.0) in case.outputkgen
      if( (abs(RIXX(II))==(1.0/3.0)) .and. (abs(RIYY(II))==(1.0/3.0)) .and. (RIZZ(II)==0.0) .and. KPF(4)==0  ) then
        KPF(4)=1;SGF=4
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," K:(1/3,1/3,0.0)"
        !write(4,*) 'K_E K_C2 K_base'
        write(4,*) 'K_E, K_C2'
        write(8,*) 'K_E_BG, K_C2_BG'
        write(3,*) 'CG_K_E, CG_K_C2'
      end if
      ! LAMBDA (1/6, 1/6, 0.0) in case.klist, No.5 -> (1/6, 1/6, 0.0) in case.outputkgen
      if( (abs(RIXX(II))==(1.0/6.0)) .and. (abs(RIYY(II))==(1.0/6.0)) .and. (RIZZ(II)==0.0) .and. KPF(5)==0  ) then
        KPF(5)=1;SGF=5
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," LAMBDA:(1/6,1/6,0.0)"
        !write(4,*) 'LAM_E LAM_C2 LAM_base'
        write(4,*) 'LAM_E, LAM_C2'
        write(8,*) 'LAM_E_BG, LAM_C2_BG'
        write(3,*) 'CG_LAM_E, CG_LAM_C2'
      end if
      ! DELTA (0.0, 0.0, 0.25) in case.klist, No.6 -> (0.0, 0.0, 0.25) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.25) .and. KPF(6)==0  ) then
        KPF(6)=1;SGF=6
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," DELTA:(0.0,0.0,0.25)"
        !write(4,*) 'DEL_E DEL_C2 DEL_base'
        write(4,*) 'DEL_E, DEL_C2'
        write(8,*) 'DEL_E_BG, DEL_C2_BG'
        write(3,*) 'CG_DEL_E, CG_DEL_C2'
      end if
      ! A (0.0, 0.0, 0.5) in case.klist, No.7 -> (0.0, 0.0, 0.5) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.5) .and. KPF(7)==0  ) then
        KPF(7)=1;SGF=7
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," A:(0.0,0.0,0.5)"
        !write(4,*) 'A_E A_C2 A_base'
        write(4,*) 'A_E, A_C2'
        write(8,*) 'A_E_BG, A_C2_BG'
        write(3,*) 'CG_A_E, CG_A_C2'
      end if
    !
    else
    ! (Simple) Cubic case
      ! R (0.5, 0.5, 0.5) in case.klist, No.1 -> (0.5, 0.5, 0.5) in case.outputkgen
      if( (RIXX(II)==0.5) .and. (RIYY(II)==0.5) .and. (RIZZ(II)==0.5) .and. KPF(1)==0 ) then
        KPF(1)=1;SGF=1
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," R:(0.5,0.5,0.5)"
        !write(4,*) 'R_E R_C2 R_base'
        write(4,*) 'R_E, R_C2'
        write(8,*) 'R_E_BG, R_C2_BG'
        write(3,*) 'CG_R_E, CG_R_C2'
      end if
      ! LAMBDA (11/40, 11/40, 11/40) in case.klist, No.2 -> (11/40, 11/40, 11/40) in case.outputkgen
      if( (abs(RIXX(II))==(11.0/40.0)) .and. (abs(RIYY(II))==(11.0/40.0)) &
          .and. (abs(RIZZ(II))==(11.0/40.0)) .and. KPF(2)==0 ) then
        KPF(2)=1;SGF=2
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," LAMBDA:(11/40,11/40,11/40)"
        !write(4,*) 'LAM_E LAM_C2 LAM_base'
        write(4,*) 'LAM_E, LAM_C2'
        write(8,*) 'LAM_E_BG, LAM_C2_BG'
        write(3,*) 'CG_LAM_E, CG_LAM_C2'
      end if
      ! GAMMA (0.0, 0.0, 0.0) in case.klist, No.3 -> (0.0, 0.0, 0.0) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.0) .and. KPF(3)==0  ) then
        KPF(3)=1;SGF=3
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," G:(0.0,0.0,0.0)"
        !write(4,*) 'G_E G_C2 G_base'
        write(4,*) 'G_E, G_C2'
        write(8,*) 'G_E_BG, G_C2_BG'
        write(3,*) 'CG_G_E, CG_G_C2'
      end if
      ! DELTA (0.25, 0.0, 0.0) in case.klist, No.4 -> (0.0, 0.0, 0.25) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.25) .and. KPF(4)==0  ) then
        KPF(4)=1;SGF=4
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," DELTA:(0.25,0.0,0.0)"
        !write(4,*) 'DEL_E DEL_C2 DEL_base'
        write(4,*) 'DEL_E, DEL_C2'
        write(8,*) 'DEL_E_BG, DEL_C2_BG'
        write(3,*) 'CG_DEL_E, CG_DEL_C2'
      end if
      ! X (0.5, 0.0, 0.0) in case.klist, No.5 -> (0.0, 0.0, 0.5) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.0) .and. (RIZZ(II)==0.5) .and. KPF(5)==0  ) then
        KPF(5)=1;SGF=5
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," X:(0.5,0.0,0.0)"
        !write(4,*) 'X_E X_C2 X_base'
        write(4,*) 'X_E, X_C2'
        write(8,*) 'X_E_BG, X_C2_BG'
        write(3,*) 'CG_X_E, CG_X_C2'
      end if
      ! Z (0.5, 0.25, 0.0) in case.klist, No.6 -> (0.0, 0.25, 0.5) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.25) .and. (RIZZ(II)==0.5) .and. KPF(6)==0  ) then
        KPF(6)=1;SGF=6
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," Z:(0.5,0.25,0.0)"
        !write(4,*) 'Z_E Z_C2 Z_base'
        write(4,*) 'Z_E, Z_C2'
        write(8,*) 'Z_E_BG, Z_C2_BG'
        write(3,*) 'CG_Z_E, CG_Z_C2'
      end if
      ! M (0.5, 0.5, 0.0) in case.klist, No.7 -> (0.0, 0.5, 0.5) in case.outputkgen
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.5) .and. (RIZZ(II)==0.5) .and. KPF(7)==0  ) then
        KPF(7)=1;SGF=7
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," M:(0.5,0.5,0.0)"
        !write(4,*) 'M_E M_C2 M_base'
        write(4,*) 'M_E, M_C2'
        write(8,*) 'M_E_BG, M_C2_BG'
        write(3,*) 'CG_M_E, CG_M_C2'
      end if
      ! SIGMA (0.25, 0.25, 0.0) No.8
      if( (RIXX(II)==0.0) .and. (RIYY(II)==0.25) .and. (RIZZ(II)==0.25) .and. KPF(8)==0  ) then
        KPF(8)=1;SGF=8
        I=IEX(II)
        write(6,'(i7,1x,3(f6.3,1x),A)') II,RIXX(II),RIYY(II),RIZZ(II)," SIGMA:(0.25,0.25,0.0)"
        !write(4,*) 'SIG_E SIG_C2 SIG_base'
        write(4,*) 'SIG_E, SIG_C2'
        write(8,*) 'SIG_E_BG, SIG_C2_BG'
        write(3,*) 'CG_SIG_E, CG_SIG_C2'
      end if
    end if
    !
    if( (IO/=I) .or. (SGFO/=SGF) )then
      if(IO==I)then
        write(6,*) "     (previous k point is equal to now k point)"
      end if
    !
      do J=1,MMIN
        !AKKKK=AKKK(1,I,J)
        NN=0
        do JJ=1,MMIN
          do LL=1,LMAX
            ! AD=abs( AKKK(LL,I,JJ)-AKKK(1,I,J) )
            if( abs( AKKK(LL,I,JJ)-AKKK(1,I,J) ) < 0.0001 ) then
              NN=NN+1
              AMA2(J,NN)=AMA1(I,JJ)
              BKKK(J,NN)=AKKK(LL,I,JJ)
              BFFF(J,NN)=AFFF(LL,I,JJ)
              NMAX(J)=NN
              !write(6,*) NN,AMA2(J,NN),BKKK(J,NN),BFFF(J,NN),NMAX(J)
            end if
          end do
        end do
      end do
      ! write(4,*) 'NaAuGa'
      do J=1,MMIN
        do N=1,NMAX(J)
          X1=AMA2(J,N)-DDX
          !X1=AMA2(J,N)+DDX
          Y1=BKKK(J,N)
          Y2=BFFF(J,N)
          if( Y2 < 0.00001 ) then
          else
            write(4,'(f14.8,",",a4)') X1,"NaN"
            write(4,'(f14.8,",",f14.8)') X1,Y1
            write(4,'(f14.8,",",f14.8)') X1,Y1+Y2*TIMES
          end if
          !
          do M=1,MMAX
            if( CKKK(M)==-1.0 )then
              FLC=0
              do MM=1,MMAX
                if( CKKK(MM)==Y1 )then
                  FLC=1
                end if
              end do
              if(FLC==0)then
                CKKK(M)=Y1
                AMA3(M)=X1
                AMA4(M)=X1
                NMMAX=MMAX+1
              end if
            end if
          end do
          MMAX=NMMAX
          !
          do M=1,MMAX
            if( CKKK(M)==Y1 )then
              if( AMA3(M)>=X1 )then
                AMA3(M)=X1
              end if
              if ( AMA4(M)<=X1 )then
                AMA4(M)=X1
              end if
            end if
          end do
          !
        end do
        !write(4,'(f14.8,2(1x,a4))') X1,"NaN","NaN"
      end do
      do M=1,MMAX-1
        write(8,'(f14.8,",",f14.8)') AMA3(M),CKKK(M)
        write(8,'(f14.8,",",f14.8)') AMA4(M),CKKK(M)
        write(8,'(f14.8,",",a4)') AMA4(M),"NaN"
      end do
      MMAX=9999
      do M=1,MMAX
        CKKK(M)=-1.0
        AMA3(M)=0.0
        AMA4(M)=0.0
      end do
      MMAX=1
      !
      do J=1,MMIN
        EAV=0.0
        FAV=0.0
        do N=1,NMAX(J)
          FAA=BFFF(J,N)
          !if(AMA2(J,N) < -3.0) FAA=0.0
          EAV = EAV + AMA2(J,N) * FAA
          FAV = FAV + BFFF(J,N)
        end do
        EAV = EAV / FAV
        Y=BKKK(J,1)
        !if( (Y>=100) .and. (Y<=120) ) then
          write(3,'(4f14.8,",",i5)') EAV,Y
        !end if
      end do
    !
    end if
    IO=I
    SGFO=SGF
  end do
  close(4)
  close(8)
  close(3)
  !
  write(6,'(A)') '==================DEALLOCATE=================='
  write(6,'(A)') 'deallocate IEX,IWW,IXX,IYY,IZZ,RIXX,RIYY,RIZZ'
  deallocate(IEX)
  deallocate(IWW)
  deallocate(IXX)
  deallocate(IYY)
  deallocate(IZZ)
  deallocate(RIXX)
  deallocate(RIYY)
  deallocate(RIZZ)
  write(6,'(A)') 'deallocate MAXII,IAB,ES'
  deallocate(IAB)
  deallocate(ES)
  deallocate(AMA1)
  deallocate(AKKK)
  deallocate(AFFF)
  write(6,'(A)') 'deallocate BKKK,BFFF,AMA2,NMAX,CKKK,AMA3,AMA4'
  deallocate(BKKK)
  deallocate(BFFF)
  deallocate(AMA2)
  deallocate(NMAX)
  deallocate(CKKK)
  deallocate(AMA3)
  deallocate(AMA4)
  write(6,'(A)') '============CALCURATE NFE-HR END=============='
end program calcff
