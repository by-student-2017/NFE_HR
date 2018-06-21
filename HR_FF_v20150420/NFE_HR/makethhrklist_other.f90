!----------------------------------------------------------------------
! klist for As3Ca(Str=No2) 2012.6.4 original H. Sato
! 2013.10.07 F90 type: M. Inukai
!----------------------------------------------------------------------
! Compile option for gfortran
! gfortran -o makethhrklist -Wall -pedantic -std=f95 -fbounds-check -O
! -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace
! -g makethhrklist.f90 -free -m64
!----------------------------------------------------------------------
program makethhrklist_other
  implicit none
  integer(4)::I,J,K,L
  integer(4)::X,Y,Z, RX, RY, RZ
  integer(4)::I2
  integer(4)::LEX(50000)
  character(7)::ENDC
  integer(4)::LEND,RTEND
  integer(4)::NP,NP2,NPX,NPY,NPZ
  integer(4)::No(4000), Nh(4000),Nk(4000),Nl(4000), DIV
  real(8)::WTT(4000)
  integer(4)::R, NUM_RT
  integer(4)::RT_MATRIX(48,3,3), SameSym(0:99,0:99,0:99)
  integer(4) READ_NR(48)
  !
  ! bug fix
  NUM_RT = 0
! initialization
  do I=1, 4000
    No(I)=0
    Nh(I)=0
    Nk(I)=0
    Nl(I)=0
    DIV=0
    WTT(I)=0.0
  end do
  !
  do I=1,48
    do J=1,3
      do K=1,3
        RT_MATRIX(I,J,K) = 0
      end do
    end do
    READ_NR(I) = 0
  end do
  !
  do I=1, 99
    do J=1, 99
      do K=1, 99
        SameSym(I,J,K)=0
      end do
    end do
  end do
  !
! read: klist
  open(10,file='f03')
  read(10,'(85x,3i3)') NPX, NPY, NPZ
  L=1
  do while(ENDC /= 'END')
    read(10,'(a7)') ENDC
    L=L+1
  end do
  LEND=L-1
  rewind(10)
  do L=1, LEND
    read(10,'(5i10,1x,f4.1)') No(L), Nh(L), Nk(L), Nl(L), DIV, WTT(L)
    !write(6,'(5i10,1x,f4.1)') No(L), Nh(L), Nk(L), Nl(L), DIV, WTT(L)
    Nh(L) = Nh(L) * NPX/DIV 
    Nk(L) = Nk(L) * NPY/DIV
    Nl(L) = Nl(L) * NPZ/DIV
  end do
  close(10)

! read: outputkgen
  open(1,file='f17')
  ! open(2,file='WIEN2k_klist.dat')
  do while(ENDC /= '  DEPEN')
    read(1,'(a7)') ENDC
  end do
  !
  L=0
  do while(ENDC /= '    G1 ')
    read(1,'(a7)') ENDC
    L=L+1
  end do
  RTEND=L-1
  rewind(1)
  !
  do while(ENDC /= '  DEPEN')
    read(1,'(a7)') ENDC
  end do
  !
  R=1
  do L=1,RTEND/4
    read(1,'(1x,4(3x,19x,i3))') &
      READ_NR(1+4*(L-1)),READ_NR(2+4*(L-1)),READ_NR(3+4*(L-1)),READ_NR(4+4*(L-1))
    !write(6,*) READ_NR(1+4*(L-1)),READ_NR(2+4*(L-1)),READ_NR(3+4*(L-1)),READ_NR(4+4*(L-1))
    read(1,'(1x,4(3x,3i5,7x))') &
      RT_MATRIX(R  ,1,1), RT_MATRIX(R  ,1,2), RT_MATRIX(R  ,1,3), &
      RT_MATRIX(R+1,1,1), RT_MATRIX(R+1,1,2), RT_MATRIX(R+1,1,3), &
      RT_MATRIX(R+2,1,1), RT_MATRIX(R+2,1,2), RT_MATRIX(R+2,1,3), &
      RT_MATRIX(R+3,1,1), RT_MATRIX(R+3,1,2), RT_MATRIX(R+3,1,3)
    read(1,'(1x,4(3x,3i5,7x))') &
      RT_MATRIX(R  ,2,1), RT_MATRIX(R  ,2,2), RT_MATRIX(R  ,2,3), &
      RT_MATRIX(R+1,2,1), RT_MATRIX(R+1,2,2), RT_MATRIX(R+1,2,3), &
      RT_MATRIX(R+2,2,1), RT_MATRIX(R+2,2,2), RT_MATRIX(R+2,2,3), &
      RT_MATRIX(R+3,2,1), RT_MATRIX(R+3,2,2), RT_MATRIX(R+3,2,3)
    read(1,'(1x,4(3x,3i5,7x))') &
      RT_MATRIX(R  ,3,1), RT_MATRIX(R  ,3,2), RT_MATRIX(R  ,3,3), &
      RT_MATRIX(R+1,3,1), RT_MATRIX(R+1,3,2), RT_MATRIX(R+1,3,3), &
      RT_MATRIX(R+2,3,1), RT_MATRIX(R+2,3,2), RT_MATRIX(R+2,3,3), &
      RT_MATRIX(R+3,3,1), RT_MATRIX(R+3,3,2), RT_MATRIX(R+3,3,3)
    R=R+4
  end do
  !
  R=1
  do L=1,RTEND/4
    write(6,'(1x,4(3x,19x,i3))') &
      READ_NR(1+4*(L-1)),READ_NR(2+4*(L-1)),READ_NR(3+4*(L-1)),READ_NR(4+4*(L-1))
    !write(6,*) READ_NR(1+4*(L-1)),READ_NR(2+4*(L-1)),READ_NR(3+4*(L-1)),READ_NR(4+4*(L-1))
    write(6,'(1x,4(3x,3i5,7x))') &
      RT_MATRIX(R  ,1,1), RT_MATRIX(R  ,1,2), RT_MATRIX(R  ,1,3), &
      RT_MATRIX(R+1,1,1), RT_MATRIX(R+1,1,2), RT_MATRIX(R+1,1,3), &
      RT_MATRIX(R+2,1,1), RT_MATRIX(R+2,1,2), RT_MATRIX(R+2,1,3), &
      RT_MATRIX(R+3,1,1), RT_MATRIX(R+3,1,2), RT_MATRIX(R+3,1,3)
    write(6,'(1x,4(3x,3i5,7x))') &
      RT_MATRIX(R  ,2,1), RT_MATRIX(R  ,2,2), RT_MATRIX(R  ,2,3), &
      RT_MATRIX(R+1,2,1), RT_MATRIX(R+1,2,2), RT_MATRIX(R+1,2,3), &
      RT_MATRIX(R+2,2,1), RT_MATRIX(R+2,2,2), RT_MATRIX(R+2,2,3), &
      RT_MATRIX(R+3,2,1), RT_MATRIX(R+3,2,2), RT_MATRIX(R+3,2,3)
    write(6,'(1x,4(3x,3i5,7x))') &
      RT_MATRIX(R  ,3,1), RT_MATRIX(R  ,3,2), RT_MATRIX(R  ,3,3), &
      RT_MATRIX(R+1,3,1), RT_MATRIX(R+1,3,2), RT_MATRIX(R+1,3,3), &
      RT_MATRIX(R+2,3,1), RT_MATRIX(R+2,3,2), RT_MATRIX(R+2,3,3), &
      RT_MATRIX(R+3,3,1), RT_MATRIX(R+3,3,2), RT_MATRIX(R+3,3,3)
    R=R+4
  end do
  !
  read(1,'(a20)') ENDC
  write(6,*) ENDC, LEND, RTEND
  close(1)
  !
  write(6,*) ' NPX=',NPX
  write(6,*) ' NPY=',NPY
  write(6,*) ' NPZ=',NPZ
  write(6,*) ' NPO=',LEND
  !
  do R=1, 48
    if(READ_NR(R) == 0 )then
      NUM_RT = R - 1
      write(6,*) "number of rotation matrix = ", NUM_RT
      exit
    end if
  end do
  !
  do L=1, LEND
    do R=1, NUM_RT
      X = RT_MATRIX(R,1,1)*Nh(L) + RT_MATRIX(R,1,2)*Nk(L) + RT_MATRIX(R,1,3)*Nl(L)
      Y = RT_MATRIX(R,2,1)*Nh(L) + RT_MATRIX(R,2,2)*Nk(L) + RT_MATRIX(R,2,3)*Nl(L)
      Z = RT_MATRIX(R,3,1)*Nh(L) + RT_MATRIX(R,3,2)*Nk(L) + RT_MATRIX(R,3,3)*Nl(L)
      if(X<0) X = NPX + X
      if(Y<0) Y = NPY + Y
      if(Z<0) Z = NPZ + Z
      SameSym(X,Y,Z) = No(L)
    end do
  end do
  
  NP=0
  do I=0,NPX
    do J=0,NPY
      do K=0,NPZ
        NP=NP+1
        RX=I; RY=J; RZ=K
        if(I==NPX) RX=0
        if(J==NPY) RY=0
        if(K==NPZ) RZ=0
        NP2 = SameSym(RX,RY,RZ) 
        !if(NP2==0) write(6,*) ' error', NP
        LEX(NP)=NP2
      end do
    end do
  end do
  !
  NP=0
  open(2,file='THHR_klist.dat')
  do I=0,NPX
    do J=0,NPY
      do K=0,NPZ
        NP=NP+1
        I2=LEX(NP)
        write(2,'(2I5,f10.5,3I5)') NP,I2,WTT(I2),I,J,K
      end do
    end do
  end do
  close(2)
  !
  open(13,file='f13')
  write(13,'(a4, i6)') "NPX=",NPX
  write(13,'(a4, i6)') "NPY=",NPY
  write(13,'(a4, i6)') "NPZ=",NPZ
  write(13,'(a4, i10)') "NPT=",NP
  write(6, *) "NPX=", NPX,",NPY=", NPY,",NPZ=", NPZ, ": NPT=",NP
  close(13)
end program makethhrklist_other
