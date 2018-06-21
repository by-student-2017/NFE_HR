!----------------------------------------------------------------------
! klist for As3Ca(Str=No2) 2012.6.4 original H. Sato
! 2013.06.15 F90 type: M. Inukai
!----------------------------------------------------------------------
! Compile option for gfortran
! gfortran -o makethhrklist -Wall -pedantic -std=f95 -fbounds-check -O
! -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace
! -g makethhrklist.f90 -free -m64
!----------------------------------------------------------------------
program makethhrklist
  implicit none
  integer(4)::I,J,K,L,X,Y,Z
  integer(4)::NL,I1,J1,K1
  integer(4)::I2,J2,K2
  integer(4)::IXX(99999),IYY(99999),IZZ(99999),LEX(99999)
  real(8)::WTT(99999)
  character(7)::ENDC
  integer(4)::LEND,KLEND
  integer(4)::NUMOKP(99999),IMNUMOKP,IOKPNUM
  integer(4)::NP,NP2,NPX,NPY,NPZ
  !
  open(10,file='f03')
  read(10,'(85x,3i3)') NPX, NPY, NPZ
  L=1
  do while(ENDC /= 'END')
    read(10,'(a7)') ENDC
    L=L+1
  end do
  LEND=L-1
  close(10)

! read: outputkgen
  open(1,file='f17')
  ! open(2,file='WIEN2k_klist.dat')
  do while(ENDC /= '  point')
    read(1,'(a7)') ENDC
  end do
  !
  L=0
  do while(ENDC /= '  weigh')
    read(1,'(a7)') ENDC
    L=L+1
  end do
  KLEND=L-1
  rewind(1)
  !
  do while(ENDC /= '  point')
    read(1,'(a7)') ENDC
  end do
  !
  do L=1,KLEND
    read(1,'(i5,5x,3i4,4x,i4,f10.5)') NL,I1,J1,K1,NUMOKP(L),WTT(L)
!    write(6,*) '  ',L,' ',I1,' ',J1,' ',K1
    IXX(L)=I1
    IYY(L)=J1
    IZZ(L)=K1
  end do
  !
  read(1,'(a20)') ENDC
  write(6,*) ENDC, LEND, KLEND
  do L=1,LEND
    read(1,'(2i8)') IMNUMOKP,IOKPNUM
!    write(6,*) INUMOKP,IOKPNUM
    do K=1, KLEND
      if( NUMOKP(K) == IOKPNUM ) then
        NUMOKP(K) = IMNUMOKP
      endif
    end do
  end do
  ! close(2)
  close(1)
  !
  write(6,*) ' NPX=',NPX
  write(6,*) ' NPY=',NPY
  write(6,*) ' NPZ=',NPZ
  write(6,*) ' NPO=',LEND

  NP=0
  do I=0,NPX
    do J=0,NPY
      do K=0,NPZ
        X=I
        Y=J
        Z=K
!        if(X == NPX) X=0
!        if(Y == NPY) Y=0
!        if(Z == NPZ) Z=0
!        XM=NPX-X
!        YM=NPY-Y
!        ZM=NPZ-Z
!        if(XM == NPX) XM=0
!        if(YM == NPY) YM=0
!        if(ZM == NPZ) ZM=0
!        write(6,*) ' ',X,' ',Y,' ',Z
        NP=NP+1
!        write(6,103) NP,X,Y,Z
        NP2=0
        do L=1,KLEND
          I2=abs(IXX(L))
          J2=abs(IYY(L))
          K2=abs(IZZ(L))
          if( (X==I2).and.(Y==J2).and.(Z==K2) ) then
            NP2=L
!            write(6,'(i5,6i3,i5)') L,X,I2,Y,J2,Z,K2,NUMOKP(L)
          end if
!          else if
!          if( (XM==I2).and.(YM==J2).and.(ZM==K2) ) then
!            NP2=L
!            write(6,103) NP,X,Y,Z,NUMOKP(L),I2,J2,K2,2
!          end if
        end do
        if(NP2==0) write(6,*) ' error',NP
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
        write(2,'(2I5,f10.5,3I5)') NP,NUMOKP(I2),WTT(I2),IXX(I2),IYY(I2),IZZ(I2)
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
end program makethhrklist
