c klist for As3Ca(Str=No2) 2012.6.4
c
      implicit real*8 (A-H,O-Z)
      parameter(NPMAX=1367,NPX=30)
      integer*4 I,J,K,L,II,JJ,KK,III,JJJ,KKK,IX,IY,IZ,X,Y,Z
      integer*4 XM,YM,ZM
c     dimension IC(0:NPYM,0:NPYM,0:NPYM),IC2(0:NPYM,0:NPYM,0:NPYM)
      dimension IXX(99999),IYY(99999),IZZ(99999),WTT(99999),LEX(99999)
      character*3 ENDC
      character*10 NUMOKP(99999)
c
      open(1,file='wien.klist')
      open(2,file='ca3c2.dat')
      L=0
    1 CONTINUE
      read(1,102) ENDC
      if(ENDC.eq.'END') then
        goto 2
      endif
        L=L+1
       goto 1
    2  LEND=L
        rewind(1)
      do 10 L=1,LEND
      read(1,101) NUMOKP(L),I1,J1,K1,KB,WTT(L)
c      write(6,*) '  ',L,' ',I1,' ',J1,' ',K1
      IXX(L)=I1
      IYY(L)=J1
      IZZ(L)=K1
       LWW=WTT(L)
      write(2,103) L,LWW
   10 continue
      close(1)
      close(2)
c 100 format(i10,4i5,i3)
 101  format(a10,4i10,f5.2)
 102  format(a3)
c
      write(6,*) ' NPT=',LEND
      NP=0
      NPY=NPX
      NPZ=NPX
c
      DO 33 I=0,NPX
      DO 33 J=0,NPY
      DO 33 K=0,NPZ
        X=I
        Y=J
        Z=K
         if(X.eq.NPX) X=0
         if(Y.eq.NPY) Y=0
         if(Z.eq.NPZ) Z=0    
        XM=NPX-X
        YM=NPY-Y
        ZM=NPZ-Z
         if(XM.eq.NPX) XM=0
         if(YM.eq.NPY) YM=0
         if(ZM.eq.NPZ) ZM=0       
c        write(6,*) ' ',X,' ',Y,' ',Z
       NP=NP+1
c      write(6,103) NP,X,Y,Z
          NP2=0
        DO 21 L=1,LEND
        I2=IXX(L)
        J2=IYY(L)
        K2=IZZ(L)
        IF(X.eq.I2.and.Y.eq.J2.and.Z.eq.K2) then
          NP2=L
          write(6,103) NP,X,Y,Z,L,I2,J2,K2,1
          goto 21
        endif
        IF(XM.eq.I2.and.YM.eq.J2.and.ZM.eq.K2) then
            NP2=L
          write(6,103) NP,X,Y,Z,L,I2,J2,K2,2
          goto 21
        ENDIF
   21   CONTINUE
         if(NP2.eq.0) write(6,*) ' error',NP
         LEX(NP)=NP2
   20 CONTINUE
 33   CONTINUE
c
      NP=0
      open(2,file='ca3c1.dat')
      DO 30 I=0,NPX
      DO 30 J=0,NPY
      DO 30 K=0,NPZ
        NP=NP+1
      I2=LEX(NP)
      write(2,103) NP,I2
   30 CONTINUE
      close(2)
  103 format(4i5,3x,5i5)
      stop
      END
