      SUBROUTINE LOPW(NAT)
!
      use matrices, only: HSROWS, KZZ, XK, YK, ZK
      use lolog, only : nlo, ilo
      use lstapw, only  : NV
      use rotmat, only: ROTIJ, ROTLOC
      use struk, only : POS, MULT, NDF
      use parallel, only: abort_parallel
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Scalar Arguments
!
      INTEGER            NAT
!
!     ..................................................................
!
!        generates the LAPW (K+G)-vector for local orbitals
!
!     ..................................................................
!
!        Locals
!
      INTEGER            IA1, IEQ, IIX, INDEX, J, K, KOFF, L, LM, LMDN
      INTEGER            LMUP, LMX, N, NATX, NATXX, NB, NBM, ipass
      INTEGER            JLO
      DOUBLE PRECISION   HL, RKGM, SX, TPI
      DOUBLE PRECISION   ROTV1(3), ROTV2(3), VEC(3)
      COMPLEX*16         CC
!      COMPLEX*16         HH((LOMAX+1)**2*NDF,(LOMAX+1)**2*NDF)
      COMPLEX*16         HH((2*LOMAX+1)*48,(2*LOMAX+1)*48)
      COMPLEX*16         SF(NDF), YL(0:(LOMAX+1)**2) !,nv:HSROWS)
!
!        External Subroutines
!
      EXTERNAL           ROTATE, YLM
!
!        Intrinsic Functions
!
      INTRINSIC          ATAN, DCMPLX, DCONJG, EXP, SQRT
!
!     ** Maybe Experiment **
      DOUBLE PRECISION  VEC2(3), TMP1, TMP2
!
      TPI = 8.0D+0*ATAN(1.0D+0)
!
      ipass=0
 1    continue
      KOFF = NV
      IA1 = 0
      DO 140 N = 1, NAT
         DO 130 L = 0, LOMAX
!            IF (LOOR(L,N)) THEN
            do jlo=1,ilo(l,n)
               LMDN = L*L + 1
               LMUP = (L+1)*(L+1)
               INDEX = 0
               NB = 0
               NBM = MULT(N)*(1+LMUP-LMDN)
               DO 120 IEQ = 1, MULT(N)
                  DO 110 LM = LMDN, LMUP
                     NB = NB + 1
                     K = KOFF + NB
   10                CONTINUE
                     INDEX = INDEX + 1
                     IF (INDEX .GT. NV) GOTO 900
!                  WRITE (6,*) 'INDEX,K,N,L,IEQ,LM',INDEX,K,N,L,IEQ,LM
                     KZZ(1,K) = KZZ(1,INDEX)
                     KZZ(2,K) = KZZ(2,INDEX)
                     KZZ(3,K) = KZZ(3,INDEX)
                     XK(K) = XK(INDEX)
                     YK(K) = YK(INDEX)
                     ZK(K) = ZK(INDEX)
                     RKGM = SQRT(XK(K)*XK(K)+YK(K)*YK(K)+ZK(K)*ZK(K))
                     IF (NBM .NE. 1) THEN
                        DO 20 NATX = 1, MULT(N)
                           NATXX = IA1 + NATX
                           SX = KZZ(1,K)*POS(1,NATXX) + &
                                KZZ(2,K)*POS(2,NATXX) + &
                                KZZ(3,K)*POS(3,NATXX)
!                          SF(NATX) = EXP(DCMPLX(0.0D+0,TPI*SX))
                           SF(NATX) = DCMPLX(DCOS(TPI*SX),DSIN(TPI*SX))
   20                   CONTINUE
                        IIX = 0
                        DO 50 NATX = 1, MULT(N)
                           IF (RKGM .LE. 1.0D-5) THEN
                              DO 30 LMX = LMDN, LMUP
!                                YL(LMX-1,K) = (0.0D+0,0.0D+0)
                                 YL(LMX-1) = 0.0D0
   30                         CONTINUE
!                             YL(0,K) = (1.0D+0,0.0D+0)
                              YL(0) = 1.D0 
                           ELSE
                              VEC(1) = XK(K)
                              VEC(2) = YK(K)
                              VEC(3) = ZK(K)
                              CALL ROTATE(VEC,ROTIJ(1,1,IA1+NATX),ROTV1)
                              CALL ROTATE(ROTV1,ROTLOC(1,1,N),ROTV2)
                              CALL YLM(ROTV2,LOMAX,YL(0)) !,K))
                           ENDIF
                           DO 40 LMX = LMDN, LMUP
                              IIX = IIX + 1
                              HH(IIX,NB) = SF(NATX)*YL(LMX-1) !,K)
   40                      CONTINUE
   50                   CONTINUE
                        IF (NB .NE. 1) THEN
                           DO 80 J = 1, NB - 1
                              CC = (0.0D+0,0.0D+0)
                              DO 60 IIX = 1, NBM
                                 CC = CC + HH(IIX,NB)*DCONJG(HH(IIX,J))
   60                         CONTINUE
                              DO 70 IIX = 1, NBM
                                 HH(IIX,NB) = HH(IIX,NB) - CC*HH(IIX,J)
   70                         CONTINUE
   80                      CONTINUE
                        ENDIF
                        HL = 0.0D+0
                        DO 90 IIX = 1, NBM
                           HL = HL + DCONJG(HH(IIX,NB))*HH(IIX,NB)
   90                   CONTINUE
!
!       Change here to increase test so now it is RMS > 0.1D0
!
                        IF (HL .LE. dble(NBM)*1.0D-2) then
                          if(ipass.eq.2) then
        WRITE (6,6001) n,l,ieq,index,K, RKGM, (KZZ(J,K),J=1,3),hl
 6001   format(' atom',i3,' L',i2,' equiv',i2,' PW',i5,i5,f10.4,3i4,e15.5)
                          IF (HL .LE. dble(NBM)*5.0D-3)    GOTO 10  !! PB Change here, was 1.0D-3
                          else
                          GOTO 10  !! PB Change here, was 1.0D-3
                          endif
                        endif
!
!       *** Maybe Experiment **
!       Ensure that the K vectors themselves are not linearly dependent
!       This can be done quicker by storing the VEC2 values....
                        if(ipass.eq.0) then
                          IF(NB .GT. 2) THEN
                                TMP1=KZZ(1,K)*KZZ(1,K)+KZZ(2,K)*KZZ(2,K)+KZZ(3,K)*KZZ(3,K)
                                IF(TMP1 .lt. 1D-15)TMP1=1.D0
                                VEC(1:3)=KZZ(1:3,K)/sqrt(TMP1)
                                DO J = 1, NB-1
                                        IIX=KOFF+J
                                        TMP2=KZZ(1,IIX)*KZZ(1,IIX)+KZZ(2,IIX)*KZZ(2,IIX)+KZZ(3,IIX)*KZZ(3,IIX)
                                        IF(TMP2 .lt. 1D-15)TMP2=1.D0
                                        VEC2(1:3)=KZZ(1:3,IIX)/sqrt(TMP2)

                                        if((abs(VEC(1)-VEC2(1)).lt. 1D-10) .and. &
                                          ( abs(VEC(2)-VEC2(2)).lt. 1D-10) .and. &
                                          ( abs(VEC(3)-VEC2(3)).lt. 1D-10) ) GOTO 10

                                        if((abs(VEC(1)+VEC2(1)).lt. 1D-10) .and. &
                                          ( abs(VEC(2)+VEC2(2)).lt. 1D-10) .and. &
                                          ( abs(VEC(3)+VEC2(3)).lt. 1D-10) ) GOTO 10
                                ENDDO
                          ENDIF
                        ENDIF
!       *** End of Experiment **

                     ELSE
                        GOTO 110
                     ENDIF
                     HL = 1.0D+0/SQRT(HL)
                     DO 100 IIX = 1, NBM
                        HH(IIX,NB) = HH(IIX,NB)*HL
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
!               KOFF = KOFF + MULT(N)*(LDIFF(N)+1)**2
               KOFF = KOFF + MULT(N)*(2*L + 1)
            ENDdo
  130    CONTINUE
         IA1 = IA1 + MULT(N)
  140 CONTINUE
!
!      DO 150 K=nv+1, NV+NLO
!        WRITE (6,6000) K, RKGM, (KZZ(J,K),J=1,3), XK(K), YK(K), ZK(K)
! 150  CONTINUE
6000  FORMAT (5X,I5,F15.8,3I5,10X,3F10.5,5X)
!
      RETURN
!
!        Error messages
!
  900 if(ipass.eq.0) then
        print*, 'IPASS set to 1 in lopw.f'
        ipass=1
        goto 1
      elseif(ipass.eq.1) then
        print*, 'IPASS set to 2 in lopw.f'
        ipass=2
        goto 1
      endif
      CALL OUTERR('LOPW','Plane waves exhausted ') 
      call abort_parallel
      STOP 'LOPW - Error'
!
!        End of 'LOPW'
!
      END
