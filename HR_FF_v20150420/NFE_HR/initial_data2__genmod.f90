        !COMPILER-GENERATED INTERFACE MODULE: Wed Mar 12 10:27:54 2014
        MODULE INITIAL_DATA2__genmod
          INTERFACE 
            SUBROUTINE INITIAL_DATA2(NPXD,NPYD,NPZD,IC,NPT,N1,N2,N3,N4, &
     &N5,N6,N7,IEX,IWW)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: NPZD
              INTEGER(KIND=4), INTENT(IN) :: NPYD
              INTEGER(KIND=4), INTENT(IN) :: NPXD
              INTEGER(KIND=4), INTENT(INOUT) :: IC(NPXD,NPYD,NPZD)
              INTEGER(KIND=4), INTENT(INOUT) :: N1(NPT)
              INTEGER(KIND=4), INTENT(INOUT) :: N2(NPT)
              INTEGER(KIND=4), INTENT(INOUT) :: N3(NPT)
              INTEGER(KIND=4), INTENT(INOUT) :: N4(NPT)
              INTEGER(KIND=4), INTENT(INOUT) :: N5(NPT)
              INTEGER(KIND=4), INTENT(INOUT) :: N6(NPT)
              INTEGER(KIND=4), INTENT(INOUT) :: N7(NPT)
              INTEGER(KIND=4), INTENT(INOUT) :: IEX(NPT)
              REAL(KIND=8), INTENT(INOUT) :: IWW(NPT)
            END SUBROUTINE INITIAL_DATA2
          END INTERFACE 
        END MODULE INITIAL_DATA2__genmod
