      FUNCTION PROBKS(ALAM)
      PARAMETER (EPS1=0.001, EPS2=1.E-8)
      A2=-2.*ALAM**2
      FAC=2.
      PROBKS=0.
      TERMBF=0.
      DO 11 J=1,100
        TERM=FAC*EXP(A2*J**2)
        PROBKS=PROBKS+TERM
        IF(ABS(TERM).LT.EPS1*TERMBF.OR.ABS(TERM).LT.EPS2*PROBKS)RETURN
        FAC=-FAC
        TERMBF=ABS(TERM)
11    CONTINUE
      PROBKS=1.
      RETURN
      END