      SUBROUTINE FTEST(DATA1,N1,DATA2,N2,F,PROB)
      DIMENSION DATA1(N1),DATA2(N2)
      CALL AVEVAR(DATA1,N1,AVE1,VAR1)
      CALL AVEVAR(DATA2,N2,AVE2,VAR2)
      IF(VAR1.GT.VAR2)THEN
        F=VAR1/VAR2
        DF1=N1-1
        DF2=N2-1
      ELSE
        F=VAR2/VAR1
        DF1=N2-1
        DF2=N1-1
      ENDIF
      PROB = BETAI(0.5*DF2,0.5*DF1,DF2/(DF2+DF1*F))
     *    +(1.-BETAI(0.5*DF1,0.5*DF2,DF1/(DF1+DF2/F)))
      RETURN
      END
