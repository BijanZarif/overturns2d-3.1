C        Generated by TAPENADE     (INRIA, Tropics team)
C  Tapenade 3.6 (r4343) - 10 Feb 2012 10:52
C
C  Differentiation of lim_minmod in reverse (adjoint) mode:
C   gradient     of useful results: lim_minmod
C   with respect to varying inputs: x y
C
C*************************************************************************
      SUBROUTINE LIM_MINMOD_BQ(x, xb, y, yb, lim_minmodb)
      IMPLICIT NONE
C
C*************************************************************************
      REAL x, y
      REAL xb, yb
      REAL lim_minmodb
      REAL lim_minmod
      INTRINSIC SIGN
      INTRINSIC ABS
      REAL abs1
      REAL abs0
C
      IF (SIGN(1.0, x) .NE. SIGN(1.0, y)) THEN
        xb = 0.0
        yb = 0.0
      ELSE
        IF (x .GE. 0.) THEN
          abs0 = x
        ELSE
          abs0 = -x
        END IF
        IF (y .GE. 0.) THEN
          abs1 = y
        ELSE
          abs1 = -y
        END IF
        IF (abs0 .LT. abs1) THEN
          xb = lim_minmodb
          yb = 0.0
        ELSE
          yb = lim_minmodb
          xb = 0.0
        END IF
      END IF
      END
