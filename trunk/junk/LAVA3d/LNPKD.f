C	@(#)lsgefa.f	1.1	MINERVA	8/7/89	10:52:51
      SUBROUTINE LSGEFA (A, LDA, N, IPVT, INFO)
      IMPLICIT NONE
      INTEGER LDA, N, IPVT(1), INFO
      DOUBLE PRECISION A(LDA,1)
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          A [R:LDA,N] - MATRIX TO BE FACTORED
C
C          LDA [I:1] - LEADING DIMENSION OF A
C
C          N [I:1] - ORDER OF THE MATRIX  A
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          A - CONTAINS AN UPPER TRIANGULAR MATRIX AND THE
C              MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C              THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C              L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C              TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C          IPVT [I:N] - PIVOT INDICES
C
C          INFO [I:1] - COMPLETION FLAG
C              = 0, NORMAL VALUE
C              = K, IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                   CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                   INDICATE THAT LSGESL OR LSGEDI WILL DIVIDE BY ZERO
C                   IF CALLED.  USE  RCOND  IN LSGECO FOR A RELIABLE
C                   INDICATION OF SINGULARITY.
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER J, K, KP1, L, NM1
      DOUBLE PRECISION TEMP, ZERO
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
      INTEGER VSAMXI
C
C          MINERVA INTEGER FUNCTION - VSAMXI
C          MINERVA SUBROUTINES - VSAXPY, VSXTSA
C
C     ..E..  CONSTANTS  .....
C
      DATA ZERO /0.0D0/
C
C  ...3...  NOTES / REMARKS  .......
C
C     ..A..  USAGE  .....
C
C          LSGEFA IS USUALLY CALLED BY LSGECO, BUT IT CAN BE CALLED
C          DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C          (TIME FOR LSGECO) = (1 + 9/N)*(TIME FOR LSGEFA) .
C
C     ..B..  ALGORITHM  .....
C
C     ..C..  PORTABILITY  .....
C
C     ..D..  REFERENCES  .....
C
C          (1) J.J. DONGARRA, C.B. MOLER, J.R. BUNCH, AND G.W. STEWART,
C          LINPACK USER'S GUIDE, SIAM, PHILADELPHIA, 1979
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - SGEFA - CLEVE MOLER, LINPACK, 08/14/78
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/5/80
C
C ****
      INFO = 0
      IF (N .EQ. 1) GO TO 60
C ---               --------
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C ---
      NM1 = N - 1
      DO 50 K = 1, NM1
         KP1 = K + 1
C    ---
C     A  FIND L = PIVOT INDEX
C    ---
         L = K - 1 + VSAMXI (N-K+1, A(K,K), 1)
         IPVT(K) = L
C    ---
C     B  ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C    ---
         IF (A(L,K) .NE. ZERO) GO TO 10
            INFO = K
            GO TO 50
C           --------
   10    CONTINUE
C    ---
C     C  INTERCHANGE IF NECESSARY
C    ---
         IF (L .EQ. K) GO TO 20
            TEMP   = A(L,K)
            A(L,K) = A(K,K)
            A(K,K) = TEMP
   20    CONTINUE
C    ---
C     D  COMPUTE MULTIPLIERS
C    ---
         CALL VSXTSA (N-K, -1.0D0/A(K,K), A(K+1,K), 1)
C    ---
C     E  ROW ELIMINATION WITH COLUMN INDEXING
C    ---
         DO 40 J = KP1, N
            TEMP = A(L,J)
            IF (L .EQ. K) GO TO 30
               A(L,J) = A(K,J)
               A(K,J) = TEMP
   30       CONTINUE
            CALL VSAXPY (N-K, TEMP, A(K+1,K), 1, A(K+1,J), 1)
   40       CONTINUE
   50    CONTINUE
C ---------------
   60 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. ZERO) INFO = N
      RETURN
      END
C	@(#)lsgesl.f	1.1	MINERVA	8/7/89	10:52:52
      SUBROUTINE LSGESL (A, LDA, N, IPVT, B, JOB)
      IMPLICIT NONE
      INTEGER LDA, N, IPVT(1), JOB
      DOUBLE PRECISION A(LDA,1), B(1)
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       SOLVES THE REAL SYSTEM  A * X = B  OR  TRANS(A) * X = B
C       USING THE FACTORS COMPUTED BY LSGEFA
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          A [R:LDA,N] - FACTORS PRODUCED BY LSGEFA
C
C          LDA [I:1] - LEADING DIMENSION OF A
C
C          N [I:1] - ORDER OF THE MATRIX  A
C
C          IPVT [I:N] - PIVOT VECTOR FROM LSGEFA
C
C          B [R:N] - RIGHT HAND SIDE VECTOR
C
C          JOB [I:1] - JOB PATH PARAMETER
C              = 00, TO SOLVE  A*X = B,
C              = 01, TO SOLVE  L*X = B,
C              = 10, TO SOLVE  TRANS(A)*X = B,
C              = 11, TO SOLVE  TRANS(U)*X = B,  WHERE
C                    TRANS(.)  IS THE TRANSPOSE.
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          B - SOLUTION VECTOR  X
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER K, KB, L, NM1
      INTEGER JOB1,JOB2
      DOUBLE PRECISION TEMP
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
      DOUBLE PRECISION VSDOTY
C
C          MINERVA REAL FUNCTION - VSDOTY
C          MINERVA SUBROUTINE - VSAXPY
C          INTRINSIC FUNCTION - MOD
C
C     ..E..  CONSTANTS  .....
C
C  ...3...  NOTES / REMARKS  .......
C
C     ..A..  USAGE  .....
C
C          (1) ERROR CONDITION:
C
C          A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C          ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C          BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C          SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C          CALLED CORRECTLY AND IF LSGECO HAS SET RCOND .GT. 0.0
C          OR LSGEFA HAS SET INFO .EQ. 0 .
C
C          (2) TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C          WITH  P  COLUMNS:
C
C               CALL LSGECO (A, LDA, N, IPVT, RCOND, Z)
C               IF (RCOND IS TOO SMALL) GO TO ...
C               DO 10 J = 1, P
C                  CALL LSGESL (A, LDA, N, IPVT, C(1,J), 0)
C            10    CONTINUE
C
C     ..B..  ALGORITHM
C
C     ..C..  PORTABILITY  .....
C
C     ..D..  REFERENCES  .....
C
C          (1) J.J. DONGARRA, C.B. MOLER, J.R. BUNCH, AND G.W. STEWART,
C          LINPACK USER'S GUIDE, SIAM, PHILADELPHIA, 1979
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - SGESL - LINPACK - CLEVE MOLER, 08/14/78
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/5/80
C
C ****
      NM1 = N - 1
      JOB1 = JOB / 10
      JOB2 = MOD (JOB, 10)
      IF (JOB1 .EQ. 1) GO TO 50
C ---                  --------
C  1  SOLVE  A * X = B
C ---
         IF (N .EQ. 1) GO TO 30
C    ---               --------
C     2  FIRST SOLVE  L*Y = B
C    ---
         DO 20 K = 1, NM1
            L = IPVT(K)
            TEMP = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = TEMP
   10       CONTINUE
            CALL VSAXPY (N-K, TEMP, A(K+1,K), 1, B(K+1), 1)
   20       CONTINUE
   30    CONTINUE
         IF (JOB2 .EQ. 1) GO TO 80
C    ---                  --------
C     B  NOW SOLVE  U*X = Y
C    ---
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / A(K,K)
            CALL VSAXPY (K-1, -B(K), A(1,K), 1, B(1), 1)
   40       CONTINUE
         GO TO 80
C ---------------
   50 CONTINUE
C ---
C  1' SOLVE  TRANS(A) * X = B
C ---
C    ---
C     A  FIRST SOLVE  TRANS(U)*Y = B
C    ---
         DO 60 K = 1, N
            TEMP = VSDOTY (K-1, A(1,K), 1, B(1), 1)
            B(K) = (B(K) - TEMP) / A(K,K)
   60       CONTINUE
         IF (JOB2 .EQ. 1) GO TO 80
C    ---                  --------
C     B  NOW SOLVE TRANS(L)*X = Y
C    ---
         IF (N .EQ. 1) GO TO 80
C                      --------
         DO 70 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + VSDOTY (N-K, A(K+1,K), 1, B(K+1), 1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               TEMP = B(L)
               B(L) = B(K)
               B(K) = TEMP
   70       CONTINUE
C ------------------
   80 CONTINUE
      RETURN
      END
C	@(#)lsgbfa.f	1.1	MINERVA	8/7/89	10:52:46
      SUBROUTINE LSGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
      IMPLICIT NONE
      INTEGER LDA, N, ML, MU, IPVT(1), INFO
      DOUBLE PRECISION ABD(LDA,1)
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          ABD [R:LDA,N] - CONTAINS THE MATRIX IN BAND STORAGE.  THE
C              COLUMNS OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD
C              AND THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C              ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C              SEE THE COMMENTS BELOW FOR DETAILS.
C
C          LDA [I:1] - LEADING DIMENSION OF ABD
C              LDA MUST BE .GE. 2*ML + MU + 1 .
C
C          N [I:1] - ORDER OF THE MATRIX  A
C
C          ML [I:1] - NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL
C              0 .LE. ML .LT. N .
C
C          MU [I:1] - NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL
C              0 .LE. MU .LT. N .
C              MORE EFFICIENT IF  ML .LE. MU .
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          ABD - CONTAINS AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C              THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C              THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C              L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C              TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C          IPVT [I:N] - PIVOT INDICES
C
C          INFO [I:1] - COMPLETION FLAG
C              = 0, NORMAL VALUE.
C              = K, IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                   CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                   INDICATE THAT LSGBSL WILL DIVIDE BY ZERO IF
C                   CALLED.  USE  RCOND  IN LSGBCO FOR A RELIABLE
C                   INDICATION OF SINGULARITY.
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER J, JU, JZ, J0, J1, K, KP1, L, LM, M, MM, NM1
      DOUBLE PRECISION TEMP, ZERO
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
      INTEGER VSAMXI
C
C          MINERVA INTEGER FUNCTION - VSAMXI
C          MINERVA SUBROUTINES - VSAXPY, VSFILL, VSXTSA
C          INTRINSIC FUNCTIONS - MAX0, MIN0
C
C     ..E..  CONSTANTS  .....
C
      DATA ZERO /0.0D0/
C
C  ...3...  NOTES / REMARKS  .......
C
C     ..A..  USAGE  .....
C
C          (1) GENERAL:
C
C          LSGBFA IS USUALLY CALLED BY LSGBCO, BUT IT CAN BE CALLED
C          DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C          (2) BAND STORAGE:
C
C          IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C          WILL SET UP THE INPUT.
C
C               ML = (BAND WIDTH BELOW THE DIAGONAL)
C               MU = (BAND WIDTH ABOVE THE DIAGONAL)
C               M = ML + MU + 1
C               DO 20 J = 1, N
C                  I1 = MAX0(1, J-MU)
C                  I2 = MIN0(N, J+ML)
C                  DO 10 I = I1, I2
C                     K = I - J + M
C                     ABD(K,J) = A(I,J)
C            10       CONTINUE
C            20    CONTINUE
C
C          THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C          IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C          ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C          THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C          THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C          ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     ..B..  ALGORITHM  .....
C
C     ..C..  PORTABILITY  .....
C
C          ANSI 66 - NONSTANDARD SUBSCRIPTS
C
C     ..D..  REFERENCES  .....
C
C          (1) J.J. DONGARRA, C.B. MOLER, J.R. BUNCH, AND G.W. STEWART,
C          LINPACK USER'S GUIDE, SIAM, PHILADELPHIA, 1979
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - SGBFA - LINPACK, CLEVE MOLER, 08/14/78
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/6/80
C
C ****
      INFO = 0
      M = ML + MU + 1
      IF (N .EQ. 1) GO TO 80
C ---
C  1  ZERO INITIAL FILL-IN COLUMNS
C ---
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 20
         DO 10 JZ = J0, J1
            CALL VSFILL (JZ+1-J0, ZERO, ABD(ML+J0-JZ,JZ), 1)
   10       CONTINUE
   20 CONTINUE
C ---
C  2  GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C ---
      NM1 = N - 1
      JZ = J1
      JU = 0
      DO 70 K = 1, NM1
C    ---
C     A  ZERO NEXT FILL-IN COLUMN
C    ---
         JZ = JZ + 1
         IF (JZ .LE. N .AND. ML .GE. 1) CALL VSFILL(ML,ZERO,ABD(1,JZ),1)
C    ---
C     B  FIND L = PIVOT INDEX
C    ---
         LM = MIN0 (ML, N-K)
         L = M - 1 + VSAMXI (LM+1, ABD(M,K), 1)
         IPVT(K) = L + K - M
C    ---
C     C  ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C    ---
         IF (ABD(L,K) .NE. ZERO) GO TO 30
            INFO = K
            GO TO 70
C           --------
   30    CONTINUE
C    ---
C     D  INTERCHANGE IF NECESSARY
C    ---
         IF (L .EQ. M) GO TO 40
            TEMP     = ABD(L,K)
            ABD(L,K) = ABD(M,K)
            ABD(M,K) = TEMP
   40    CONTINUE
C    ---
C     E  COMPUTE MULTIPLIERS
C    ---
         CALL VSXTSA (LM, -1.0D0/ABD(M,K), ABD(M+1,K), 1)
C    ---
C     F  ROW ELIMINATION WITH COLUMN INDEXING
C    ---
         KP1 = K + 1
         JU = MIN0 (MAX0(JU,MU+IPVT(K)), N)
         IF (JU .LT. KP1) GO TO 70
C                         --------
            MM = M
            DO 60 J = KP1, JU
               L = L - 1
               MM = MM - 1
               TEMP = ABD(L,J)
               IF (L .EQ. MM) GO TO 50
                  ABD(L,J)  = ABD(MM,J)
                  ABD(MM,J) = TEMP
   50          CONTINUE
               CALL VSAXPY(LM, TEMP, ABD(M+1,K),1, ABD(MM+1,J),1)
   60          CONTINUE
   70    CONTINUE
C ---------------
   80 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. ZERO) INFO = N
      RETURN
      END
C	@(#)lsgbsl.f	1.1	MINERVA	8/7/89	10:52:47
      SUBROUTINE LSGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      IMPLICIT NONE
      INTEGER LDA, N, ML, MU, IPVT(1), JOB
      INTEGER JOB1,JOB2
      DOUBLE PRECISION ABD(LDA,1), B(1)
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       SOLVES THE REAL BAND SYSTEM  A * X = B  OR  TRANS(A) * X = B
C       USING THE FACTORS COMPUTED BY LSGBFA
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          ABD [R:LDA,N] - FACTORS PRODUCED BY LSGBFA
C
C          LDA [I:1] - LEADING DIMENSION OF ABD
C
C          N [I:1] - ORDER OF THE MATRIX  A
C
C          ML [I:1] - NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL
C
C          MU [I:1] - NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL
C
C          IPVT [I:N] - PIVOT VECTOR FROM LSGBFA
C
C          B [R:N] - RIGHT HAND SIDE VECTOR
C
C          JOB [I:1] - JOB PATH PARAMETER
C              = 00, TO SOLVE  A*X = B,
C              = 01, TO SOLVE  L*X = B,
C              = 10, TO SOLVE  TRANS(A)*X = B,
C              = 11, TO SOLVE  TRANS(U)*X = B, WHERE
C                    TRANS(.)  IS THE TRANSPOSE.
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          B - SOLUTION VECTOR  X
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER K, KB, L, LM, M, NM1
      DOUBLE PRECISION TEMP
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
      DOUBLE PRECISION VSDOTY
C
C          MINERVA REAL FUNCTION -  VSDOTY
C          MINERVA SUBROUTINE - VSAXPY
C          INTRINSIC FUNCTION - MIN0, MOD
C
C     ..E..  CONSTANTS  .....
C
C  ...3...  NOTES / REMARKS  .......
C
C     ..A..  USAGE  .....
C
C          (1) ERROR CONDITION:
C
C          A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C          ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C          BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C          SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C          CALLED CORRECTLY AND IF LSGBCO HAS SET RCOND .GT. 0.0
C          OR LSGBFA HAS SET INFO .EQ. 0 .
C
C          (2) TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C          WITH  P  COLUMNS:
C
C               CALL LSGBCO (ABD, LDA, N, ML, MU, IPVT, RCOND, Z)
C               IF (RCOND IS TOO SMALL) GO TO ...
C               DO 10 J = 1, P
C                  CALL LSGBSL (ABD, LDA, N, ML, MU, IPVT, C(1,J), 0)
C            10    CONTINUE
C
C     ..B..  ALGORITHM  .....
C
C     ..C..  PORTABILITY  .....
C
C          ANSI 66 - NONSTANDARD SUBSCRIPTS
C
C     ..D..  REFERENCES  .....
C
C          (1) J.J. DONGARRA, C.B. MOLER, J.R. BUNCH, AND G.W. STEWART,
C          LINPACK USER'S GUIDE, SIAM, PHILADELPHIA, 1979
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - SGBSL - LINPACK, CLEVE MOLER, 08/14/78
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/6/80
C
C ****
      M = MU + ML + 1
      NM1 = N - 1
      JOB1 = JOB / 10
      JOB2 = MOD (JOB, 10)
      IF (JOB1 .EQ. 1) GO TO 50
C ---                  --------
C  1  SOLVE  A * X = B
C ---
         IF (ML .EQ. 0) GO TO 30
         IF (N  .EQ. 1) GO TO 30
C    ---                --------
C     A  FIRST SOLVE L*Y = B
C    ---
         DO 20 K = 1, NM1
            LM = MIN0 (ML, N-K)
            L = IPVT(K)
            TEMP = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = TEMP
   10       CONTINUE
            CALL VSAXPY (LM, TEMP, ABD(M+1,K), 1, B(K+1), 1)
   20       CONTINUE
   30    CONTINUE
         IF (JOB2 .EQ. 1) GO TO 80
C    ---                  --------
C     B  NOW SOLVE  U*X = Y
C    ---
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / ABD(M,K)
            LM = MIN0(K,M) - 1
            CALL VSAXPY (LM, -B(K), ABD(M-LM,K), 1, B(K-LM), 1)
   40       CONTINUE
         GO TO 80
C ---------------
   50 CONTINUE
C ---
C  1' SOLVE  TRANS(A) * X = B
C ---
C    ---
C     A  FIRST SOLVE  TRANS(U)*Y = B
C    ---
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            TEMP = VSDOTY (LM, ABD(M-LM,K), 1, B(K-LM), 1)
            B(K) = (B(K) - TEMP) / ABD(M,K)
   60       CONTINUE
         IF (JOB2 .EQ. 1) GO TO 80
C    ---                  --------
C     B   NOW SOLVE TRANS(L)*X = Y
C    ---
         IF (ML .EQ. 0) GO TO 80
         IF (N  .EQ. 1) GO TO 80
C                       --------
         DO 70 KB = 1, NM1
            K = N - KB
            LM = MIN0 (ML, N-K)
            B(K) = B(K) + VSDOTY (LM, ABD(M+1,K), 1, B(K+1), 1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               TEMP = B(L)
               B(L) = B(K)
               B(K) = TEMP
   70       CONTINUE
C ------------
   80 CONTINUE
      RETURN
      END
C	@(#)vsaxpy.f	1.1	MINERVA	8/7/89	10:54:13
      SUBROUTINE VSAXPY (N, SA, SX, INCX, SY, INCY)
      IMPLICIT NONE
      INTEGER INCX, INCY, N
      DOUBLE PRECISION SA, SX(1), SY(1)
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       REAL CONSTANT TIMES A REAL VECTOR PLUS A REAL VECTOR
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          N [I:1] - LENGTH OF VECTORS
C
C          SA [R:1] - REAL CONSTANT
C
C          SX, SY [R:N] - REAL VECTORS
C
C          INCX, INCY [I:1] - INCREMENTS OF STORAGE
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          SY - RESULTANT VECTOR
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER I, IX, IY, M, MP1
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
C          INTRINSIC FUNCTION - MOD
C
C     ..E..  CONSTANTS  .....
C
C  ...3...  NOTES / REMARKS  ......
C
C     ..A..  USAGE  .....
C
C     ..B..  ALGORITHM  .....
C
C          USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE
C
C     ..C..  PORTABILITY  .....
C
C     ..D..  REFERENCES  .....
C
C          (1) J.J. DONGARRA, C.B. MOLER, J.R. BUNCH, AND G.W. STEWART,
C          LINPACK USERS' GUIDE, SIAM, PHILADELPHIA, 1979
C
C          (2) C. LAWSON, R. HANSON, D. KINCAID, AND F. KROGH,
C          "BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE",
C          ACM TRANS. MATH. SOFTWARE 5(1979), PP308-325
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - SAXPY - JACK DONGARRA, LINPACK, 3/11/78
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/1/80
C
C ****
      IF (N .LE. 0)      GO TO 60
      IF (SA .EQ. 0.0D0) GO TO 60
C                        --------
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C ---
C  1  UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
C ---
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         SY(IY) = SY(IY) + SA*SX(IX)
         IX = IX + INCX
         IY = IY + INCY
   10    CONTINUE
      GO TO 60
C     --------
   20 CONTINUE
C ---
C  1  BOTH INCREMENTS EQUAL TO 1
C ---
      M = MOD (N, 4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1, M
         SY(I) = SY(I) + SA*SX(I)
   30    CONTINUE
      IF (N .LT. 4) GO TO 60
C                   --------
   40 CONTINUE
      MP1 = M + 1
      DO 50 I = MP1, N, 4
         SY(I)   = SY(I)   + SA*SX(I)
         SY(I+1) = SY(I+1) + SA*SX(I+1)
         SY(I+2) = SY(I+2) + SA*SX(I+2)
         SY(I+3) = SY(I+3) + SA*SX(I+3)
   50    CONTINUE
C ---------------
   60 CONTINUE
      RETURN
      END
C	@(#)vsxtsa.f	1.1	MINERVA	8/7/89	10:54:26
      SUBROUTINE VSXTSA (N, SA, SX, INCX)
      IMPLICIT NONE
      INTEGER INCX, N
      DOUBLE PRECISION SA, SX(1)
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       SCALES A REAL VECTOR BY A REAL CONSTANT
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          N [I:1] - LENGTH OF VECTOR
C
C          SA [R:1] - REAL CONSTANT
C
C          SX [R:N] - REAL VECTOR
C
C          INCX [I:1] - INCREMENT OF STORAGE
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          SX - SCALED VECTOR
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER I, IX, M, MP1
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
C          INTRINSIC FUNCTION - MOD
C
C     ..E..  CONSTANTS  .....
C
C  ...3...  NOTES / REMARKS  ......
C
C     ..A..  USAGE  .....
C
C     ..B..  ALGORITHM  .....
C
C          USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1
C
C     ..C..  PORTABILITY  .....
C
C     ..D..  REFERENCES  .....
C
C          (1) J.J. DONGARRA, C.B. MOLER, J.R. BUNCH, AND G.W. STEWART,
C          LINPACK USERS' GUIDE, SIAM, PHILADELPHIA, 1979
C
C          (2) C. LAWSON, R. HANSON, D. KINCAID, AND F. KROGH,
C          "BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE",
C          ACM TRANS. MATH. SOFTWARE 5(1979), PP308-325
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - SSCAL - JACK DONGARRA, LINPACK, 3/11/78
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/1/80
C
C ****
      IF (N .LE. 0) GO TO 60
C                   --------
      IF (INCX .EQ. 1) GO TO 20
C ---
C  1  INCREMENT NOT EQUAL TO 1
C ---
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1, N
         SX(IX) = SA*SX(IX)
         IX = IX + INCX
   10    CONTINUE
      GO TO 60
C     --------
   20 CONTINUE
C ---
C  1  INCREMENT EQUAL TO 1
C ---
      M = MOD (N, 5)
      IF (M .EQ. 0) GO TO 40
C
      DO 30 I = 1, M
         SX(I) = SA*SX(I)
   30    CONTINUE
      IF (N .LT. 5) GO TO 60
C                   --------
   40 CONTINUE
      MP1 = M + 1
      DO 50 I = MP1, N, 5
         SX(I)   = SA*SX(I)
         SX(I+1) = SA*SX(I+1)
         SX(I+2) = SA*SX(I+2)
         SX(I+3) = SA*SX(I+3)
         SX(I+4) = SA*SX(I+4)
   50    CONTINUE
C ---------------
   60 CONTINUE
      RETURN
      END
C	@(#)vsfill.f	1.1	MINERVA	8/7/89	10:54:15
      SUBROUTINE VSFILL (N, SA, SX, INCX)
      IMPLICIT NONE
      INTEGER N, INCX
      DOUBLE PRECISION SX(N), SA
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       FILLS A REAL VECTOR WITH A CONSTANT VALUE
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          N [I:1] - LENGTH OF VECTOR
C
C          SA [R:1] - SCALAR CONSTANT
C
C          INCX [I:1] - INCREMENT OF STORAGE
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          SX [R:N] - VECTOR FILLED WITH CONSTANT
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER I, IX, M, MP1
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
C          INTRINSIC FUNCTION - MOD
C
C     ..E..  CONSTANTS  .....
C
C  ...3...  NOTES / REMARKS  ......
C
C     ..A..  USAGE  .....
C
C     ..B..  ALGORITHM  .....
C
C          USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1
C
C     ..C..  PORTABILITY  .....
C
C     ..D..  REFERENCES  .....
C
C          (1) R.G. GRIMES, D.R. KINCAID, AND D.M. YOUNG,
C          ITPACK 2.0 USER'S GUIDE, CNA-150, U. TEXAS,
C          AUSTIN, TEXAS, 1979
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - VFILL - ITPACK
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/3/80
C          MODIFIED - ADDED INCX TO ARGUMENT LIST; PARALLELS
C             VSXTSA - ESM, 6/3/80
C
C ****
      IF (N .LE. 0) GO TO 60
C                   --------
      IF (INCX .EQ. 1) GO TO 20
C ---
C  1  INCREMENT NOT EQUAL TO 1
C ---
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1, N
         SX(IX) = SA
         IX = IX + INCX
   10    CONTINUE
      GO TO 60
C     --------
   20 CONTINUE
C ---
C  1  INCREMENT EQUAL TO 1
C ---
      M = MOD (N, 9)
      IF (M .EQ. 0) GO TO 40
C
      DO 30 I = 1, M
         SX(I) = SA
   30    CONTINUE
      IF (N .LT. 9) GO TO 60
C                   --------
   40 CONTINUE
      MP1 = M + 1
      DO 50 I = MP1, N, 9
         SX(I)   = SA
         SX(I+1) = SA
         SX(I+2) = SA
         SX(I+3) = SA
         SX(I+4) = SA
         SX(I+5) = SA
         SX(I+6) = SA
         SX(I+7) = SA
         SX(I+8) = SA
   50    CONTINUE
C ---------------
   60 CONTINUE
      RETURN
      END
C	@(#)vsamxi.f	1.1	MINERVA	8/7/89	10:54:11
      INTEGER FUNCTION VSAMXI (N, SX, INCX)
      IMPLICIT NONE
      INTEGER INCX, N
      DOUBLE PRECISION SX(1)
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       FINDS INDEX OF THE FIRST ELEMENT HAVING MAXIMUM ABSOLUTE VALUE
C       OF A REAL VECTOR
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          N [I:1] - LENGTH OF VECTOR
C
C          SX [R:N] - REAL VECTOR
C
C          INCX [I:1] - INCREMENT OF STORAGE
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          VSAMXI [I:1] - INDEX OF LARGEST ELEMENT
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER I, IX, ITEMP
      DOUBLE PRECISION SMAX
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
C          INTRINSIC FUNCTION - ABS
C
C     ..E..  CONSTANTS  .....
C
C  ...3...  NOTES / REMARKS  ......
C
C     ..A..  USAGE  .....
C
C     ..B..  ALGORITHM  .....
C
C     ..C..  PORTABILITY  .....
C
C     ..D..  REFERENCES  .....
C
C          (1) J.J. DONGARRA, C.B. MOLER, J.R. BUNCH, AND G.W. STEWART,
C          LINPACK USERS' GUIDE, SIAM, PHILADELPHIA, 1979
C
C          (2) C. LAWSON, R. HANSON, D. KINCAID, AND F. KROGH,
C          "BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE",
C          ACM TRANS. MATH. SOFTWARE 5(1979), PP308-325
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - ISAMAX - JACK DONGARRA, LINPACK, 3/11/78
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/1/80
C
C ****
      ITEMP = 0
      IF (N .LE. 0) GO TO 50
C                   --------
      ITEMP = 1
      IF (N .EQ. 1) GO TO 50
C                   --------
      IF (INCX .EQ. 1) GO TO 30
C ---
C  1  INCREMENT NOT EQUAL TO 1
C ---
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      SMAX = ABS(SX(IX))
      IX = IX + INCX
      DO 20 I = 2, N
         IF (ABS(SX(IX)) .LE. SMAX) GO TO 10
            ITEMP = I
            SMAX = ABS(SX(IX))
   10    CONTINUE
         IX = IX + INCX
   20    CONTINUE
      GO TO 50
C     --------
   30 CONTINUE
C ---
C  1  INCREMENT EQUAL TO 1
C ---
      SMAX = ABS(SX(1))
      DO 40 I = 2, N
         IF (ABS(SX(I)) .LE. SMAX) GO TO 40
            ITEMP = I
            SMAX = ABS(SX(I))
   40    CONTINUE
C ---------------
   50 CONTINUE
      VSAMXI = ITEMP
      RETURN
      END
C	@(#)vsdoty.f	1.1	MINERVA	8/7/89	10:54:15
      DOUBLE PRECISION FUNCTION VSDOTY (N, SX, INCX, SY, INCY)
      IMPLICIT NONE
      INTEGER INCX, INCY, N
      DOUBLE PRECISION SX(1), SY(1)
C ****
C
C  ...1...  FUNCTION / PURPOSE  .......
C
C       FORMS THE DOT PRODUCT OF TWO REAL VECTORS
C
C  ...2...  ARGUMENTS / PARAMETERS / CONSTANTS  .......
C
C     ..A..  INPUT ARGUMENTS  .....
C
C          N [I:1] - LENGTH OF VECTORS
C
C          SX, SY [R:N] - REAL VECTORS
C
C          INCX, INCY [I:1] - INCREMENTS OF STORAGE
C
C     ..B..  OUTPUT ARGUMENTS  .....
C
C          VSDOTY - DOT PRODUCT
C
C     ..C..  INTERNAL / LOCAL VARIABLES  .....
C
      INTEGER I, IX, IY, M, MP1
      DOUBLE PRECISION STEMP
C
C     ..D..  EXTERNAL / SUBPROGRAM REFERENCES  .....
C
C          INTRINSIC FUNCTION - MOD
C
C     ..E..  CONSTANTS  .....
C
C  ...3...  NOTES / REMARKS  ......
C
C     ..A..  USAGE  .....
C
C     ..B..  ALGORITHM  .....
C
C          USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE
C
C     ..C..  PORTABILITY  .....
C
C     ..D..  REFERENCES  .....
C
C          (1) J.J. DONGARRA, C.B. MOLER, J.R. BUNCH, AND G.W. STEWART,
C          LINPACK USERS' GUIDE, SIAM, PHILADELPHIA, 1979
C
C          (2) C. LAWSON, R. HANSON, D. KINCAID, AND F. KROGH,
C          "BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE",
C          ACM TRANS. MATH. SOFTWARE 5(1979), PP308-325
C
C     ..E..  HISTORY  .....
C
C          ORIGIN - SDOT - JACK DONGARRA, LINPACK, 3/11/78
C          MODIFIED - HEADER, STRUCTURE - E.S. MARWIL, 6/1/80
C
C ****
      STEMP = 0.0D0
      IF (N .LE. 0) GO TO 60
C                   --------
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C ---
C  1  UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
C ---
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         STEMP = STEMP + SX(IX)*SY(IY)
         IX = IX + INCX
         IY = IY + INCY
   10    CONTINUE
      GO TO 60
C     --------
   20 CONTINUE
C ---
C  1  BOTH INCREMENTS EQUAL TO 1
C ---
      M = MOD (N, 5)
      IF (M .EQ. 0) GO TO 40
C
      DO 30 I = 1, M
         STEMP = STEMP + SX(I)*SY(I)
   30    CONTINUE
      IF (N .LT. 5) GO TO 60
C                   --------
   40 CONTINUE
      MP1 = M + 1
      DO 50 I = MP1, N, 5
         STEMP = STEMP + SX(I)  *SY(I)   + SX(I+1)*SY(I+1)
     *                 + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3)
     *                 + SX(I+4)*SY(I+4)
   50    CONTINUE
C ---------------
   60 CONTINUE
      VSDOTY = STEMP
      RETURN
      END
C
