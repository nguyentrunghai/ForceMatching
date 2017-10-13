c########################
c  solve the least square
c  problems
c########################
      subroutine least_quare(tot_matr_size, nr_char_tp
     *, tot_matr, tot_targ
     *, solut)
c
      implicit none
      integer tot_matr_size, nr_char_tp
     *, m, n, np1, k, i,j
      double precision tot_matr, tot_targ, solut
     *, sm, srsmsq, dtmp
      dimension tot_matr( tot_matr_size, nr_char_tp)
     *, tot_targ (tot_matr_size), solut( nr_char_tp)
c
      double precision, dimension (:,:), allocatable :: a
      double precision, dimension (:), allocatable :: b, h
c
      write(*,*)''
      write(*,*)'--------------------------------'
      write(*,*)'solving the least square problem'
      write(*,*)'Reference:' 
      write(*,*)'Charles L. Lawson and Richard J. Hanson.'
      write(*,*)'SOLVING LEAST SQUARES PROBLEMS
     *, Prentice-HalL, 1974.'
      write(*,*)''
c
      m = tot_matr_size
      n = nr_char_tp
c
      allocate(a(m,n), b(m), h(n))
c
      do i=1,m
        do j=1,n
          a(i,j) =  tot_matr( i, j)
        enddo
        b(i) =  tot_targ (i)
      enddo
c
c#############################
c  apply algorithm HFT to factor the matrix a
c
      do j=1,n
      call H12(1,j,j+1,m,a(1,j),1,h(j),a(1,min(j+1,n)),1
     *,m,n-j)
      enddo
c
c  apply algorithm HS1
c  Apply the transformations  q(n)...q(1) = q to b
c replacing the previous contents of the array, b 
c
      do j=1,n
      call H12(2,j,j+1,m,a(1,j),1,h(j),b, 1, 1, 1)
      enddo
c  Solve the triangular system for the solution X.
c Store X in the array B.
c
      np1 = n+1
c
      do k=1,n
      i = np1-k
      sm = dot_product ( a(i,i+1:n), b(i+1:n) )
c
                   if ( a(i,i) == 0.0D+00 ) then
                   write ( *, '(a)' ) ' '
                   write ( *, '(a)' ) '  Terminating this case.'
                   write ( *, '(a)' ) '  A divisor is exactly zero.'
                   stop
                 end if
c
      b(i) = ( b(i) - sm ) / a(i,i)
c
      enddo
c
c  Compute the norm of the residual.
c
      srsmsq = 0.0
      do j = 1, m-n
      srsmsq = srsmsq+b(n+j)**2
      end do
      srsmsq = sqrt ( srsmsq )
c
      write(*,*)'LS Solution:'
      dtmp=0.0d0
      do i=1,n
       solut(i) = b(i)
       dtmp=dtmp+b(i)
      write(*,"(I8,F15.5)")i,solut(i)
      enddo
c      write(*,"(a10,f15.10)")'total: ',dtmp
c
      deallocate(a,b,h)
c
      write(*,*)''
      write(*,*)'solving the least square problem done!'
      write(*,*)''
c
      return
      end
c
c I copied H12 from http://www.netlib.org/lawson-hanson/all 
c
C###################################################
C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C   
C  CONSTRUCTION AND/OR APPLICATION OF A SINGLE   
C  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B   
C   
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
c                     Subroutine Arguments
c
C     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
c            Householder transformation, or Algorithm H2 to apply a
c            previously constructed transformation.
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. 
C     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO   
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M     
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
c            vector.  IUE is the storage increment between elements.  
c            On exit when MODE = 1, U() and UP contain quantities
c            defining the vector U of the Householder transformation.
c            on entry with MODE = 2, U() and UP should contain
c            quantities previously computed with MODE = 1.  These will
c            not be modified during the entry with MODE = 2.   
C     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
c            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
c            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
c            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().  
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().  
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0  
C            NO OPERATIONS WILL BE DONE ON C(). 
C     ------------------------------------------------------------------
      SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C     ------------------------------------------------------------------
      integer I, I2, I3, I4, ICE, ICV, INCR, IUE, J
      integer L1, LPIVOT, M, MODE, NCV
      double precision B, C(*), CL, CLINV, ONE, SM
c     double precision U(IUE,M)
      double precision U(IUE,*)
      double precision UP
      parameter(ONE = 1.0d0)
C     ------------------------------------------------------------------
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN    
      CL=abs(U(1,LPIVOT))   
      IF (MODE.EQ.2) GO TO 60   
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
          DO 10 J=L1,M  
   10     CL=MAX(abs(U(1,J)),CL)  
      IF (CL) 130,130,20
   20 CLINV=ONE/CL  
      SM=(U(1,LPIVOT)*CLINV)**2   
          DO 30 J=L1,M  
   30     SM=SM+(U(1,J)*CLINV)**2 
      CL=CL*SQRT(SM)   
      IF (U(1,LPIVOT)) 50,50,40     
   40 CL=-CL
   50 UP=U(1,LPIVOT)-CL 
      U(1,LPIVOT)=CL    
      GO TO 70  
C            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C   
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN  
      B= UP*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C   
      IF (B) 80,130,130 
   80 B=ONE/B   
      I2=1-ICV+ICE*(LPIVOT-1)   
      INCR=ICE*(L1-LPIVOT)  
          DO 120 J=1,NCV
          I2=I2+ICV     
          I3=I2+INCR    
          I4=I3 
          SM=C(I2)*UP
              DO 90 I=L1,M  
              SM=SM+C(I3)*U(1,I)
   90         I3=I3+ICE 
          IF (SM) 100,120,100   
  100     SM=SM*B   
          C(I2)=C(I2)+SM*UP
              DO 110 I=L1,M 
              C(I4)=C(I4)+SM*U(1,I)
  110         I4=I4+ICE 
  120     CONTINUE  
  130 RETURN
      END

