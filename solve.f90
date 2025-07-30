module Solver
	implicit none

contains

!------------------------------------------------------------------------------------------------------
!-----	WILL_ICCG, WILL_SICCG, WILL_SSORCG, ooo_ssor are provided in this module -----
!-----------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------------------
!********************************************************************************
!****	 ICCG solving A x = b 	< Wu X. P., 2003, Geophysical Prospecting >	*****
!********************************************************************************

!------------------	Parameters Instruction 	 -----------------------------
! Input parameter:
! n -- number of unknows
! nSA -- size of SA 
! maxIter -- maximum iteration number of this solver
! SA, IJA -- compressed storage of A
! b -- right-hand item of the equation

! Output parameter:
! x -- unknows that to be solved
! iter -- actual iteration number according to a certain tolerance
!-----------------------------------------------------------------------------

!---------------------     Special Notification ----------------------------------------------------
! 1. Constant parameter "tolerance" in subroutine "WILL_ICCG" controls convergence. You can alter it.
! 2. Comments by "@1@" and "@2@" indicate two different convergent norm. You can choose any one.
!    @1@:              2        2             @2@:           
!				|AXi+1|2 - |AXi|2                        |b-AXi|2
!              _____________________                    __________
!                           2
!                      |AXi|2                            |b-AX0|2              
!---------------------------------------------------------------------------------------------------    
	
subroutine WILL_ICCG( n, nSA, maxIter, SA, IJA, b, x, iter )
      implicit none
      integer :: n, nSA, maxIter, IJA(*), iter
      real(kind=8) :: SA(*), b(n), x(n)

      integer :: i, j, k, negtiveDi
      real(kind=8), parameter :: tolerance = 1.0D-8		!$$$$$ Tolerance controling $$$$$
      real(kind=8) :: D(n), r(n), p(n), q(n), q0, q1, alpha, beta, sr, srr, spq
      
!      real(kind=8) :: qq, qqq, qqqqq	! @1@
      real(kind=8) :: s(n), L2, L2r0	! @2@

	  write(*,*)	"WILL_ICCG solving...( L2-norm, |b-Axi| / |b-Ax0| )."
      s = b				! @2@

!---  Get diagonal matrix D --------------------------------------
      D(1) = SA(1)
      do i = 2, n
	    D(i) = 0.0
	    do k = IJA(i), IJA(i+1) - 1
			j = IJA(k)
			if( j /= 0 ) then
				D(i) = D(i) + SA(k)**2 / D(j)
			end if
	    end do
	    D(i) = SA(i) - D(i)
      end do

!--- Find out how many negtive entries in matrix D ---------------
      negtiveDi = 0
      do i = 1, n
	    if( D(i)<= 0.0 ) negtiveDi = negtiveDi + 1
      end do
      write(*,*)	"	There are",real(negtiveDi)/n*100,"% negtive entries of the diagonal of matrix D."
      read(*,*)

!	  Initializing
      iter = 0

      x = 0.
      r = b

      q0 = 0.
      do i = 1, n
	    q0 = q0 + b(i)**2
      end do
      
!      qqq = q0		! @1@
      L2r0 = dsqrt(q0)	! @2@

      call ICCG_LLTx( n, SA, IJA, D, r, p )		! p <= (LL[T])(-1) r
      ! r changed after calling LLTx
      r = b

      call ICCG_Ax( n, SA, IJA, p, q )		! q <= A p
      sr = 0.
      do i = 1, n
	    sr = sr + r(i) * p(i)
      end do

      do while( .true. )
	    iter = iter + 1
	    if( iter >= maxIter ) then
		  write(*,*) ( x(i), i = 1, n )
		  return
	    end if
	    spq = 0.
	    do i = 1, n
		  spq = spq + p(i) * q(i)
	    end do
	    alpha = sr / spq

	    do i = 1, n
		  x(i) = x(i) + alpha * p(i)
	    end do
	    call ICCG_Ax( n, SA, IJA, x, b )		! b <= A x

	    q1 = 0.
	    do i = 1, n
		  q1 = q1 + b(i)**2
	    end do

	    L2 = 0.				!@2@
	    do i = 1, n				!@2@
		  L2 = L2 + ( s(i) - b(i) )**2	!@2@
	    end do				!@2@
	    L2 = dsqrt(L2)			!@2@	

!	    qq = dabs( ( q1-q0 ) / q0 )		!@1@
!	    qqqqq = dabs( ( q1-qqq ) / qqq )	!@1@
!	    write(*,*)	iter, qq, qqqqq		!@1@
	    write(*,*) iter, L2, L2 / L2r0	!@2@
!	    if( dsqrt(qq) <= tolerance ) then	!@1@
	    if( ( L2 / L2r0 ) <= tolerance ) then	!@2@
		  write(*,*) "WILL_ICCG Finished."
		  write(*,*)
		  return
	    end if
	    
	    do i = 1, n
		  r(i) = r(i) - alpha * q(i)
		  b(i) = r(i)
	    end do
	    call ICCG_LLTx( n, SA, IJA, D, r, q )

	    r = b
	    srr = 0.
	    do i = 1, n
		  srr = srr + r(i) * q(i)
	    end do
	    beta = srr / sr

	    do i = 1, n
		  p(i) = q(i) + beta * p(i)
	    end do
	    call ICCG_Ax( n, SA, IJA, p, q )
      
	    sr = srr
	    q0 = q1
      end do
end subroutine WILL_ICCG

subroutine ICCG_Ax( n, SA, IJA, x, b )	! b <= A x
      implicit none
      integer :: n, IJA(*)
      real(kind=8) :: SA(*), x(n), b(n)

      integer :: i, j, k

      if( IJA(1) .ne. (n+2) ) then
	    write(*,*) "Mismatch vector and matrix"
	    stop
      end if

      do i = 1, n
	    b(i) = SA(i) * x(i)
	    do k = IJA(i), IJA(i+1) - 1
			if( IJA(k) /= 0 ) then
				b(i) = b(i) + SA(k) * x(IJA(k))
			end if
	    end do
      end do

      do i = 1, n
	    do k = IJA(i), IJA(i+1) - 1
		  j = IJA(k)
		  if( j /= 0 ) then
				b(j) = b(j) + SA(k) * x(i)
		  end if	
	    end do
      end do

end subroutine ICCG_Ax

subroutine ICCG_LLTx( n, SA, IJA, D, b, x )		! x <= (LL[T])(-1) b
      implicit	none
      integer :: n, IJA(*)
      real(kind=8) :: SA(*), D(n), b(n), x(n)

      integer :: i, j, k

      if( IJA(1) .ne. (n+2) ) then
	    write(*,*) "Mismatch vector and matrix"
	    stop
      end if

      x(1) = b(1) / D(1)
      do i = 1, n
	    do k = IJA(i), IJA(i+1) - 1
			if(IJA(k) /= 0 ) then
				b(i) = b(i) - SA(k) * x(IJA(k))
			end if
	    end do
	    x(i) = b(i) / D(i)
      end do

      do i = 1, n
            b(i) = x(i) * d(i)
      end do

      x(n) = b(n) / D(n)
      do i = n, 2, -1
	    do k = IJA(i), IJA(i+1) - 1
		  j = IJA(k)
		  if( j /= 0 ) then
				b(j) = b(j) - SA(k) * x(i)
		  end if
	    end do
	    x(i-1) = b(i-1) / D(i-1)
      end do

end subroutine ICCG_LLTx 

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!********************************************************************************
!****	 SICCG solving A x = b 	< Wu X. P., 2003, Geophysical Prospecting >	*****
!---- SICCG is based on ICCG, just multiply (1+MU) on diagonal of A, or divide 
!	  (1+MU) on off-diagonal of A
!********************************************************************************

!------------------	Parameters Instruction 	 -----------------------------
! Input parameter:
! n -- number of unknows
! nSA -- size of SA 
! maxIter -- maximum iteration number of this solver
! SA, IJA -- compressed storage of A
! b -- right-hand item of the equation

! Output parameter:
! x -- unknows that to be solved
! iter -- actual iteration number according to a certain tolerance
!-----------------------------------------------------------------------------

!---------------------     Special Notification ----------------------------------------------------
! 1. Constant parameter "tolerance" in subroutine "WILL_SICCG" controls convergence. You can alter it.
! 2. Comments by "@1@" and "@2@" indicate two different convergent norm. You can choose any one.
!    @1@:              2        2             @2@:           
!				|AXi+1|2 - |AXi|2                        |b-AXi|2
!              _____________________                    __________
!                           2
!                      |AXi|2                            |b-AX0|2              
!--------------------------------------------------------------------------------------------------- 
subroutine WILL_SICCG( n, nSA, maxIter, SA, IJA, b, x, iter )
      implicit none
      integer :: n, nSA, maxIter, IJA(*), iter
      real(kind=8) :: SA(*), b(n), x(n)

      integer :: i, j, k, negtiveDi
      real(kind=8), parameter :: tolerance = 1.0D-8		!$$$$$ Tolerance controling $$$$$
	  REAL(KIND=8), PARAMETER :: MU = 0.4
      real(kind=8) :: D(n), r(n), p(n), q(n), q0, q1, alpha, beta, sr, srr, spq
      
!      real(kind=8) :: qq, qqq, qqqqq	! @1@
      real(kind=8) :: s(n), L2, L2r0	! @2@

	  write(*,*)	"WILL_SICCG solving...( L2-norm, |b-Axi| / |b-Ax0| )."
      s = b				! @2@

      D(1) = SA(1)
      do i = 2, n
	    D(i) = 0.0
	    do k = IJA(i), IJA(i+1) - 1
			j = IJA(k)
			if( j /= 0)	then
				D(i) = D(i) + (SA(k)/(1+MU))**2 / D(j)
			end if
	    end do
	    D(i) = SA(i) - D(i)
      end do

      negtiveDi = 0
      do i = 1, n
	    if( D(i)<= 0.0 ) negtiveDi = negtiveDi + 1
      end do
      write(*,*)	"	There are",real(negtiveDi)/n*100,"% negtive entries of the diagonal of matrix D."
      read(*,*)

      ! Initializing
      iter = 0

      x = 0.
      r = b

      q0 = 0.
      do i = 1, n
	    q0 = q0 + b(i)**2
      end do
      
!      qqq = q0		! @1@
      L2r0 = dsqrt(q0)	! @2@

      call SICCG_LLTx( n, SA, IJA, D, r, p, MU )		! p <= (LL[T])(-1) r
      ! r changed after calling LLTx
      r = b

      call SICCG_Ax( n, SA, IJA, p, q )		! q <= A p
      sr = 0.
      do i = 1, n
	    sr = sr + r(i) * p(i)
      end do

      do while( .true. )
	    iter = iter + 1
	    if( iter >= maxIter ) then
		  write(*,*) ( x(i), i = 1, n )
		  return
	    end if
	    spq = 0.
	    do i = 1, n
		  spq = spq + p(i) * q(i)
	    end do
	    alpha = sr / spq

	    do i = 1, n
		  x(i) = x(i) + alpha * p(i)
	    end do
	    call SICCG_Ax( n, SA, IJA, x, b )		! b <= A x

	    q1 = 0.
	    do i = 1, n
		  q1 = q1 + b(i)**2
	    end do

	    L2 = 0.				!@2@
	    do i = 1, n				!@2@
		  L2 = L2 + ( s(i) - b(i) )**2	!@2@
	    end do				!@2@
	    L2 = dsqrt(L2)			!@2@	

!	    qq = dabs( ( q1-q0 ) / q0 )		!@1@
!	    qqqqq = dabs( ( q1-qqq ) / qqq )	!@1@
!	    write(*,*)	iter, qq, qqqqq		!@1@
	    write(*,*) iter, L2, L2 / L2r0	!@2@
!	    if( dsqrt(qq) <= tolerance ) then	!@1@
		if( ( L2 / L2r0 ) <= tolerance ) then	!@2@
		  write(*,*) "WILL_SICCG Finished."
		  write(*,*)
		  return
	    end if
	    
	    do i = 1, n
		  r(i) = r(i) - alpha * q(i)
		  b(i) = r(i)
	    end do
	    call SICCG_LLTx( n, SA, IJA, D, r, q, MU )

	    r = b
	    srr = 0.
	    do i = 1, n
		  srr = srr + r(i) * q(i)
	    end do
	    beta = srr / sr

	    do i = 1, n
		  p(i) = q(i) + beta * p(i)
	    end do
	    call SICCG_Ax( n, SA, IJA, p, q )
      
	    sr = srr
	    q0 = q1
      end do
end subroutine WILL_SICCG

subroutine SICCG_Ax( n, SA, IJA, x, b )	! b <= A x
      implicit none
      integer :: n, IJA(*)
      real(kind=8) :: SA(*), x(n), b(n)

      integer :: i, j, k

      if( IJA(1) .ne. (n+2) ) then
	    write(*,*) "Mismatch vector and matrix"
	    stop
      end if

      do i = 1, n
	    b(i) = SA(i) * x(i)
	    do k = IJA(i), IJA(i+1) - 1
			if(IJA(k) /= 0)		then
				b(i) = b(i) + SA(k) * x(IJA(k))
			end if
	    end do
      end do

      do i = 1, n
	    do k = IJA(i), IJA(i+1) - 1
		  j = IJA(k)
		  if(j/=0) then
			b(j) = b(j) + SA(k) * x(i)
		  end if
	    end do
      end do

end subroutine SICCG_Ax

subroutine SICCG_LLTx( n, SA, IJA, D, b, x, MU )		! x <= (LL[T])(-1) b
      implicit	none
      integer :: n, IJA(*)
      real(kind=8) :: SA(*), D(n), b(n), x(n), MU

      integer :: i, j, k

      if( IJA(1) .ne. (n+2) ) then
	    write(*,*) "Mismatch vector and matrix"
	    stop
      end if

      x(1) = b(1) / D(1)
      do i = 1, n
	    do k = IJA(i), IJA(i+1) - 1
			if(IJA(k)/=0) then
				b(i) = b(i) - (SA(k)/(1+MU)) * x(IJA(k))
			end if
	    end do
	    x(i) = b(i) / D(i)
      end do

      do i = 1, n
            b(i) = x(i) * d(i)
      end do

      x(n) = b(n) / D(n)
      do i = n, 2, -1
	    do k = IJA(i), IJA(i+1) - 1
		  j = IJA(k)
		  if(j/=0) then
				b(j) = b(j) - (SA(k)/(1+MU)) * x(i)
		  end if
	    end do
	    x(i-1) = b(i-1) / D(i-1)
      end do

end subroutine SICCG_LLTx 
 
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!********************************************
!****	 SSOR-PCG solving A x = b 	*****
! < according to Saad Y., 2000, Iterative methods for sparse linear systems >
! The preconditioner is M = (D+wE) D[-1] (D+wF), where D, E, F are diagonal of A, strict lower and upper
! part of A, respectively.
!********************************************

!------------------	Parameters Instruction 	 -----------------------------
! Input parameter:
! n -- number of unknows
! nSA -- size of SA 
! maxIter -- maximum iteration number of this solver
! SA, IJA -- compressed storage of A
! b -- right-hand item of the equation

! Output parameter:
! x -- unknows that to be solved
! iter -- actual iteration number according to a certain tolerance
!-----------------------------------------------------------------------------

!---------------------     Special Notification ----------------------------------------------------
! 1. Constant parameter "tolerance" in subroutine "WILL_SSORCG" controls convergence. You can alter it.
! 2. Comments by "@1@" and "@2@" indicate two different convergent norm. You can choose any one.
!    @1@:              2        2             @2@:           
!				|AXi+1|2 - |AXi|2                        |b-AXi|2
!              _____________________                    __________
!                           2
!                      |AXi|2                            |b-AX0|2              
!---------------------------------------------------------------------------------------------------    

subroutine WILL_SSORCG( n, nSA, maxIter, SA, IJA, b, x, iter )
      implicit none
      integer :: n, nSA, maxIter, IJA(*), iter
      real(kind=8) :: SA(*), b(n), x(n)

      integer :: i, j, k
      real(kind=8), parameter :: tolerance = 1.0D-8		!$$$$$ Tolerance controling $$$$$
      real(kind=8) :: D(n), r(n), p(n), q(n), q0, q1, alpha, beta, sr, srr, spq, w
    
!      real(kind=8) :: qq, qqq, qqqqq	! @1@
      real(kind=8) :: s(n), L2, L2r0	! @2@

	  write(*,*)	"WILL_SSORCG solving...( L2-norm, |b-Axi| / |b-Ax0| )."
      s = b				! @2@

      w = 1.4

      do i = 1, n
	    D(i) = SA(i)
      end do
  
      ! Initializing
      iter = 0

      x = 0.
      r = b

      q0 = 0.
      do i = 1, n
	    q0 = q0 + b(i)**2
      end do
      
!      qqq = q0		! @1@
      L2r0 = dsqrt(q0)	! @2@

      call ssor_LLTx( n, w, SA, IJA, D, r, p )		! p <= (LL[T])(-1) r
      ! r changed after calling LLTx
      r = b

      call ssor_Ax( n, SA, IJA, p, q )		! q <= A p
      sr = 0.
      do i = 1, n
	    sr = sr + r(i) * p(i)
      end do

      do while( .true. )
	    iter = iter + 1
	    if( iter >= maxIter ) then
		  write(*,*) ( x(i), i = 1, n )
		  return
	    end if
	    spq = 0.
	    do i = 1, n
		  spq = spq + p(i) * q(i)
	    end do
	    alpha = sr / spq

	    do i = 1, n
		  x(i) = x(i) + alpha * p(i)
	    end do
	    call ssor_Ax( n, SA, IJA, x, b )		! b <= A x

	    q1 = 0.
	    do i = 1, n
		  q1 = q1 + b(i)**2
	    end do

	    L2 = 0.				!@2@
	    do i = 1, n				!@2@
		  L2 = L2 + ( s(i) - b(i) )**2	!@2@
	    end do				!@2@
	    L2 = dsqrt(L2)			!@2@	

!	    qq = dabs( ( q1-q0 ) / q0 )		!@1@
!	    qqqqq = dabs( ( q1-qqq ) / qqq )	!@1@
!	    write(*,*)	iter, qq, qqqqq		!@1@
!	    write(*,*) iter, L2, L2 / L2r0	!@2@	-----$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	    if( dsqrt(qq) <= tolerance ) then	!@1@
		if( ( L2 / L2r0 ) <= tolerance ) then	!@2@
		  write(*,*) "WILL_SSORCG Finished."
		  write(*,*)
		  return
	    end if
	    
	    do i = 1, n
		  r(i) = r(i) - alpha * q(i)
		  b(i) = r(i)
	    end do
	    call ssor_LLTx( n, w, SA, IJA, D, r, q )

	    r = b
	    srr = 0.
	    do i = 1, n
		  srr = srr + r(i) * q(i)
	    end do
	    beta = srr / sr

	    do i = 1, n
		  p(i) = q(i) + beta * p(i)
	    end do
	    call ssor_Ax( n, SA, IJA, p, q )
      
	    sr = srr
	    q0 = q1
      end do
end subroutine WILL_SSORCG

subroutine ssor_Ax( n, SA, IJA, x, b )	! b <= A x
      implicit none
      integer :: n, IJA(*)
      real(kind=8) :: SA(*), x(n), b(n)

      integer :: i, j, k

      if( IJA(1) .ne. (n+2) ) then
	    write(*,*) "Mismatch vector and matrix"
	    stop
      end if

      do i = 1, n
	    b(i) = SA(i) * x(i)
	    do k = IJA(i), IJA(i+1) - 1
			if( IJA(k) /= 0 ) then				!!!!!!!$$$$$$$$$$$$$$$$$$$$$$$
				b(i) = b(i) + SA(k) * x(IJA(k))
			end if
	    end do
      end do

      do i = 1, n
	    do k = IJA(i), IJA(i+1) - 1
		  j = IJA(k)
		  if( j /= 0 ) then						!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			b(j) = b(j) + SA(k) * x(i)
		  end if
	    end do
      end do

end subroutine ssor_Ax

subroutine ssor_LLTx( n, w, SA, IJA, D, b, x )		! x <= (LL[T])(-1) b
      implicit	none
      integer :: n, IJA(*)
      real(kind=8) :: w, SA(*), D(n), b(n), x(n)

      integer :: i, j, k

      if( IJA(1) .ne. (n+2) ) then
	    write(*,*) "Mismatch vector and matrix"
	    stop
      end if

      x(1) = b(1) / D(1)
      do i = 1, n
	    do k = IJA(i), IJA(i+1) - 1
			if( IJA(k) /= 0 ) then
				b(i) = b(i) - w * SA(k) * x(IJA(k))
			end if
	    end do
	    x(i) = b(i) / D(i)
      end do

      do i = 1, n
            b(i) = x(i) * d(i)
      end do

      x(n) = b(n) / D(n)
      do i = n, 2, -1
	    do k = IJA(i), IJA(i+1) - 1
		  j = IJA(k)
		  if( j /= 0 )	then
			b(j) = b(j) - w * SA(k) * x(i)
		  end if
	    end do
	    x(i-1) = b(i-1) / D(i-1)
      end do

end subroutine ssor_LLTx 



 
!------------------------------------------------------------------------------------------------------
!	Old version of SSOR-CG according to Klaus S. and Wu X.P., the preconditioner is M = (I+wE)(I+wF),
!where I is the indentity matrix, E=F[T].
!	The original code of 'ooo_ssor' is from Wu X.P. Little revision is made.
!------------------------------------------------------------------------------------------------------    

subroutine ooo_ssor( n, nSA, maxIter, SA, IJA, b, x, iter )
      implicit none
      integer :: n, nSA, maxIter, IJA(*), iter
      real(kind=8) :: SA(*), b(n), x(n)

      integer :: i, j, k
      real(kind=8), allocatable :: SB(:)
      real(kind=8), parameter :: tolerance = 1.0D-8		!$$$$$ Tolerance controling $$$$$
      real(kind=8) :: D(n), r(n), p(n), q(n), q0, q1, alpha, beta, sr, srr, spq
      
!      real(kind=8) :: qq, qqq, qqqqq	! @1@
      real(kind=8) :: s(n), L2, L2r0	! @2@
      s = b				! @2@

      allocate( SB(nSA) )
      SB = 1.4 * SA(1:nSA)			! Relaxation factor w = 1.4

      do i = 1, n
	    D(i) = 1.4 * ( SA(i) - 1.0 ) / 2.0 + 1.0  ! Relaxation factor w = 1.4
      end do
  
      ! Initializing
      iter = 0

      x = 0.
      r = b

      q0 = 0.
      do i = 1, n
	    q0 = q0 + b(i)**2
      end do
      
!      qqq = q0		! @1@
      L2r0 = dsqrt(q0)	! @2@

      call ooo_LLTx( n, SB, IJA, D, r, p )		! p <= (LL[T])(-1) r
      ! r changed after calling LLTx
      r = b

      call ooo_Ax( n, SA, IJA, p, q )		! q <= A p
      sr = 0.
      do i = 1, n
	    sr = sr + r(i) * p(i)
      end do

      do while( .true. )
	    iter = iter + 1
	    if( iter >= maxIter ) then
		  write(*,*) ( x(i), i = 1, n )
		  return
	    end if
	    spq = 0.
	    do i = 1, n
		  spq = spq + p(i) * q(i)
	    end do
	    alpha = sr / spq

	    do i = 1, n
		  x(i) = x(i) + alpha * p(i)
	    end do
	    call ooo_Ax( n, SA, IJA, x, b )		! b <= A x

	    q1 = 0.
	    do i = 1, n
		  q1 = q1 + b(i)**2
	    end do

	    L2 = 0.				!@2@
	    do i = 1, n				!@2@
		  L2 = L2 + ( s(i) - b(i) )**2	!@2@
	    end do				!@2@
	    L2 = dsqrt(L2)			!@2@	

!	    qq = dabs( ( q1-q0 ) / q0 )		!@1@
!	    qqqqq = dabs( ( q1-qqq ) / qqq )	!@1@
!	    write(*,*)	iter, qq, qqqqq		!@1@
	    write(*,*) iter, L2, L2 / L2r0	!@2@
!	    if( dsqrt(qq) <= tolerance ) then	!@1@
	    if( ( L2 / L2r0 ) <= tolerance ) then	!@2@
		  write(*,*) "Finished."
		  write(*,*)
		  return
	    end if
	    
	    do i = 1, n
		  r(i) = r(i) - alpha * q(i)
		  b(i) = r(i)
	    end do
	    call ooo_LLTx( n, SA, IJA, D, r, q )

	    r = b
	    srr = 0.
	    do i = 1, n
		  srr = srr + r(i) * q(i)
	    end do
	    beta = srr / sr

	    do i = 1, n
		  p(i) = q(i) + beta * p(i)
	    end do
	    call ooo_Ax( n, SA, IJA, p, q )
      
	    sr = srr
	    q0 = q1
      end do
end subroutine ooo_ssor

subroutine ooo_Ax( n, SA, IJA, x, b )	! b <= A x
      implicit none
      integer :: n, IJA(*)
      real(kind=8) :: SA(*), x(n), b(n)

      integer :: i, j, k

      if( IJA(1) .ne. (n+2) ) then
	    write(*,*) "Mismatch vector and matrix"
	    stop
      end if

      do i = 1, n
	    b(i) = SA(i) * x(i)
	    do k = IJA(i), IJA(i+1) - 1
			if(IJA(k) /= 0 ) then
				b(i) = b(i) + SA(k) * x(IJA(k))
            end if
	    end do
      end do

      do i = 1, n
	    do k = IJA(i), IJA(i+1) - 1
		  j = IJA(k)
		  if(j/=0) then
			 b(j) = b(j) + SA(k) * x(i)
		  end if
	    end do
      end do

end subroutine ooo_Ax

subroutine ooo_LLTx( n, SA, IJA, D, b, x )		! x <= (LL[T])(-1) b
      implicit	none
      integer :: n, IJA(*)
      real(kind=8) :: SA(*), D(n), b(n), x(n)

      integer :: i, j, k

      if( IJA(1) .ne. (n+2) ) then
	    write(*,*) "Mismatch vector and matrix"
	    stop
      end if

      x(1) = b(1) / D(1)
      do i = 1, n
	    do k = IJA(i), IJA(i+1) - 1
			if(IJA(k) /= 0 )  then
				b(i) = b(i) - SA(k) * x(IJA(k))
			end if
	    end do
	    x(i) = b(i) / D(i)
      end do

      x(n) = b(n) / D(n)
      do i = n, 2, -1
	    do k = IJA(i), IJA(i+1) - 1
		  j = IJA(k)
		  if(j/=0)  then
				b(j) = b(j) - SA(k) * x(i)
		  end if
	    end do
	    x(i-1) = b(i-1) / D(i-1)
      end do

end subroutine ooo_LLTx


end module Solver