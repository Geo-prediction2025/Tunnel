module ElemCpt
	implicit none
	real(kind=8), parameter :: pi = 3.141592653589793
contains

!-----------------------------------------------------------------------------------------------------------------------------
!	Elemental matrix: Ke -----
!-----------------------------
	subroutine GetKe(tenSigma,tenSigma_S,X,Y,Z,Volume,Ke,Pe)
		real(kind=8) :: tenSigma(3,3),tenSigma_s(3,3),X(4),Y(4),Z(4),Volume,Ke(4,4),Pe(4,4)

		integer :: i,j,k
		real(kind=8) :: a(3), b(3), c(3), d(4), e(4), f(4), DEFG(4,3),DTMJacobi, coef, temp1(9,10), temp2(9), temp(10) 
		
		do i = 1, 3
			a(i) = X(i) - X(4)
			b(i) = Y(i) - Y(4)
			c(i) = Z(i) - Z(4)
		end do
		do i = 1, 3
			j = mod(i+1,3)
			k = mod(i+2,3)
			if(j==0) j = 3
			if(k==0) k = 3
			d(i) = b(j)*c(k) - b(k)*c(j)
			e(i) = c(j)*a(k) - c(k)*a(j)
			f(i) = a(j)*b(k) - a(k)*b(j)
		end do
		d(4) = 0.
		e(4) = 0.
		f(4) = 0.
		do i = 1, 3
			d(4) = d(4) - d(i)
			e(4) = e(4) - e(i)
			f(4) = f(4) - f(i)
		end do
		do i = 1, 4
			DEFG(i,1) = d(i)
			DEFG(i,2) = e(i)
			DEFG(i,3) = f(i)
		end do
		DTMJacobi = a(1)*b(2)*c(3) + b(1)*c(2)*a(3) + c(1)*a(2)*b(3) & 
				  - c(1)*b(2)*a(3) - c(2)*b(3)*a(1) - c(3)*a(2)*b(1)
		coef = Volume / DTMJacobi / DTMJacobi
		!----- ke(4*4)的下三角共10个元素，用temp2(9)*temp1(9,10)得到
		!----- temp2(9)为电阻率张量按列存储的一维数组.
		do i = 1, 3
			temp1(i,1) = DEFG(1,i) * d(1)
			temp1(i+3,1) = DEFG(1,i) * e(1)
			temp1(i+6,1) = DEFG(1,i) * f(1)	! 11

			temp1(i,2) = DEFG(2,i) * d(1)
			temp1(i+3,2) = DEFG(2,i) * e(1)
			temp1(i+6,2) = DEFG(2,i) * f(1)	! 21

			temp1(i,3) = DEFG(3,i) * d(1)
			temp1(i+3,3) = DEFG(3,i) * e(1)
			temp1(i+6,3) = DEFG(3,i) * f(1)	! 31

			temp1(i,4) = DEFG(4,i) * d(1)
			temp1(i+3,4) = DEFG(4,i) * e(1)
			temp1(i+6,4) = DEFG(4,i) * f(1)	! 41

			temp1(i,5) = DEFG(2,i) * d(2)
			temp1(i+3,5) = DEFG(2,i) * e(2)
			temp1(i+6,5) = DEFG(2,i) * f(2)	! 22

			temp1(i,6) = DEFG(3,i) * d(2)
			temp1(i+3,6) = DEFG(3,i) * e(2)
			temp1(i+6,6) = DEFG(3,i) * f(2)	! 32

			temp1(i,7) = DEFG(4,i) * d(2)
			temp1(i+3,7) = DEFG(4,i) * e(2)
			temp1(i+6,7) = DEFG(4,i) * f(2)	! 42

			temp1(i,8) = DEFG(3,i) * d(3)
			temp1(i+3,8) = DEFG(3,i) * e(3)
			temp1(i+6,8) = DEFG(3,i) * f(3)	! 33

			temp1(i,9) = DEFG(4,i) * d(3)
			temp1(i+3,9) = DEFG(4,i) * e(3)
			temp1(i+6,9) = DEFG(4,i) * f(3)	! 43

			temp1(i,10) = DEFG(4,i) * d(4)
			temp1(i+3,10) = DEFG(4,i) * e(4)
			temp1(i+6,10) = DEFG(4,i) * f(4)	! 44
		end do

		!----- Ke -----
		do j = 1, 3
			do i = 1, 3
				k = (j-1)*3 + i
				temp2(k) = tenSigma(i,j)
			end do
		end do
		temp = matmul(temp2,temp1) * coef

		k=1
		do i=1,4			 
			do j=i,4		
				Ke(j,i)=temp(k)	  !----- down-triangle
				Ke(i,j)=Ke(j,i)	  !----- up-triangle <-- down-triangle
				k=k+1
			end do
		end do

		!----- Pe -----
		do j = 1, 3
			do i = 1, 3
				k = (j-1)*3 + i
				temp2(k) = tenSigma_s(i,j)
			end do
		end do
		temp = matmul(temp2,temp1) * coef

		k=1
		do i=1,4			 
			do j=i,4		
				Pe(j,i)=temp(k)	  !----- down-triangle
				Pe(i,j)=Pe(j,i)	  !----- up-triangle <-- down-triangle
				k=k+1
			end do
		end do

	end subroutine GetKe

!------------------------------------------------------------------------------------------------------------------------
!	Compute Vp -----
!------------------------------
	subroutine AnisoHalfSpace( current0, tenRes, x_0_p, y_0_p, z_0_p, x, y, z, v )
		implicit none
		real(kind=8) :: current0, x_0_p, y_0_p, z_0_p, x, y, z, v, tenRes(3,3)

		integer :: i, j
		real(kind=8) ::  current, r1, r2, x_0_r, y_0_r, z_0_r, XX_0_p, YY_0_p, ZZ_0_p, &
						XX_0_r, YY_0_r, ZZ_0_r, eta_0_r, eta_0_p, deno, detTenRes

		detTenRes = tenRes(1,1) * tenRes(2,2) * tenRes(3,3) + tenRes(1,2) * tenRes(2,3) * tenRes(3,1) + &
					tenRes(1,3) * tenRes(3,2) * tenRes(2,1) - tenRes(1,3) * tenRes(2,2) * tenRes(3,1) - &
					tenRes(2,3) * tenRes(3,2) * tenRes(1,1) - tenRes(1,2) * tenRes(2,1) * tenRes(3,3)    
		if(detTenRes<=1.e-8) then
			write(*,*) "Subroutine AnisoHalfSpace -> Det <= 0 "
			stop
		end if
		current = current0 * dsqrt(detTenRes)

		deno = tenRes(2,2)*tenRes(1,1) - tenRes(1,2)*tenRes(1,2)
		if( dabs(deno) <= 1.e-8 ) then
			write(*,*)	"Subroutine AnisoHalfSpace -> deno == 0 "
			stop
		end if
		r1 = ( tenRes(2,2)*tenRes(1,3) - tenRes(1,2)*tenRes(2,3) ) / deno
		r2 = ( tenRes(1,1)*tenRes(2,3) - tenRes(1,2)*tenRes(1,3) ) / deno
!$$$$$
!write(*,*) r1,r2,r1+r2
!stop
		x_0_r = 2 * z_0_p * r1 + x_0_p
		y_0_r = 2 * z_0_p * r2 + y_0_p
		z_0_r = - z_0_p

		XX_0_p = x - x_0_p
		YY_0_p = y - y_0_p
		ZZ_0_p = z - z_0_p
	
		XX_0_r = x - x_0_r
		YY_0_r = y - y_0_r
		ZZ_0_r = z - z_0_r

		eta_0_p = dsqrt( tenRes(1,1)*(XX_0_p)**2 + tenRes(2,2)*(YY_0_p)**2 + tenRes(3,3)*(ZZ_0_p)**2 + &
					     2*tenRes(1,2)*XX_0_p*YY_0_p + 2*tenRes(1,3)*XX_0_p*ZZ_0_p + 2*tenRes(2,3)*YY_0_p*ZZ_0_p )
   
		eta_0_r = dsqrt( tenRes(1,1)*(XX_0_r)**2 + tenRes(2,2)*(YY_0_r)**2 + tenRes(3,3)*(ZZ_0_r)**2 + &
						2*tenRes(1,2)*XX_0_r*YY_0_r + 2*tenRes(1,3)*XX_0_r*ZZ_0_r + 2*tenRes(2,3)*YY_0_r*ZZ_0_r )

		v = current / 4.0 / pi * ( 1.0/eta_0_p + 1.0/eta_0_r )
	
	end subroutine AnisoHalfSpace

	subroutine AnisoFullSpace( current, tenRes, x_0, y_0, z_0, x, y, z, v )
		implicit none
		real(kind=8) :: current, x_0, y_0, z_0, x, y, z, v, tenRes(3,3)

		real(kind=8) :: detTenRes, coef, B

		detTenRes = tenRes(1,1) * tenRes(2,2) * tenRes(3,3) + tenRes(1,2) * tenRes(2,3) * tenRes(3,1) + &
					tenRes(1,3) * tenRes(3,2) * tenRes(2,1) - tenRes(1,3) * tenRes(2,2) * tenRes(3,1) - &
					tenRes(2,3) * tenRes(3,2) * tenRes(1,1) - tenRes(1,2) * tenRes(2,1) * tenRes(3,3)    
		if(detTenRes<=1.e-8) then
			write(*,*) "Subroutine AnisoHalfSpace -> Det <= 0 "
			stop
		end if
		coef = current * dsqrt(detTenRes) / 4.0 / pi

		B = dsqrt( tenRes(1,1)*(x-x_0)**2 + tenRes(2,2)*(y-y_0)**2 + tenRes(3,3)*(z-z_0)**2 + &
					     2*tenRes(1,2)*(x-x_0)*(y-y_0) + 2*tenRes(1,3)*(x-x_0)*(z-z_0) + 2*tenRes(2,3)*(y-y_0)*(z-z_0) )

		v = coef / B

!		write(*,*) tenRes
!		write(*,*) current, x_0, y_0, z_0, x, y, z, coef, B, v
!		pause

	end subroutine AnisoFullSpace

!-----------------------------------------------------------------------------------------------------------------------------
!	Transfer primary resistivity to tensor resistivity -----
!-----------------------------------------------------------
	subroutine pR2tR(priRho,tenRho)			 !-----	Rho = D Rho0 (D)T
		implicit none
		real(kind=8) :: a,b,c,priRho(6),tenRho(3,3),rho(3,3),D(3,3)
		integer :: i

		a = priRho(4)
		b = priRho(5)
		c = priRho(6) 
		rho = 0.
		do i = 1, 3
			rho(i,i) = priRho(i)
		end do
		D(1,1) = dcos(c) * dcos(b)
		D(2,1) = dsin(c) * dcos(a) + dcos(c) * dsin(b) * dsin(a)
		D(3,1) = dsin(c) * dsin(a) - dcos(c) * dsin(b) * dcos(a)
		D(1,2) = - dsin(c) * dcos(b)
		D(2,2) = dcos(c) * dcos(a) - dsin(c) * dsin(b) * dsin(a)
		D(3,2) = dcos(c) * dsin(a) + dsin(c) * dsin(b) * dcos(a)
		D(1,3) = dsin(b)
		D(2,3) = - dcos(b) * dsin(a)
		D(3,3) =  dcos(b) * dcos(a)

		tenRho = matmul(D,rho)
		D = transpose(D)
		tenRho = matmul(tenRho,D)
	end subroutine pR2tR

!-----------------------------------------------------------------------------------------------------------------------------
!	Get B & Bp emerged in boundary conditions	 ------
!------------------------------------------------------
	real(kind=8) function CptB(tRho,x,y,z,SourceXYZ)
		implicit none
		real(kind=8) :: tRho(3,3),x,y,z,SourceXYZ(3)
		CptB= tRho(1,1) * (x-SourceXYZ(1)) * (x-SourceXYZ(1)) + 2 * tRho(2,1) * (x-SourceXYZ(1)) * (y-SourceXYZ(2)) &
			+ 2 * tRho(3,1) * (x-SourceXYZ(1)) * (z-SourceXYZ(3)) + tRho(2,2) * (y-SourceXYZ(2)) * (y-SourceXYZ(2)) &
			+ 2 * tRho(3,2) * (y-SourceXYZ(2)) * (z-SourceXYZ(3)) + tRho(3,3) * (z-SourceXYZ(3)) * (z-SourceXYZ(3))	
	end function CptB

!-----------------------------------------------------------------------------------------------------------------------------
!	Area of a triangle   -----
!-----------------------------
	subroutine GetArea(SurX,SurY,SurZ,area)
		implicit none
		real(kind=8) :: SurX(3),SurY(3),SurZ(3),area, b(3), c(3), ri, rj, rk
		
		b(1) = SurX(1)-SurX(3)
		b(2) = SurY(1)-SurY(3)
		b(3) = SurZ(1)-SurZ(3)

		c(1) = SurX(1)-SurX(2)
		c(2) = SurY(1)-SurY(2)
		c(3) = SurZ(1)-SurZ(2)

		ri = c(2)*b(3)-b(2)*c(3)
		rj = c(3)*b(1)-b(3)*c(1)
		rk = c(1)*b(2)-b(1)*c(2)

		area = dsqrt( ri**2 + rj**2 + rk**2 ) / 2.

	end subroutine GetArea

!-----------------------------------------------------------------------------------------------------------------------------
!	Volume of a tetrahedron -----
!--------------------------------	
	subroutine GetVolume(pt_i,pt_j,pt_m,pt_l,volume)
		implicit none
		real(kind=8) :: pt_i(3), pt_j(3), pt_m(3), pt_l(3), volume
		volume = Det(pt_j,pt_m,pt_l) - Det(pt_i,pt_m,pt_l) + Det(pt_i,pt_j,pt_l) - Det(pt_i,pt_j,pt_m)
		volume = dabs(volume) / 6.
	end subroutine GetVolume

	real(kind=8) function Det(r1,r2,r3)
		implicit none
		real(kind=8) :: r1(3), r2(3), r3(3)
		Det = r1(1)*r2(2)*r3(3) + r1(2)*r2(3)*r3(1) + r1(3)*r2(1)*r3(2) - &
			  r1(3)*r2(2)*r3(1) - r2(3)*r3(2)*r1(1) - r3(3)*r2(1)*r1(2)
	end function Det

end module ElemCpt

!!!矩阵求逆，输入矩阵A，返回A的逆，仍存放在A中
 SUBROUTINE gjcp(A,N)
    integer N
    real(kind=8) A(N,N)   !存放矩阵A
      INTEGER IP(N)      !记录主列号
    REAL(kind=8) P          !工作单元，放主元
    INTEGER  I0,R      !工作单元，放主列号
    EPS=0.001

                             !     write(*,*)'gjcp  ok0000000'
       DO K=1,N
          P=0
          I0=K
          IP(K)=K
          
          DO I=K,N
             IF(ABS(A(I,K)).GT.ABS(P))THEN
                P=A(I,K)
                I0=I
                IP(K)=I
             ENDIF
          ENDDO
          
          IF(ABS(P).LE.EPS)THEN
             WRITE(*,*)'DET=0'
             stop !后来加的
             GOTO 10
          ENDIF
          
          IF(I0.NE.K)THEN
             DO J=1,N
                S=A(K,J)
                A(K,J)=A(I0,J)
                A(I0,J)=S
             ENDDO
          ENDIF
         
          A(K,K)=1./P
          
          DO I=1,N
             IF(I.NE.K)THEN
                A(I,K)=-A(I,K)*A(K,K)
                DO J=1,N
                   IF(J.NE.K)THEN
                      A(I,J)=A(I,J)+A(I,K)*A(K,J)
                     ENDIF
                ENDDO
             ENDIF
          ENDDO
         
          DO J=1,N
             IF(J.NE.K)THEN
                A(K,J)=A(K,K)*A(K,J)
             ENDIF
          ENDDO
       ENDDO
  !  write(*,*)'gjcp  ok01'    
       DO K=N-1,1,-1
          R=IP(K)
          IF(R.NE.K)THEN
             DO I=1,N
                S=A(I,R)
                A(I,R)=A(I,K)
                A(I,K)=S
             ENDDO
           ENDIF
        ENDDO
10   RETURN
     END