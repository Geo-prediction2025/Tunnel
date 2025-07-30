module ResProc
	implicit none
	real(kind=8) , parameter :: pi = 3.141592653589793
contains

!-------------------------------------------------------------------------------------------------------------------------------
!	Declaration -----
!--------------------
	subroutine Declaration()
		implicit none

		write( *, "(6X,'|-------------   USGDC4.1 - By WangWei  ---------------|')" )
		write(*,*)
	
    end subroutine Declaration

    subroutine fanchangoutput(current, NNode, coor, SourceXYZ, Vp, v)	
		implicit none
		integer :: NNode
		real(kind=8) :: current, coor(3,Nnode), SourceXYZ(3),Vp(NNode), v(nNode), K, Acoor(3), Mcoor(3), Ncoor(3), AM, AN, MN
		integer :: i

        Acoor = SourceXYZ
        MN = 20
        
		open(unit=2,file="fcresult.txt")
		do i = 10, 19
            Mcoor = coor(:,i)
            Ncoor = coor(:,i+1)
            AM = dsqrt( sum( ( Acoor - Mcoor ) * ( Acoor - Mcoor ) ) )
			AN = dsqrt( sum( ( Acoor - Ncoor ) * ( Acoor - Ncoor ) ) )
			K = 4 * pi * AM * AN / MN
            write(2,"(8F15.5)") K * (Vp(i)-Vp(i+1)) / current
        enddo
        
            write(2,*) 
            
		do i = 21, 30
            Mcoor = coor(:,i)
            Ncoor = coor(:,i+1)
            AM = dsqrt( sum( ( Acoor - Mcoor ) * ( Acoor - Mcoor ) ) )
			AN = dsqrt( sum( ( Acoor - Ncoor ) * ( Acoor - Ncoor ) ) )
			K = 4 * pi * AM * AN / MN
            write(2,"(8F15.5)") K * (Vp(i)-Vp(i+1)) / current
        enddo
        close(2)
    end subroutine fanchangoutput

    subroutine cqres(current, NNode, coor, sourcePt, Vp, v, beta,kamn,res)	
		implicit none
		integer :: NNode,sourcePt
		real(kind=8) :: current, coor(3,Nnode),Vp(NNode), v(nNode), beta(nNode), kamn(:), res(:)
		integer :: i

		open(unit=2,file="cqres.txt")
		do i = sourcePt+1,sourcePt+99   ! tunnel:x=x+8  || no tunnel:x=x-8
            res(i-sourcePt) = kamn(i-sourcePt)*beta(i-sourcePt)*((Vp(i)+v(i))-(Vp(i+1)+v(i+1)))*10
            write(2,"(8F15.5)") res(i-sourcePt)
        enddo
        close(2)
    end subroutine cqres
    
    subroutine betaoutput(current, NNode, coor, sourcePt, Vp, v, beta)	
		implicit none
		integer :: NNode,sourcePt
		real(kind=8) :: current, coor(3,Nnode), Vp(NNode), v(nNode), beta(:)
		integer :: i

		open(unit=2,file="beta.txt")
		do i = sourcePt+1,sourcePt+99   ! tunnel:x=x+8  || no tunnel:x=x-8
            beta(i-sourcePt) = (Vp(i)-Vp(i+1))/((Vp(i)+v(i))-(Vp(i+1)+v(i+1)))
            write(2,"(8F15.5)") beta(i-sourcePt)
        enddo
        close(2)
    end subroutine betaoutput
    
    subroutine cqtcoutput(current, NNode, coor, sourcePt, Vp, v)	
		implicit none
		integer :: NNode,sourcePt
		real(kind=8) :: current, coor(3,Nnode),Vp(NNode), v(nNode)
		integer :: i

		open(unit=2,file="cqresult.txt")
		do i = sourcePt+1,sourcePt+100   ! tunnel:x=x+8  || no tunnel:x=x-8
            write(2,"(8F15.5)") (Vp(i)+v(i))*10
            !write(2,"(8F15.5)") Vp(i)*10
        enddo
        close(2)
    end subroutine cqtcoutput
    
!-------------------------------------------------------------------------------------------------------------------------------
!	Time & Iteration -----
!-------------------------
subroutine TimeConsumed( startTime, endTime, INT4, iter )
	implicit none

	integer :: startTime, endTime, diffTime, INT4, iter

	diffTime = endTime - startTime
	open(unit=1, file = "time_iter.txt")
	write(1,*)	"Iteration:", iter, ";		Time consumed:", diffTime/INT4, "seconds"
	write(*,*)	"Iteration:", iter, ";		Time consumed:", diffTime/INT4, "seconds"
end subroutine TimeConsumed

!-------------------------------------------------------------------------------------------------------------------------------
!	nodes' potentials on the ground -----
!----------------------------------------
	subroutine SurfaceP(current, NNode, coor, SourceXYZ, Vp, Vs)
		implicit none
		integer :: NNode
		real(kind=8) :: x, y, AM, current, coef, ra, coor(3,Nnode), SourceXYZ(3), Vp(NNode), Vs(nNode)
		integer :: i

		open(unit=2,file="surfaceP.dat")
		coef = 2 * pi / current
		do i = 1, Nnode
			if( coor(3,i)==0 ) then
				x = coor(1,i)
				y = coor(2,i)
				AM = dsqrt( (x-SourceXYZ(1))**2 + (y-SourceXYZ(2))**2 )
				ra = coef * AM * ( Vp(i) + Vs(i) )
				write(2,"(6F20.10)") x, y, Vp(i), Vs(i), ( Vp(i)+Vs(i) ), ra
			end if
		end do
	end subroutine SurfaceP

!-------------------------------------------------------------------------------------------------------------------------------
!	nodes' potentials in interested area on the ground -----
!-----------------------------------------------------------
	subroutine InterestedSurfacePotential(current, nNode, coor, SourceXYZ, Vp, Vs)
		implicit none
		integer :: nNode, i
		real(kind=8) :: x, y, z, AM, current, coef, ra, coor(3,nNode), sourceXYZ(3), Vp(nNode), Vs(nNode), V

		open(unit=2, file='ISP.dat')
		open(unit=3, file='IXP.dat')
		open(unit=4, file='IYP.dat')
		coef = 2 * pi / current
		 
		do i = 773,3173
			x = coor(1,i)
			y = coor(2,i)
			z = coor(3,i)
			AM = dsqrt( (x-SourceXYZ(1))**2 + (y-SourceXYZ(2))**2 )
			V = Vp(i) + Vs(i)
			ra = coef * AM * V
			write(2,"(7F20.10)") x, y, z, Vp(i), Vs(i), V,  ra
			if(y==0) then
				write(3,"(7F20.10)") x, y, z, Vp(i), Vs(i), V, ra
			end if
			if(x==0) then
				write(4,"(7F20.10)") x, y, z, Vp(i), Vs(i), V, ra
			end if
		end do
	end subroutine InterestedSurfacePotential

	subroutine DplDpl(nE, potential1, potential2, potential3, potential4)
		implicit none
		integer :: nE
		real(kind=8) :: potential1(nE,nE), potential2(nE,nE), potential3(nE,nE), potential4(nE,nE)

		integer :: i
		real(kind=8) :: AM, AN, coef(4), U, appR, depth

		open(unit=1, file="BH2.dat")
		open(unit=2, file="BH3.dat")
		open(unit=3, file="BH4.dat")
		open(unit=4, file="BH5.dat")

		AM = 110
		AN = sqrt(AM**2+1**2)
		coef(1) = 2*pi/(1/AM-1/AN)

		AM = 55*sqrt(2.0)
		AN = sqrt(AM**2+1**2)
		coef(2) = 2*pi/(1/AM-1/AN)
		coef(3) = coef(2)

		AM = 55
		AN = sqrt(AM**2+1**2)
		coef(4) = 2*pi/(1/AM-1/AN)

		do i = 1, nE-1
			U = potential1(i,i) - potential1(i+1,i) + potential1(i+1,i+1) - potential1(i,i+1)
			depth = 900.5 + (i-1)*1.0
			appR = coef(1) * U
	!		write(1,*) potential1(i,i) , potential1(i+1,i) , potential1(i+1,i+1) , potential1(i,i+1)
			write(1,*) depth, appR
		end do

		do i = 1, nE-1
			U = potential2(i,i) - potential2(i+1,i) + potential2(i+1,i+1) - potential2(i,i+1)
			depth = 900.5 + (i-1)*1.0
			appR = coef(2) * U
			write(2,*) depth, appR
		end do

		do i = 1, nE-1
			U = potential3(i,i) - potential3(i+1,i) + potential3(i+1,i+1) - potential3(i,i+1)
			depth = 900.5 + (i-1)*1.0
			appR = coef(3) * U
			write(3,*) depth, appR
		end do

		do i = 1, nE-1
			U = potential4(i,i) - potential4(i+1,i) + potential4(i+1,i+1) - potential4(i,i+1)
			depth = 900.5 + (i-1)*1.0
			appR = coef(4) * U
			write(4,*) depth, appR
		end do

	end subroutine DplDpl

end module ResProc

