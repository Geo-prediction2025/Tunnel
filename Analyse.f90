
module	ANALYSE
	use ElemCpt
	implicit none
contains

!-----------------------------------------------------------------------------------------------------------------------------
!	Assemble the system matrix ----------
!----------------------------------------
	subroutine GetSigmaK(current,nR,tenRho_temp,Nnode,coor,Nelem,elem,Nface,face,Vp,SigmaP,SA,IJA,rowWidth,sourcePt,SourceXYZ)
		implicit none

		integer :: nR, Nnode, Nelem, Nface, elem(5,Nelem), face(4,Nface), IJA(*), rowWidth, sourcePt
		real(kind=8) :: current, rho(6,nR), coor(3,Nnode), Vp(nNode) ,SigmaP(Nnode), SA(*)

		real(kind=8) :: VolX(4),VolY(4),VolZ(4),SurX(3),SurY(3),SurZ(3),B,tenRho(3,3),tenRho_temp(3,3,nR),&
						tenSigma(3,3),tenSigma_temp(3,3,nR),Ke(4,4),area,Volume,xx,yy,zz,SourceXYZ(3),RR, tenRho_p(3,3), &
						tenSigma_p(3,3), tenSigma_s(3,3), x, y, z,P(4),eVp(4),P2(3),Pe2(3,3),eVp2(3),coef, coef1, coef2, &
						Pe(4,4), Bp
		integer :: i,j,m,mm,nn,node1,node2,pos(Nnode),NT, Attr, layerID, faceID

		real(kind=8) :: MirrorSourceXYZ(3), deno, s, t, BB, q, qp, BBp

		do i = 1, nR
			tenSigma_temp(:,:,i) = tenRho_temp(:,:,i)
			call gjcp(tenSigma_temp(:,:,i),3)
		end do
		tenRho_p = tenRho_temp(:,:,1)
		tenSigma_p = tenSigma_temp(:,:,1)

		!----- Output all the resistivity the program used
		do i = 1, nR
		!do i = 1, 1
			write(*,"(' Resistivities and Conductivities', I4)") i
			write(*,"(6F16.2)") ((tenRho_temp(j,:,i),tenSigma_temp(j,:,i)),j=1,3)
		end do
		write(*,*)


		!----- Forming right-hand item of the FE linear equations
		SigmaP = 0.0			
!		sourcePt = 1483					! $$$$$$ PAY ATTENTION TO IT $$$$$$$

		do i=1,3
			SourceXYZ(i) = coor(i,sourcePt)
		end do
		write(*, "(' Source Position:', 3F15.3)")	SourceXYZ
!		pause

		!----- Compute the location of mirror source
		deno = tenRho_p(2,2)*tenRho_p(1,1) - tenRho_p(1,2)**2
		if ( deno < 1.0e-6 ) return 
		s = ( tenRho_p(2,2)*tenRho_p(1,3) - tenRho_p(1,2)*tenRho_p(2,3) ) / deno
		t = ( tenRho_p(1,1)*tenRho_p(2,3) - tenRho_p(1,2)*tenRho_p(1,3) ) / deno
		MirrorSourceXYZ(3) = - SourceXYZ(3)
		MirrorSourceXYZ(1) = 2 * SourceXYZ(3) * s + SourceXYZ(1)
		MirrorSourceXYZ(2) = 2 * SourceXYZ(3) * t + SourceXYZ(2)

		!	Compute Vp(nNode)
		do i = 1, nNode
			x = coor(1,i)
			y = coor(2,i)
			z = coor(3,i)
			if(dabs( x - SourceXYZ(1) ) < 1.e-8 .and. dabs( y - SourceXYZ(2) ) < 1.e-8 &
			   .and. dabs( z - SourceXYZ(3) ) < 1.e-8 ) then
				Vp(i) = 100000.0
               else
                    !!!!!
				call AnisoHalfSpace( current, tenRho_p, SourceXYZ(1), SourceXYZ(2), SourceXYZ(3), x, y, z, Vp(i) )
			end if
		end do

!----- Volume elements		
		do i = 1, Nelem

			m = - elem(5,i)

			if(m==0) m = 1		!*****************************

			tenRho = tenRho_temp(:,:,m)
			tenSigma = tenSigma_temp(:,:,m)
			tenSigma_s = tenSigma_temp(:,:,m) - tenSigma_p

			Ke = 0.0
			Pe = 0.0

			do j = 1, 4
				NT = elem(j,i)
				VolX(j) = coor(1,NT)
				VolY(j) = coor(2,NT)
				VolZ(j) = coor(3,NT)
				eVp(j) = Vp(NT)
			end do
			
			call GetVolume(coor(:,elem(1,i)),coor(:,elem(2,i)),coor(:,elem(3,i)),coor(:,elem(4,i)),Volume)

			call GetKe(tenSigma,tenSigma_s,VolX,VolY,VolZ,Volume,Ke,Pe)
			P = matmul(Pe,eVp)

			do mm=1,4
				node1 = elem(mm,i)
				do nn=1,4
					node2 = elem(nn,i)
					if( node1 .ge. node2 ) then
						call AddInSK( Nnode, rowWidth, node1, node2, Ke(mm,nn), SA, IJA )
					end if
				end do
				SigmaP(node1) = SigmaP(node1) + P(mm)
			end do

		end do

!----- Face elements
		do i = 1, Nface

			attr = - face(4,i)
			if( attr.ne.0 ) then

				do j = 1, 3
					NT = face(j,i)
					SurX(j) = coor(1,NT)
					SurY(j) = coor(2,NT)
					SurZ(j) = coor(3,NT)
					eVp2(j) = Vp(NT)
				end do
				
				layerID = mod(Attr,10)
				if( layerID > nR ) then
					write(*,*)	"facet",i,": layerID > nR, check!"
					stop
				end if
				
				tenRho = tenRho_temp(:,:,layerID)

				faceID = Attr / 10			! 从attr中提取“面”信息，进而可以对不同的面分别进行边界积分
				if( faceID > 3 ) then
					write(*,*) "facet",i,": wrong faceID, check!"
					stop
				end if

				!----- front and back surface	! 前后面属性1
				if( faceID == 1 ) then
					xx = SurX(1)
					yy = ( SurY(1) + SurY(2) + SurY(3) ) / 3.
					zz = ( SurZ(1) + SurZ(2) + SurZ(3) ) / 3.
					B = CptB(tenRho,xx,yy,zz,SourceXYZ)	
					BB = CptB(tenRho,xx,yy,zz,MirrorSourceXYZ)	
					call GetArea(SurX,SurY,SurZ,area)
					q = ( dabs(xx-SourceXYZ(1)) / (B**1.5) + dabs(xx-MirrorSourceXYZ(1)) / (BB**1.5) ) / ( 1/dsqrt(B) + 1/dsqrt(BB) )
					coef = q * area / 12.0

					Bp = CptB(tenRho_p,xx,yy,zz,SourceXYZ)
					BBp = CptB(tenRho_p,xx,yy,zz,MirrorSourceXYZ)
					qp = ( dabs(xx-SourceXYZ(1)) / (Bp**1.5) + dabs(xx-MirrorSourceXYZ(1)) / (BBp**1.5) ) / ( 1/dsqrt(Bp) + 1/dsqrt(BBp) )
					coef1 = (q-qp) * area / 12.0 
				end if

				!----- left and right surface	! 左右面属性2
				if( faceID == 2 ) then
					xx = ( SurX(1) + SurX(2) + SurX(3) ) / 3.
					yy = SurY(1) 
					zz = ( SurZ(1) + SurZ(2) + SurZ(3) ) / 3.
					B = CptB(tenRho,xx,yy,zz,SourceXYZ)	
					BB = CptB(tenRho,xx,yy,zz,MirrorSourceXYZ)	
					call GetArea(SurX,SurY,SurZ,area)
					q = ( dabs(yy-SourceXYZ(2)) / (B**1.5) + dabs(yy-SourceXYZ(2)) / (BB**1.5) ) / ( 1/dsqrt(B) + 1/dsqrt(BB) )
					coef = q * area / 12.0

					Bp = CptB(tenRho_p,xx,yy,zz,SourceXYZ)
					BBp = CptB(tenRho_p,xx,yy,zz,MirrorSourceXYZ)
					qp = ( dabs(yy-SourceXYZ(2)) / (Bp**1.5) + dabs(yy-SourceXYZ(2)) / (BBp**1.5) ) / ( 1/dsqrt(Bp) + 1/dsqrt(BBp) )
					coef1 = (q-qp) * area / 12.0
				end if

				!----- bottom surface	! 底面属性3
				if( faceID == 3 ) then
					xx = ( SurX(1) + SurX(2) + SurX(3) ) / 3.
					yy = ( SurY(1) + SurY(2) + SurY(3) ) / 3. 
					zz = SurZ(1)
					B = CptB(tenRho,xx,yy,zz,SourceXYZ)	
					BB = CptB(tenRho,xx,yy,zz,MirrorSourceXYZ)	
					call GetArea(SurX,SurY,SurZ,area)
					q = ( dabs(zz-SourceXYZ(3)) / (B**1.5) + dabs(zz-SourceXYZ(3)) / (BB**1.5) ) / ( 1/dsqrt(B) + 1/dsqrt(BB) )
					coef = q * area / 12.0

					Bp = CptB(tenRho_p,xx,yy,zz,SourceXYZ)
					BBp = CptB(tenRho_p,xx,yy,zz,MirrorSourceXYZ)
					qp = ( dabs(zz-SourceXYZ(3)) / (Bp**1.5) + dabs(zz-SourceXYZ(3)) / (BBp**1.5) ) / ( 1/dsqrt(Bp) + 1/dsqrt(BBp) )
					coef1 = (q-qp) * area / 12.0
				end if
				
				
				Pe2 = 1.
				do j = 1, 3
					Pe2(j,j) = 2.
				end do
				Pe2 = coef1 * Pe2
				P2 = matmul(Pe2,eVp2)

				do mm=1,3
					node1 = face(mm,i)
					do nn=1,3
						node2 = face(nn,i)
						if( node1 .gt. node2 ) then
							call AddInSK( Nnode, rowWidth, node1, node2, coef, SA, IJA )
						end if
						if( node1 .eq. node2 ) then
							coef2 = 2 * coef
							call AddInSK( Nnode, rowWidth, node1, node2, coef2, SA, IJA )
						end if
					end do
					SigmaP(node1) = SigmaP(node1) + P2(mm)
				end do

			end if

		end do

	
	end subroutine GetSigmaK
	
!-------------------------------------------------------------------------------------------------------------------------
!	Add an entry into 'SA' ------
!--------------------------------
	subroutine AddInSK( Nnode, rowWidth, rowOfA, columnOfA, num, SA, IJA )
		integer :: rowOfA, columnOfA, rowWidth, Nnode, IJA(*)
		real(kind=8) :: num, SA(*)
		integer :: i, m, existColumn, lctBegin, lctEnd
		
		if( rowOfA == columnOfA ) then
			sa(rowOfA) = sa(rowOfA) + num
		else 
			lctBegin = nNode + 1 + rowWidth * ( rowOfA - 1 ) + 1
			lctEnd = lctBegin + rowWidth - 1
			do i = lctBegin, lctEnd
				existColumn = IJA(i)
				if( existColumn == 0 ) then
					IJA(i) = columnOfA
					SA(i) = num
					exit
				else if( existColumn == columnOfA ) then
					SA(i) = SA(i) + num
					exit
				end if
			end do
			if( i == ( lctEnd + 1 ) ) then
				write(*,*)	"Error: subroutine AddInSK, rowWidth too small !"
				stop
			end if

		end if

	end subroutine AddInSK

end module ANALYSE 


