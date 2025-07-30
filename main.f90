
!----------------------------------------------------------------------------------------------------------------------------
!										USGDC < VERSION 4.1>
!	PROGRAM DESCRIPTION:
!		- THREE-D DC ANISOTROPIC RESISTIVITY MODELLING USING UNSTRUCTURED GRID.
!		- 4 INPUT FILES: RES.TXT, *.1.NODE, *.1.FACE, *.1.ELE
!		- THE LOCATION OF THE SOURCE SHOULD BE SPECIFIED IN SUBROUTINE 'GETSIGMAK()', WHICH IS IN FILE 'ANALYSE.F90'.
!		- Interested area should be revised in 'ResProc.f90'
!
!		- THE SOLVING DOMAIN IS A HALF-SPHERE.	
!
!		- HIGHTLIGHT: DO NOT NEED THE RELATIONSHIP TABLE ANY MORE!	
!
!		- Grid is produced by Tetgen (Si Hang). The command line should be somehow like this: tetgen.exe -pq1.2A, 
!		  that is to say, attributes of elements should be specified. The attribute of background should always be -1.
!
!		- PLATFORM: WINDOWS XP
!															
!-----------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------
!	在以上基础上修改了此程序：
!		- 求解域换为长方体形状；
!
!	注意，此程序用于测井，源在地下一定深度，故在全空间中计算。在tetgen的poly文件中地表的属性不可以再为0了，要参与边界积分。
!-----------------------------------------------------------------------------------------------------------------------------
!	logging1r2在LoggingUstrVs基础上修改，主要将pole-pole装置改为dipole-dipole装置。
!-----------------------------------------------------------------------------------------------------------------------------
!	longging1r3在logging1r2基础上修改，1r2中将地表当无穷边界（混合边界条件）处理，1r3中采用Neumann边界条件。
!	在生产poly文件的时候需要将地表的属性从-31改为0
!-----------------------------------------------------------------------------------------------------------------------------
!	logging1r4在logging1r3基础上修改，1r3中在地表采用纽曼边界条件，那么推导边界条件的时候需考虑镜像的作用。事实上1r3相当于
!	忽略了镜像的微弱作用，在1r4中，我们考虑镜像。
!
!-----------------------------------------------------------------------------------------------------------------------------

program main
	use analyse		!----- Assemble the system matrix
	use ResProc		!----- Rusult processing
	use Solver
    
    use ran_mod

	implicit none

	integer, parameter :: rowWidth = 40, NMAX = 100000000
	integer :: Nnode, Nelem, Nface, nResType, iter, MEMORY, mi, tempmin
	integer :: nR, StartTime, EndTime, INT4, i, j, k, ii, sourcePt, nElectrodes
	real(kind=8) :: current, SourceXYZ(3), coef1, coef2, coef3, ranr, rang, ranr2

	integer, allocatable :: elem(:,:), face(:,:), IJA(:)
	real(kind=8), allocatable :: rho(:,:), tRho(:,:,:), coor(:,:), SigmaP(:), Vp(:), x(:), SA(:)
    
    real(kind=8) :: kamn(200),beta(200),res(200),eulara,eularb,eularc
    real(kind=8) , parameter :: myPI = 3.141592653589793

	integer :: StartPt, EndPt, pt, countmin(200), countrand(100)
	real(kind=8), allocatable :: potential1(:,:),potential2(:,:),potential3(:,:),potential4(:,:)
	

	write(*,"(' PROGRAM ACCESSING...',/)") 
	call Declaration()
	
!----------------------------------------------------------------------------------------------------------------------
!	reading data	-----
!------------------------
        
	!----- current and resistivity -----
	Open(unit=1,file='res.txt')
	write(*,*)	"Reading resistivity..."
	read(1,*) current
	read(1,*)	nR, nResType
	if( nResType == 1 ) then	! Resistivities expressed by three principal resistivities and three Euler angles
		allocate( rho(6,nR) )
		allocate( tRho(3,3,nR) )
		do i = 1, 2
			read(1,*) rho(:,i)
			call pR2tR(rho(:,i),tRho(:,:,i))
		end do
		do i = 3, nR
			rho(:,i) = rho(:,1)
			call pR2tR(rho(:,i),tRho(:,:,i))
		end do
	else if( nResType == 2 ) then	! Resistivities expressed by tensors
		allocate( tRho(3,3,nR) )
		do i = 1, nR
			read(1,*) ( tRho(j,:,i), j=1,3 )
		end do
	else 
		write(*,*) "Wrong resistivity type! Check!"
		stop
	end if
	write(*,*)	"Reading resistivity finished."
    close(1)
    
	!----- coordinates of nodes -----
	open(unit=2,file='model12.1.node')
	write(*,*)	"Reading *.node file..."
	read(2,*) Nnode
	allocate( coor(3,Nnode), SigmaP(Nnode), Vp(nNode), x(Nnode) )
	do i = 1, Nnode
		read(2,*) j, coor(:,i)
	end do 
	write(*,*)	"Reading *.node file finished."
	write(*,*)
	close(2)

	!----- codes of elements ------
	open(unit=3,file='model12.1.ele')
	write(*,*)	"Reading *.ele file..."
	read(3,*) Nelem
	allocate( elem(5,Nelem) )
	do i = 1, Nelem
		read(3,*) j, elem(:,i)
	end do 
	write(*,*)	"Reading *.ele file finished."
	write(*,*)
	close(3)

	!----- boundary surface -------
	open(unit=4,file='model12.1.face')
	write(*,*)	"Reading *.face file..."
	read(4,*) Nface
	allocate( face(4,Nface) )
	do i = 1, Nface
		read(4,*) j, face(:,i)
	end do 
	write(*,*)	"Reading *.face file finished."
	write(*,*)
	close(4)

	write(*,"(' Nodes:',I15)") Nnode
	write(*,"(' Elements:',I15,/)") Nelem

	MEMORY = 5*nNode*8 + 5*nElem*4 + 4*nFace*4 + NMax*8 + NMax*4
	write(*,*)	"Memory allocating..."
	write(*,*)	" ", MEMORY, "bytes will be allocated!"
	allocate(SA(NMax),IJA(NMax))
	write(*,*)	"Finished."
	write(*,*)

call system_clock(startTime,INT4)

    call random_seed()
    call random_number (ranr)
    write(*,*) nR, ranr
    
	sourcePt = 1348    ! now 142... cqtx:25 fanchang:9 cqtx no tunnel:17
	write(222,"('Source coordinate:',3F20.10,'--------------------------------')")	coor(:,sourcePt)
    
    countmin = 0
    countrand = 0
        
    do i = 2, 198, 2
        kamn(i/2) = 2*myPI*i*(i+2)
        !write(*,*) i, i+2, kamn(i/2)
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get beta of the background
    

    !    do mi = 1, 1
            
        eulara = 0
        eularb = 0
        eularc = 0       
        rho(1,1) = 500
        rho(2,1) = 500
        rho(3,1) = 500
        rho(4,1) = eulara
        rho(5,1) = eularb
        rho(6,1) = eularc
        call pR2tR(rho(:,1),tRho(:,:,1)) 
            
            
            
            
            
            
            
            
            
            
            
        !!!!! 围岩各向异性
        !call random_number (ranr)
        !eulara = myPI*2*ranr
        !call random_number (ranr)
        !eularb = myPI*ranr
        !call random_number (ranr)
        !eularc = myPI*2*ranr
        !call random_number (ranr)
        !ranr = 1 + ranr*2
        !rho(1,1) = 500*ranr
        !rho(2,1) = 500*ranr
        !rho(3,1) = 500/ranr
        !rho(4,1) = eulara
        !rho(5,1) = eularb
        !rho(6,1) = eularc
        !call pR2tR(rho(:,1),tRho(:,:,1))
        
        
        
        
 !	Assemble system matrix	-----
!--------------------------------
    
	call InitSK( Nnode, NMAX, IJA, SA, rowWidth )
	call GetSigmaK(current,nR,tRho,Nnode,coor,Nelem,elem,Nface,face,Vp,SigmaP,SA,IJA,rowWidth,sourcePt,SourceXYZ)


    !   Solving A x = b  -----
!-------------------------  
    
	SigmaP = - SigmaP
	call WILL_SSORCG( nNode, NMax, 2000, SA, IJA, sigmaP, x, iter )
    
    call betaoutput(current, NNode, coor,  sourcePt, Vp, x, beta)
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
do mi = 1, 5
    
        !!!!!! gexiangyixing
        !call random_number (ranr)
        !eulara = myPI*2*ranr
        !call random_number (ranr)
        !eularb = myPI*ranr
        !call random_number (ranr)
        !eularc = myPI*2*ranr
        
        !!!!!! gexiangtongxing
        eulara = 0
        eularb = 0
        eularc = 0
        
        do i = 3, nR
            !call random_number (ranr)
           ! rho(:,i) = ranr*480+20
           ! rho(:,i) = 20
            !!!!!! yixing
            !call random_number (ranr2)
            !ranr2 = 1 + ranr2*2
            !rho(1,i) = (ranr*480+20)*ranr2
            !rho(2,i) = (ranr*480+20)*ranr2
            !rho(3,i) = (ranr*480+20)/ranr2
            
           ! rho(4,i) = eulara
            !rho(5,i) = eularb
            !rho(6,i) = eularc
            
            
            
             
        rho(1,i) = 20
        rho(2,i) = 20
        rho(3,i) =20
        rho(4,i) = eulara
        rho(5,i) = eularb
        rho(6,i) = eularc
               
	        call pR2tR(rho(:,i),tRho(:,:,i))
        end do

	call InitSK( Nnode, NMAX, IJA, SA, rowWidth )
	call GetSigmaK(current,nR,tRho,Nnode,coor,Nelem,elem,Nface,face,Vp,SigmaP,SA,IJA,rowWidth,sourcePt,SourceXYZ)
	SigmaP = - SigmaP
        
    write(*,*) '****************************'
    write(*,*) mi
    write(*,*) '****************************'
        
	call WILL_SSORCG( nNode, NMax, 2000, SA, IJA, sigmaP, x, iter )
    
    ! Gaussian distribution
    !do i = 1, nNode
    !    rang = normal(0.001D0, 0.0003D0)
    !    countrand(rang*1000) = countrand(rang*1000) + 1
    !    x(i) = x(i)*(0.999+rang)
    !    !write(*,*) 0.95+rang
    !end do
    
    call cqtcoutput(current, NNode, coor,  sourcePt,Vp, x)
    call cqres(current, NNode, coor,  sourcePt,Vp, x, beta, kamn, res)
    
    tempmin = 1
    do i = 2, 99
        if (res(i) < res(tempmin)) then
            tempmin = i
        end if
    end do
    countmin(tempmin) = countmin(tempmin) + 1
    
end do

    open(unit=2,file="count.txt")
    do i = 1, 99
        write(2,*) i*2+1,countmin(i)
    end do
    close(2)
    
    open(unit=2,file="rand.txt")
    do i = 1, 50
        write(2,*) i,countrand(i)
    end do
    close(2)
    
	call system_clock(endTime)
	call TimeConsumed(startTime, endTime, INT4, iter)

end program main

!----------------------------------------------------------------------------------------------------------------------
!	Initialize SA(*) & IJA(*) -----
!----------------------------------
subroutine InitSK(Nnode,NMAX,IJA,SA,rowWidth)
	implicit none
	integer :: Nnode, NMAX, rowWidth, IJA(NMAX)
	real(kind=8) :: SA(NMAX)
	integer :: i

	IJA = 0
	SA = 0.0

	ija(1) = Nnode+2
	do i = 2, Nnode+1
		ija(i) = ija(i-1) + rowWidth
	end do 

end subroutine InitSK
    