!  beam.f90 
!  Clayton Smith
!  26 Mar 2018
!  http://claytonsmith.us
!  I/O and stiffness matrix code borrowed, with permission, from Kurt Gramoll, Ph.D., Univ. of Okla., 2018.
    
    ! Requires input file as part of the command string (assumes .txt)
    ! Format   
    !  line 1: number of elements
    !  line for each node (elements plus 1)
    !       node#, constraint (0=non, 1=disp, 2=rot, 3=both), loc, force, moment
    !  line for each element
    !       element#, left node#, right node#, E, I
    
! Sample, cantilever beam, load at end, 2 elements
!   2
!   1 3 0.0 0.0 0.0
!   2 0 1.0 0.0 0.0
!   3 0 2.0 1.0 0.0
!   1 1 2 100 1.0e-2
!   2 2 3 100 1.0e-2
	
    
    
    program beam
    
    implicit none
    
    !!  NOT USED, example of set dimension for array (old Fortran method)
	!!  assuming maximum nodes=50, max elements=50
    !!  real*8 :: gStiff(1:100, 1:101)  ! global stiffness matrix with RHS, size = numNodes*2, array start at 1
    !!  real*8 :: gStiff(100, 101)  ! global stiffness matrix with RHS, size = numNodes*2, array start at 1
    
    real*8, allocatable, dimension(:,:) :: gStiff, gStiffOrig   !! dynamic array method
    real*8, DIMENSION(1:4,1:4) :: elemStiff                     !! one for each element, old syntax
    
    real*8, allocatable, dimension(:,:) :: matIn, matL, matLT
    real*8, allocatable, dimension(:) :: y, x
    
    real*8, allocatable, dimension(:) :: nodeLocX, nodeLocY     !! location of node
    real*8, allocatable, dimension(:) :: load, disp, moment     !! node force (RHS) and displacement
    
    real*8, allocatable, dimension(:) :: elemLen, stress, angRad    !! element length, stress, orientation angle (in radians)
    real*8, allocatable, dimension(:) :: youngMod, area, momInertia            !! element modulus and area

    real*8 :: temp, temp1, sum          !! real = real*4 (single), real*8 (double, 64 bits), real*16 (quad)
    real*8 :: cosAng, sinAng            !! cosine and sine of each element bar
    real*8 :: elapsedTime, t1, t2       !! used for calculation time
    real*8 :: stiffnessMultiplier       !! multiplier in front of element stiffness matrix : EI/L^3
    
    integer, allocatable, dimension(:) :: nodeRest  !! type of restriction, 0=non, 1=disp x, 2=dips y, 3=both
    integer, allocatable, dimension(:) :: leftNode, rightNode    !! node order, left to right
    
    integer :: i, j, k, n, m, iRow, iCol, iElem, iNode   !! counters
    integer :: numNode, numElem, dof    !! node, element and total degree of freedom 
    
    integer :: nSize
    
    character (len=100) :: strBuffer100  !! string that holds command line argument
    character (len=30) :: fileName, fileInName, fileOutName  ! file name (no ext), file name in and out w/ ext


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Get command line input parameters, call for each item
    call GETARG(1, strBuffer100)   ! special function to read in argurments 
    read (strBuffer100,*) fileName   ! reading in file name (no extension)
    fileInName = TRIM(fileName) // ".txt"        ! file in name with extension
    fileOutName = TRIM(fileName) // "_out.txt"    ! file out name with extension
    
    !! write to console that program (1st *) has started, free format (2nd *)
    write (*,*) 'beam program start - reading in data from: ', fileInName
    
    CALL CPU_TIME(t1)     ! set timer, used to determine total time program runs
    
    !! Read in data file    ========================================================================
    open (UNIT=5, FILE=fileInName)     !  file number 5 standard for input, could use 'open (5, fileInName)'
    read (5,*), numElem        !  Just one line, needed to figure DOF (used for matrix setup
    !!  need to allocate memory for dynamic matrices, needed for load and disp vector size before reading
    dof = 2*(numElem + 1)
	allocate (gStiff(dof, dof))       !!  used when reducing global stiffness matrix to solvable set of equations (b.c. applied)
	allocate (gStiffOrig(dof, dof))
	allocate (load(dof))         ! external loads at each node, b.c. (maybe unknown)
    allocate (moment(dof))       ! external point moment at each node, b.c.
	allocate (disp(dof))         ! unknown node displacements (some known, b.c.)
    
	allocate (nodeLocX(numNode))
	allocate (nodeLocY(numNode))
	allocate (nodeRest(numNode))  ! node restriction (i.e. B.C.), 0 = none, 1 = x, 2 = y, 3 = both
    
	allocate (leftNode(numElem))  ! node 1 in local numbering
	allocate (rightNode(numElem)) ! node 2 in local numbering
	allocate (youngMod(numElem))
    allocate (momInertia(numElem))
	allocate (area(numElem))
    
	allocate (angRad(numElem))    ! calculate from nodeLocX and nodeLocY
	allocate (elemLen(numElem))   ! calculate from nodeLocX and nodeLocY
	allocate (stress(numElem))    ! to be found
    
    allocate (matIn(dof, dof))  ! temporary array for cholesky        -> matIn = gStiff
    allocate (matL(dof, dof))   ! temporary array for cholesky decomp -> Lower Triangular matrix
    allocate (matLT(dof, dof))  ! temporary array for cholesky decomp -> Lower Triangular matrix transposed
    allocate (y(dof))           ! temporary array for cholesky decomp -> L y = b
    allocate (x(dof))           ! temporary array for cholesky decomp -> LT x = y

    leftNode ( : ) = 0
    rightNode ( : ) = 0
    numNode = (numElem + 1)     ! number of nodes is number of elements + 1, because of overlap in 2D.
    
    !! Finish reading in data now that matrices sizes are set
    read (5,*) (n, nodeRest(i), nodeLocX(i), load(2*(i-1)+1), load(2*(i-1)+2), i=1, numNode)
    read (5,*) (n, leftNode(i), rightNode(i), youngMod(i), momInertia(i), i=1, numElem)
    
    !! Write out copy of inputted data     ============================================================
    open (UNIT=6, FILE=fileOutName)  !! file number 6 standard for output to file
    write (6,*) 'Number of elements: ', numElem
    write (6,*)
    write (6,*)  "Nodes"
    write (6,62)  (i, nodeRest(i), nodeLocX(i), load(2*(i-1)+1), load(2*(i-1)+2), i=1, numNode)
62  format ("  Pt  Fixed    X           Force       Moment" / (2I5, 3E12.4))    
    write (6,*)
    write (6,*)  "Elements"
    do i = 1, numElem
        elemLen(i) = (nodeLocX(rightNode(i)) - nodeLocX(leftNode(i)))
    end do
    write (6,63) (i, leftNode(i), rightNode(i), youngMod(i), momInertia(i), elemLen(i), i=1, numElem)
63  format ("  elem   i    j    Modulus     I           Length" / (3I5, 3E12.4))  
    
    !!debugging load input
    !write (6, 72) ( load(2*(i-1)+1), i=1, numNode) !load
    !write (6, 72) ( load(2*(i-1)+2), i=1, numNode) !moment
    !write (*,*) (load(i), i=1, dof) !the whole array
    
!   Initialize global stiffness matrix =================================================
    gStiff ( : , : ) = 0.0  ! or you can loop on matrix (slower)
    
 !  Construct element stiffness   ======================================================
    
    do iElem = 1, numElem
        
        stiffnessMultiplier = 0.0 !initialize the ((E * I) / (L**3)) term
        
        ! atan or datan (real*8) returns radians
        temp = nodeLocY(rightNode(iElem)) - nodeLocY(leftNode(iElem))
        temp1 = nodeLocX(rightNode(iElem)) - nodeLocX(leftNode(iElem))
        angRad(iElem) = datan( temp / temp1 )  ! put in array, use again for stress
        
        cosAng = cos(angRad(iElem))
        sinAng = sin(angRad(iElem))
        temp = area(iElem) * youngMod(iElem) / elemLen(iElem)
        
        stiffnessMultiplier = ((youngMod(iElem) * momInertia(iElem)) / (elemLen(iElem)**3))
        !! Debugging ===================================
        !write(*,69) (stiffnessMultiplier)
69      format(/"Stiffness Multiplier"  / (F9.3))
        
        elemStiff(1,1) = stiffnessMultiplier * 12
        elemStiff(1,2) = stiffnessMultiplier * (6 * elemLen(iElem))
        elemStiff(1,3) = stiffnessMultiplier * -12
        elemStiff(1,4) = stiffnessMultiplier * (6 * elemLen(iElem))
        
        elemStiff(2,1) = elemStiff(1,2)
        elemStiff(2,2) = stiffnessMultiplier * (4 * (elemLen(iElem))**2)
        elemStiff(2,3) = -elemStiff(2,1)
        elemStiff(2,4) = (2 * (elemLen(iElem)**2))    
        
        elemStiff(3,1) = -elemStiff(1,1)
        elemStiff(3,2) = -elemStiff(1,2)
        elemStiff(3,3) = -elemStiff(1,3)
        elemStiff(3,4) = -elemStiff(1,4)        
        
        elemStiff(4,1) = elemStiff(2,1)
        elemStiff(4,2) = elemStiff(2,4)
        elemStiff(4,3) = elemStiff(2,3)
        elemStiff(4,4) = elemStiff(2,2)        
        
        
        !write (*,71) ( (elemStiff(n, i), i=1, 4), n=1, 4)
71      format (/"Element Stiff"  / (4F9.3))
       
        ! put element stiffness into global matrix
        i = leftNode(iElem)
        j = rightNode(iElem)
        do n = 1, 4
            do k = 1, 4
                !! do each pair differently
                if (n < 3) iRow = 2*(i-1) + n
                if (n >= 3) iRow = 2*(j-1) + n-2
                if (k < 3) iCol = 2*(i-1) + k
                if (k >= 3) iCol = 2*(j-1) + k-2
                gStiff(iRow, iCol) = gStiff(iRow, iCol) + elemStiff(n, k)
            end do
        end do
        
    end do
    
    gStiffOrig = gStiff  !! save orginal for later to get all forces at nodes
    
    !! testing output, small matix only
    !write (6, 65) ( (gStiffOrig(n, i), i=1, dof), n=1, dof)
65  format (/"Stiffness"  /  (8E9.2) )  
   
    
!   Apply displacement boundary conditions   ================================================
!       set fixed node direction to zero (row and column) except when i=j, then 1.0
    do iNode = 1, numNode
        ! ZERO X deflection - set row and column to zero except diagonal term to 1
        if (nodeRest(iNode) == 1 .OR. nodeRest(iNode) == 3) then   ! x-dir fixed
            iRow = 2*(iNode - 1) + 1
            do iCol = 1, dof
                if (iCol /= iRow) then 
                    gStiff(iRow, iCol) = 0.0
                    gStiff(iCol, iRow) = 0.0
                else 
                    gStiff(iRow, iCol) = 1.0  !  can be anything except 0
                end if
            end do
            load(2*(iNode-1)+1) = 0.0  ! load is also set to zero (when solving, gives defl = 0)
        end if

        ! ZERO Y deflection - set row and column to zero except diagonal term to 1
        if (nodeRest(iNode) == 2 .OR. nodeRest(iNode) == 3) then   ! y-dir fixed
            iRow = 2*(iNode - 1) + 2
            do iCol = 1, dof
                if (iCol /= iRow) then 
                    gStiff(iRow, iCol) = 0.0
                    gStiff(iCol, iRow) = 0.0
                else 
                    gStiff(iRow, iCol) = 1.0  !  can be anything except 0
                end if
            end do
            load(2*(iNode-1)+2) = 0.0
        end if
    end do
    
    
    !!  Testing output    
    !write (6, 65) ( (gStiff(n, i), i=1, dof), n=1, dof)
    !write (6, 72) ( load(2*(i-1)+1), i=1, numNode)
    !write (6, 73) ( load(2*(i-1)+2), i=1, numNode)
    !write (6,*) (load(i), i=1, dof)
72  format (/"force"  /  (8F8.0))
73  format (/"moment" /  (8F8.0))
    
    
    
    !==================   Solve  ======================================
    !  Find all displacements
  
    !==================   Initialize ============
    matIn = gStiff
    nSize = dof
    matL(:,:) = 0.0
    matLT(:,:) = 0.0
    disp(:) = 0.0
    y(:) = 0.0
    x(:) = 0.0
    

    !==================   decompose  ============
    
    do k = 1, nSize !outside loop
        do i = 1, k !secondary loop
            if (i .ne. k) then
                sum = 0.0 !clear the sum var
                do j = 1, i - 1
                    sum = sum + (matL(i,j) * matL(k,j)) !performs the summation of this factor
                end do
                matL(k,i) = (matIn(k,i) - sum) / matL(i,i) !subtracts the previous factor to the original matrix, and devides it by matL(i,i)
                
            else if (i .eq. k) then
                sum = 0.0 !clear the sum var
                do j = 1, i - 1
                    sum = sum + (matL(k,j) * matL(k,j)) !performs the summation of this factor
                end do
                
                matL(k,k) = sqrt(matIn(k,k) - sum) !subtracts the previous factor to the original matrix, takes the square root
            end if
        end do
    end do
    
    !debugging ===========================================================================================
    !print *, "matL"
    !do i =1, nSize ! print out the cholesky lower triangular matrix, matL(*,*)
    !    write (*, '(*(F7.3))' ) ( matL(i,j), j=1, nSize )
    !end do
    
    !=================   sub  ===================
    
    matLT = transpose(matL) !transpose this matrix
    
    do k = 1, nSize !forward sub
        sum = 0.0 !clear the var
        do i = 1, nSize
            sum = sum + (matL(k,i) * y(i)) !sum factor
        end do
        y(k) = (load(k) - sum) / matL(k,k)
    end do
    
    do k = nSize, 1, -1 !back sub
        sum = 0.0 !clear the var
        do i = nSize, 1, -1
            sum = sum + (matLT(k,i) * x(i)) !sum factor
        end do
        x(k) = (y(k) - sum) / matLT(k,k)
    end do
    
    disp = x !write to the var Gramoll uses.
    
  
  
    !================================================  

    
    ! debugging console output
    !write (*, 74) (disp(2*(i-1) + 1), i=1, numNode)

74  format (/"displacement"  /  (8F8.3))
    
    !! get all forces (known and unknown), most will be zero
    do i = 1, dof
        load(i) = 0.0  
        do j = 1, dof
            load(i) = load(i) + gStiffOrig(i, j) * disp(j)
        end do
    end do  
    
    
 !!  Find stress in each member
 !   do iElem = 1, numElem
 !       cosAng = cos(angRad(iElem))
 !       sinAng = sin(angRad(iElem))
 !       
 !       i = leftNode(iElem)
 !       j = rightNode(iElem)
 !       temp = youngMod(iElem) / elemLen(iElem)
 !       stress(iElem) = temp * ( -cosAng*disp(2*(i-1)+1) - sinAng*disp(2*(i-1)+2) &
 !                           + cosAng*disp(2*(j-1)+1) + sinAng*disp(2*(j-1)+2) )
 !   end do
    
    
    CALL CPU_TIME(t2)
    !write (*,"('elasped time:', E14.8)") (t2-t1)


!   Results Output to file
    write (6,*)
    write (6,*)
    write (6,*) "----------- Results ---------------"
    write (6,*) "Results    "
    write (6,64) (i, nodeRest(i), nodeLocX(i), load(2*(i-1) + 1), disp(2*(i-1) + 1), load(2*(i-1) + 2), disp(2*(i-1) + 2), i=1, numNode) 
64  format (" Node  Constraint Location    Rec Load    Defl      Rec moment    Rotation" / (2I5, 5x, 5E12.4))
!64  format (" Node  Constraint Location    Rec Load    Defl      Rec moment    Rotation" / (2I5, 5x, 5E12.4))
    
    

!    write (6,*)
!    write (6,*)
!    write (6,*) "---------- Results ---------------"
!    write (6,*) "Global Node Displacement and Force"
!    write (6,64) (i, load(2*(i-1)+1), disp(2*(i-1)+1), load(2*(i-1)+2), disp(2*(i-1)+2), i = 1, numNode)
!64  format (" Node   X Load      X Disp      Y Load      Y Disp   "  /  (I5, 4E12.4))
!
!!   Results Output to file
!    write (6,*)
!    write (6,*) "Element Axial Load and Stress"
!    write (6,75) (iElem, stress(iElem)*area(iElem),  stress(iElem), iElem = 1, numElem)
!75  format ("  Element     Load         Stress   "  /  (I7, 1E14.4, 1E14.4))
!
!    
    write (*,*) 'beam End'
    
    end program beam

