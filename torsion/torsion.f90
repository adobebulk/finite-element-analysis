!  torsion.f90 
!  Clayton Smith
!  30 April 2018
!  http://claytonsmith.us
!  I/O and allocation code borrowed, with permission, from Kurt Gramoll, Ph.D., Univ. of Okla., 2018.
!  Takes in input for: ay, bx, by, cx, (number of elements along axis), shear Modulus, applied torque
!  For a 1/4 cross-section
!  Output: Total bend angle (phi), stresses at a, b, c.
    
    program torsion

    
    implicit none

    
    real*8, allocatable, dimension(:,:) :: gStiff, gStiffOrig   !! dynamic array method
    real*8, DIMENSION(1:3,1:3) :: elemStiff                     !! one for each element, old syntax
    real*8, allocatable, dimension(:, :, :, :) :: elemInfo      !! elemInfo(elem #, i loc, j loc, m loc)
    real*8, allocatable, dimension(:) :: load, oldLoad             !! load(=(node i, node j, node k values))
    real*8, allocatable, dimension(:) :: psiValues, torsionValues
    
    real*8, allocatable, dimension(:,:) :: matIn, matL, matLT
    real*8, allocatable, dimension(:) :: y, x
    
    real*8, allocatable, dimension(:) :: nodeX, nodeY     !! location of node
    real*8, allocatable, dimension(:) :: disp, moment    !! node force (RHS) and displacement
    
    real*8, allocatable, dimension(:) :: elemLen, stress, angRad, elemI, elemJ, elemM, gradMat, gradMatT   !! element length, stress, orientation angle (in radians)
    real*8, allocatable, dimension(:) :: youngMod, area, momInertia           !! element modulus and area
    real*8, allocatable, dimension(:) :: alphI, betI, gamI, alphJ, betJ, gamJ, alphM, betM, gamM
    real*8, allocatable, dimension(:) :: shearArray
    
    real*8 :: shearModulus
    real*8 :: temp, temp1, sum          !! real = real*4 (single), real*8 (double, 64 bits), real*16 (quad)
    real*8 :: cosAng, sinAng            !! cosine and sine of each element bar
    real*8 :: elapsedTime, t1, t2       !! used for calculation time
    real*8 :: stiffnessMultiplier       !! multiplier in front of element stiffness matrix : EI/L^3
    real*8 :: torque, torqueE
    real*8 :: gNum
    real*8 :: ay, bx, by, cx
    real*8 :: elemH, elemW
    real*8 :: fourA
    real*8 :: psiTot !total Psi
    real*8 :: angleTwist
    real*8 :: shearY
    real*8 :: shearX
    real*8 :: shearXY
    real*8 :: maxShear
    real*8 :: maxShearLoc
    !real*8 :: i, j, k, n, m, iRow, iCol, iElem, iNode  !! counters

    
    integer :: i, j, k, n, m, iRow, iCol, iElem, iNode  !! counters
    integer :: totNodes, totElements, dof                    !! node, element and total degree of freedom 
    integer :: nodeNum, elemNum !! Iterators, used in iteration for current element number or current node number
    integer :: nSize
    
    character*100 :: strBuffer100                       !! Command Line Arguments
    character*30 :: fileName, fileInName, fileOutName   !! file name (no ext), file name in and out w/ ext
    
    real*8 :: diagSlope = 0.0
    real*8 :: diagLength = 0.0
    real*8 :: diagSpacing = 0.0
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call CPU_TIME(t1) !begin timer
    
    call GETARG(1, strBuffer100) !read in arguments
    read (strBuffer100, *) ay
    call GETARG(2, strBuffer100)
    read (strBuffer100,*) bx
    call GETARG(3, strBuffer100)
    read (strBuffer100,*) by
    call GETARG(4, strBuffer100)
    read (strBuffer100,*) cx
    call GETARG(5, strBuffer100)
    read (strBuffer100,*) gNum
    call GETARG(6, strBuffer100)
    
    !allocate (shearModulus(20))
    
    read (strBuffer100,*) shearModulus
    call GETARG(7, strBuffer100)
    read (strBuffer100,*) torque

    
    write (*,*) 'torsion program Start - Reading in data from command line.'

    !write(*,*)'Number of Points:  ', totNodes, '; elements', totElements
    

    
    !! Allocate What we know
    totNodes = (gNum + 1)* (gNum + 1)
    totElements = (gNum * (gNum * 2))
    
    dof = totNodes
    
    write(*,*) (" ")
    write(*,*) ("---- Input ----------------")
    write(*,121) (ay)
    write(*,122) (bx)
    write(*,123) (by)
    write(*,124) (cx)
    write(*,125) (gNum)
    write(*,126) (shearModulus)
    write(*,127) (torque)
    write(*,128) (totElements)
    write(*,129) (totNodes)
    write(*,13)  !horizontal dashes. format is at bottom of code
    write(*,101) !blank line

    !write(*, 10002) (ay, bx, by, cx)
    !10002   format( ' A ', F9.3, 'Bx ', F9.3, 'By', F9.3, 'C', F9.3)    !<---- For whatever reason, ths won't compile.
101 format(/" ")
    
121 format(" A:    " (F9.3))
122 format(" Bx:   " (F9.3))
123 format(" By:   " (F9.3))
124 format(" C:    " (F9.3))    
125 format(" gNum  " (F9.3))   
126 format(" shearMod  " (E9.3))
127 format(" torque" (F9.3))
128 format(" Elements: " (i5))
129 format(" Nodes:    " (i5))    
    

    
    allocate (gStiff(dof, dof))
    gStiff(:, :)    = 0.0 !initialize gStiff
    elemStiff(:, :) = 0.0 !initialize elemStiff
    allocate (load(totElements)) !was dof
    allocate (oldLoad(dof))
    allocate (psiValues(totElements))
    allocate (torsionValues(totElements))
    allocate (shearArray(totElements))
    
    load (:) = 0.0 !initialize elemLoad
    psiValues(:) = 0.0
    psiTot = 0.0
    torqueE = 0.0
    angleTwist = 0.0
    
    allocate (nodeX(totNodes))
    allocate (nodeY(totNodes))
    
    allocate (elemI(totElements))
    allocate (elemJ(totElements))
    allocate (elemM(totElements))
    
    allocate (gradMat(totNodes))
    allocate (gradMatT(totNodes))
    
    allocate (alphI(totElements))
    allocate (betI(totElements))
    allocate (gamI(totElements))
    
    allocate (alphJ(totElements))
    allocate (betJ(totElements))
    allocate (gamJ(totElements))
    
    allocate (alphM(totElements))
    allocate (betM(totElements))
    allocate (gamM(totElements))
    
    allocate (area(totElements))
    
    allocate(matIn(dof, dof))
    allocate(matL(dof,dof))
    allocate(disp(dof))
    allocate(x(dof))
    allocate(y(dof))
    
    ! ***********************
    ! ***Grid Up the Nodes***
    ! (this indexing took forever to figure out)
    
    ! Diagonal Spacing Stuff
    diagSlope = by / bx
    diagLength = sqrt(bx**2 + by**2)
    diagSpacing = diagLength / gNum
    
    do i = 1, (gNum + 1)
        elemNum = (gNum + 2) * (i - 1) + 1
        nodeX(elemNum) = (bx / gNum) * (i - 1)
        nodeY(elemNum) = (by / gNum) * (i - 1)
    end do
    
    ! Upper Left Triangle
    elemH = ay / gNum
    elemW = bx / gNum
    
    do i = 0, (gNum - 1)
        do j = (i + 1), gNum
            nodeNum = ((gNum + 1) * i) + j + 1
            nodeX(nodeNum) = i * elemW
            nodeY(nodeNum) = nodeY(((gNum + 2) * i) + 1) + (j - i) * elemH
        end do
    end do
    
    ! Bottom Right Triangle
    elemH = by / gNum
    elemW = cx / gNum
    
    do i = 1, gNum
        do j = 0, (i - 1)
            nodeNum = (gNum + 1) * i + J + 1
            nodeX(nodeNum) = nodeX((gNum + 2) * j + 1) + (i - j) * elemW
            nodeY(nodeNum) = j * elemH
        end do
    end do
    
    !! Node Output Debugging, Only Use SMALL gNum
    !write (*,1) (i, nodeX(i), nodeY(i), i = 1, totNodes)
1   format(" pt     nodeX   nodeY" / (I2, 2F12.4))
    
    ! Find the Element Numbers
    
    !method: index, counting upward, as boxes from the left node number locations. Split the box into two triangles,
    !an upper triangle, and a lower triangle. The upper triangle is (k), the lower triangle is (k + 1).
    !the element node numbers are stored in elemI, elemJ, elemM, respectively, in a CCW fashion.
    
    !***************************
    !*** GRAPHIC OF TRIANGLES***
    !      _____
    !elemM | k / elemJ
    !      |  /
    !      | /
    !      |/
    !      elemI
    
    !        /| elemM
    !       / |
    !      /  |
    !     /k+1|
    !elemI-----elemJ
    !
    !***************************
    k = 0
    do i = 1, gNum
        do j = 0, (gNum - 1)
            
            !upper triangle
            k = k + 1
            elemI(k) = ((i - 1) * (gNum)) + i + j
            elemJ(k) = (i * gNum) + (3 + (i - 1)) + j
            elemM(k) = (((i - 1) * (gNum)) + i + j) + 1 !this is just elemI + 1, because it's just 1 up
            
            !lower triangle
            k = k + 1
            elemI(k) = ((i - 1) * (gNum)) + i + j
            elemJ(k) = (((i + 1) - 1) * (gNum)) + (i + 1) + j !this is the same as elemI, but one row over. So, just add (i+1) for i.
            elemM(k) = (i * gNum) + (3 + (i - 1)) + j !this is the same as the upper triangle elemJ
            
            !Dr. Gramoll wasn't kidding. This indexing was tough!
        end do
    end do
    
    !write (*,2) (k, elemJ(k), k = 1, totElements)
2   format("    k               elem Value" / (2F12.4))
    
    ! ***********************************
    ! Form the element stiffness matrices
    
    !******************************
    !*** Double Lookup Debugger ***
    !do k = 1, totElements
    !    write (*,*) (nodeX(elemI(k)))
    !    write (*, *) (nodeY(elemI(k)))
    !end do
    
    do k = 1, (totElements)

        !alphI(k) = (nodeY(elemJ(k)) * -nodeX(elemM(k))) - (nodeY(elemM(k)) * -nodeX(elemJ(k)))
        betI(k)  = -nodeX(elemJ(k)) + nodeX(elemM(k))
        gamI(k)  = nodeY(elemM(k)) - nodeY(elemJ(k))
        
        alphJ(k) = (nodeY(elemM(k)) * -nodeX(elemI(k))) - (nodeY(elemI(k)) * -nodeX(elemM(k)))
        betJ(k)  = -nodeX(elemM(k)) + nodeX(elemI(k))
        gamJ(k)  = nodeY(elemI(k)) - nodeY(elemM(k))
        
        !alphM(k) = (nodeY(elemI(k)) * -nodeX(elemJ(k))) - (nodeY(elemJ(k)) * -nodeX(elemI(k)))
        betM(k)  = -nodeX(elemI(k)) + nodeX(elemJ(k))
        gamM(k)  = nodeY(elemJ(k)) - nodeY(elemI(k))
        
        area(k) = 0.5 * ((nodeX(elemI(k)) * (nodeY(elemJ(k)) - nodeY(elemM(k)))) + (nodeX(elemJ(k)) * (nodeY(elemM(k)) - nodeY(elemI(k)))) + (nodeX(elemM(k)) * (nodeY(elemI(k)) - nodeY(elemJ(k)))))
        
        
        !write(*,*) (gamI(1))
        !write(*,*) (gamJ(1))
        !write(*,*) (gamM(1))
    end do
    
    !calculate 
    !       beta^T * beta + gamma^T * gamma
    !       _______________________________
    !                    4A
    ! Yields a 3x3 element stiffness matrix
    do k = 1, totElements
        fourA = 1 / (4 * area(k))
        elemStiff(1,1) = (fourA) * ((betI(k) * betI(k)) + (gamI(k) * gamI(k)))
        elemStiff(1,2) = (fourA) * ((betI(k) * betJ(k)) + (gamI(k) * gamJ(k)))
        elemStiff(1,3) = (fourA) * ((betI(k) * betM(k)) + (gamI(k) * gamM(k)))
    
        elemStiff(2,1) = (fourA) * ((betJ(k) * betI(k)) + (gamJ(k) * gamI(k)))
        elemStiff(2,2) = (fourA) * ((betJ(k) * betJ(k)) + (gamJ(k) * gamJ(k)))
        elemStiff(2,3) = (fourA) * ((betJ(k) * betM(k)) + (gamJ(k) * gamM(k)))
    
        elemStiff(3,1) = (fourA) * ((betM(k) * betI(k)) + (gamM(k) * gamI(k)))
        elemStiff(3,2) = (fourA) * ((betM(k) * betJ(k)) + (gamM(k) * gamJ(k)))
        elemStiff(3,3) = (fourA) * ((betM(k) * betM(k)) + (gamM(k) * gamM(k)))
        
        gStiff(elemI(k), elemI(k)) = gStiff(elemI(k), elemI(k)) + elemStiff(1,1)
        gStiff(elemI(k), elemJ(k)) = gStiff(elemI(k), elemJ(k)) + elemStiff(1,2)
        gStiff(elemI(k), elemM(k)) = gStiff(elemI(k), elemM(k)) + elemStiff(1,3)
        
        gStiff(elemJ(k), elemI(k)) = gStiff(elemJ(k), elemI(k)) + elemStiff(2,1)
        gStiff(elemJ(k), elemJ(k)) = gStiff(elemJ(k), elemJ(k)) + elemStiff(2,2)
        gStiff(elemJ(k), elemM(k)) = gStiff(elemJ(k), elemM(k)) + elemStiff(2,3)
        
        gStiff(elemM(k), elemI(k)) = gStiff(elemM(k), elemI(k)) + elemStiff(3,1)
        gStiff(elemM(k), elemJ(k)) = gStiff(elemM(k), elemJ(k)) + elemStiff(3,2)
        gStiff(elemM(k), elemM(k)) = gStiff(elemM(k), elemM(k)) + elemStiff(3,3)

        
        !write(*,3) ( (elemStiff(n,i), i = 1, 3), n = 1, 3)
3       format (/"elemStiff "   / (3F9.3))
        

    end do
    !Debug Stiffness Insert (not full dof bounds, only to 10).
    !write (*,4) ( (gStiff(n, i), i = 1, 10), n = 1, 10)
4   format (/"Stiffness"    / (10F9.2) )
    
    ! *********
    ! Load B.C.
    do k = 1, totElements
        load(elemI(k)) = load(elemI(k)) + 1
        load(elemJ(k)) = load(elemJ(k)) + 1
        load(elemM(k)) = load(elemM(k)) + 1
    end do
    
    do k = 1, totElements
        load(k) = load(k) * 0.01
    end do
    

    
    !oldLoad = load  
    
66  format(/"load(k) "/ (E9.3))

    do k = 1, dof
        load(k) = ((2 * shearModulus * area(k)) / 3) * load(k)
    end do

    !write(*,66)(load(i), i = 1, dof)
    
    ! Top Edge
    do i = (gNum+1), dof, (gNum+1)
        load(i) = 0
        disp(i) = 0.0
        gStiff(i, :) = 0.0
        gStiff(:, i) = 0.0
        gStiff(i, i) = 1.0
    end do
    
    ! Right Edge
    do i = ((gNum * (gNum + 1)) + 1),  dof
        load(i) = 0
        disp(i) = 0.0
        gStiff(i, :) = 0.0
        gStiff(:, i) = 0.0
        gStiff(i, i) = 1.0
    end do
    
    ! Debug Stiffness
    !write (*,4) ( (gStiff(n, i), i = 1, 10), n = 1, 10)
    
    !write(*,5) (load(k), k = 1, dof)
5   format (/"load(k) After "  / (E9.3))    
    
    !==================   Solve  ======================================
    !  Find all displacements
  
    !==================   Initialize ============
    matIn(:,:) = 0.0
    matIn = gStiff
    nSize = dof
    matL(:,:) = 0.0
    matLT(:,:) = 0.0
    disp(:) = 0.0
    y(:) = 0.0
    x(:) = 0.0
    
    !==================   decompose  ============
    
    do k = 1, nSize !outside loop
        !write(*,*)('k', (k))
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
    
    !=================   sub  ===================
    
    matLT = transpose(matL) !transpose this matrix
    
    do k = 1, nSize !forward sub
        !write(*,*)('k', (k))
        sum = 0.0 !clear the var
        do i = 1, nSize
            sum = sum + (matL(k,i) * y(i)) !sum factor
        end do
        y(k) = (load(k) - sum) / matL(k,k)
    end do
    
    do k = nSize, 1, -1 !back sub
        !write(*,*)('k', (k))
        sum = 0.0 !clear the var
        do i = nSize, 1, -1
            sum = sum + (matLT(k,i) * x(i)) !sum factor
        end do
        x(k) = (y(k) - sum) / matLT(k,k)
    end do

    
    disp = x !write to the var Gramoll uses.
    psiValues = disp
    
    
    !================================================  
    !psi debug
    !write(*,6) (disp(i), i = 1, dof)
6   format(/"Psi"/ (E9.3))    
    
    !find all torsion values
    
    do i = 1, totElements
        !write(*,*)('i',(i))
        torsionValues(i) = ((2 * area(i)) / 3 )* (psiValues(elemI(i)) + psiValues(elemJ(i)) + psiValues(elemM(i)))
    end do
    
    !write(*,*)('here')
    
    !find eq torque
    do i = 1, totElements
        torqueE = torqueE + torsionValues(i)
    end do
    
    
    
    !torque Debugging
    write(*,*)("---- Torque ---------------")
    write(*,201)(torqueE * 4)
201 format("      Tot torque " / "      ", (E10.4))    
    !write(*,*)(torqueE* 4)
    write(*,13)
    write(*,101)
    
    !find angle of twist
    angleTwist = (torque / (torqueE * 4 )) * 0.01
    
    !write(*,7)(psiTot)
7   format(/"Temp"/ (E9.3))

    temp = 0.0 
    
    !debugging
    !69  format(/"Total Torsion"/ (E9.3))
    !    !write(*,69)(torqueE)
    
    
    write(*,*) ("---- Global Info. ---------")
8   format("      Angle of Twist (rad/l)" / "      ",(E10.4))
    write(*,8)(angleTwist)
    write(*,13)
    write(*,*)(/" "/)
    
    write(*,*)("Stresses in pts A, B, C:")
    write(*,101)
    
9   format(/"area"/ (E9.3))
    !write(*,9)(area(1))
    
    ! take stress values at A

    temp = (2 * gNum) - 1
    k = temp
    shearX = angleTwist*(1 / (2 * area(1)))*((betI(k) * psiValues(elemI(k))) + (betJ(k) * psiValues(elemJ(k))) + (betM(k) * psiValues(elemM(k)))) / 0.01
    shearY = angleTwist*(1 / (2 * area(1)))*((gamI(k) * psiValues(elemI(k))) + (gamJ(k) * psiValues(elemJ(k))) + (gamM(k) * psiValues(elemM(k)))) / 0.01
    shearXY = sqrt(shearX**2 + shearY**2)
    shearArray(k) = shearXY
    
    
    write(*,*)("---- Point A --------------")
    !write(*,10)(shearX)
    !write(*,11)(shearY)
    write(*,12)(shearXY)
    write(*,13)
    write(*,101)
    
    
    ! Take shears for the whole quad
    do k = 1, totElements
        shearX = angleTwist*(1 / (2 * area(1)))*((betI(k) * psiValues(elemI(k))) + (betJ(k) * psiValues(elemJ(k))) + (betM(k) * psiValues(elemM(k)))) / 0.01
        shearY = angleTwist*(1 / (2 * area(1)))*((gamI(k) * psiValues(elemI(k))) + (gamJ(k) * psiValues(elemJ(k))) + (gamM(k) * psiValues(elemM(k)))) / 0.01
        shearXY = sqrt(shearX**2 + shearY**2)
        shearArray(k) = shearXY
    end do
    maxShearLoc = 0.0
    
    maxShear = MAXVAL(shearArray)
    maxShearLoc = MAXLOC(shearArray, DIM = 1)
    
    ! take stress values at B
    k = (totElements)
    shearX = angleTwist* (1 / (2 * area(1)))*((betI(k) * psiValues(elemI(k))) + (betJ(k) * psiValues(elemJ(k))) + (betM(k) * psiValues(elemM(k)))) / 0.01
    shearY = angleTwist*(1 / (2 * area(1)))*((gamI(k) * psiValues(elemI(k))) + (gamJ(k) * psiValues(elemJ(k))) + (gamM(k) * psiValues(elemM(k))))/ 0.01
    shearXY = sqrt(shearX**2 + shearY**2)
    temp = shearXY !temporarily store for avaerage
    !write(*,*)("k ", (k))
    !write(*,*)(temp)
    
    k = (totElements) - 1
    shearX = angleTwist*(1 / (2 * area(1)))*((betI(k) * psiValues(elemI(k))) + (betJ(k) * psiValues(elemJ(k))) + (betM(k) * psiValues(elemM(k))))/ 0.01
    shearY = angleTwist*(1 / (2 * area(1)))*((gamI(k) * psiValues(elemI(k))) + (gamJ(k) * psiValues(elemJ(k))) + (gamM(k) * psiValues(elemM(k))))/ 0.01
    shearXY = sqrt(shearX**2 + shearY**2)
    temp1 = shearXY !temporarily store for averaging
    !write(*,*)("k ", (k))
    !write(*,*)(temp1)
    !k = sqrt(temp**2 + temp1**2)
    !write(*,*) (k)
    shearXY = (temp + temp1) / 2 !average

    
    write(*,*)("---- Point B (Average) ----")
    !write(*,10)(shearX)
    !write(*,11)(shearY)
    write(*,12)(shearXY)
    write(*,13)
    write(*,101)
    
    
    ! take stress values at C
    temp = ((gNum * 2) * (gNum - 1)) + 2
    k = temp
    shearX = angleTwist*(1 / (2 * area(1)))*((betI(k) * psiValues(elemI(k))) + (betJ(k) * psiValues(elemJ(k))) + (betM(k) * psiValues(elemM(k))))/ 0.01
    shearY = angleTwist*(1 / (2 * area(1)))*((gamI(k) * psiValues(elemI(k))) + (gamJ(k) * psiValues(elemJ(k))) + (gamM(k) * psiValues(elemM(k))))/ 0.01
    shearXY = sqrt(shearX**2 + shearY**2)

    
    call CPU_TIME(t2)!call finish time
    elapsedTime = t2 - t1 !calculate delta t
    
    write(*,*)("---- Point C --------------")
    !write(*,10)(shearX)
    !write(*,11)(shearY)
    write(*,12)(shearXY)
    write(*,13)
    write(*,101)
    
10  format("      Shear X"/ (E9.3))
11  format("      Shear Y"/ (E9.3))    
12  format("      Shear XY"/     "      ", (E9.3))    
13  format(" ---------------------------")
14  format("      Max Shear", "                       " (E9.3))
15  format("      Max Shear Element Number", " " (F9.0))
    
    write(*,*)("---- Max Loc --------------")
    write(*,14)(maxShear)
    write(*,15)(maxShearLoc)
    write(*,13)
    
    !write out elapsed time
    write(*, "('elapsed time:', E14.8)") (elapsedTime)
    
    

    write (*,*) 'torsion End'
    
    end program torsion