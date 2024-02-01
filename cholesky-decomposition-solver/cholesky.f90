!  practice.f90 
!
!  FUNCTIONS:
!  practice - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: practice
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program practice

    implicit none
    
    ! Declare Variables, set array size later after matrix size is read in
    real*8, allocatable, dimension(:,:) :: matIn, matOut, matL !! input and output square 2D matrix
    !real*8, allocatable, dimension(:,:) :: b, j
    integer :: i, j, k, n, m, iRow, iCol, nSize, temp1, temp2   !! counters
   
    double precision :: aPrime, aPrimeRow, sum
    
    character (len=30) :: strBuffer                 ! string that holds command line argument, no .txt
    character (len=30) :: fileName, fileInName, fileOutName   ! actual file name, file in and out with .txt 
    
    
    ! Body of simple_matrix
    print *, 'This is the practice program'   !! free format writing, use print (write is used for formating)
    
    ! Get command line input parameters, call for each item
    call GETARG(1, strBuffer)           ! special function to read in argurments 
    read (strBuffer,*) fileName         ! reading in file name (no extension)
    
    fileInName = TRIM(fileName) // ".txt"          ! file in name with extension
    fileOutName = TRIM(fileName) // "_out.txt"     ! file out name with extension   
    
    
    !! Read in data file    ========================================================================
    open (UNIT=5, FILE=fileInName)     !! file number 5 standard for input
    read (5,*) nSize    
    
    !!  set 2D array size
    allocate (matIn(nSize, nSize))
    allocate (matOut(nSize, nSize))
    allocate (matL(nSize, nSize))
    !allocate (b(nSize, 1)) !this is on the RHS of the equation
    
    !! Finish reading in data now that matrix size is known
    read (5,*) ((matIn(i,j), j=1, nSize), i=1, nSize)

    !   Initializeout matrix =================================================
    matOut ( : , : ) = 0.0  ! or you can loop on matrix (slower)
    matL (: , : ) = 0.0
    !b (:,1) = 0.0
    !j (:,1) = 0.0
    
    aPrime = 0.0
    temp1 = 0.0
    temp2 = 0.0
    
    !!  write to input matrix to screen
    do i = 1, nSize
        write (*, '(*(F7.3))' ) ( matIn(i,j), j=1,nSize )
        !print *, i
    end do
    
    
    print *, "Now we'll print the L matrix"
!     ********************************* Cholesky Decomposition
    !do i = 1, nSize
    !    do j = 1, nSize
    !    matL(i,j) = matIn(i,j) ! put the matIn matrix in the matL matrix to prepare for L decomposition.
    !    end do
    !end do
        

    
    do k = 1, nSize !outside loop
        do i = 1, k !secondary loop
            if (k .eq. 1) then !fixes the matL(1,1) exception
                matL(k,k) = sqrt(matIn(k,k))
            else !does the rest of the loop
                sum = 0.0 !clear the sum variable
                do j = 1, i - 1
                    sum = sum + (matL(i,j) * matL(k,j)) !performs the summation of this factor
                end do
                matL(k,i) = (matIn(k,i) - sum) / matL(i,i) !subtracts the previous factor to the original matrix, and devides it by matL(i,i)
                
                sum = 0.0 !clear the sum
                do j = 1, k - 1
                    sum = sum + (matL(k,j)**2) !performs the summation of this factor 
                end do
                matL(k,k) = sqrt(matIn(k,k) - sum) !subtracts the previous factor to the origional matrix, takes the square root
            end if
        end do
    end do
    
    
    
    do i =1, nSize ! print out the cholesky lower triangular matrix, matL(*,*)
        write (*, '(*(F7.3))' ) ( matL(i,j), j=1, nSize )
    end do
    ! ******************************************* Cholesky Solve

    
    

    end program practice