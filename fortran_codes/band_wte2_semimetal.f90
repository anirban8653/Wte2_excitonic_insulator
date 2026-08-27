program band
    Use wte2_hamiltonian

    implicit none
    integer :: i, j, k, l, m, n, o, q, ii, jj
    integer, parameter :: kkm = 300
    integer, parameter :: ns = (kkm+1)**2
    integer, parameter :: num = 166464
    integer, parameter :: norb = 8
    real*8, dimension(:), allocatable :: kx, ky
    complex*16, parameter :: iota = (0.D0, 1.D0)
    real*8, parameter :: kT = 0.001
    real*8, dimension(:,:), allocatable :: nfk,nfq
    !complex*16, dimension(:,:,:), allocatable :: xi
    complex*16 :: xi 
    real*8, parameter :: omega = 0.D0
    real*8, parameter :: eta = 0.005
    real, parameter :: a = 3.477
    real, parameter :: b = 6.249
    ! complex*16, dimension(:, :), allocatable :: tensor_product
   
    Integer,Parameter :: nmax=norb
    Integer,Parameter :: Nm=nmax
    Integer,Parameter :: LDA=Nm
    Integer :: INFO
    Integer,Parameter ::LWORK=2*Nm-1
    Real*8 :: RWORK(3*Nm-2)
    Complex*16 :: Bk(ns, LDA,Nm), Bq(ns, LDA,Nm), Bks(ns, Nm**2), Bqs(ns, Nm**2)
    Complex*16 :: WORK(LWORK)
    real*8 :: Wk(ns, Nm), Wq(ns, Nm)
    External :: ZHEEV
    
    

    allocate(kx(ns))
    allocate(ky(ns))
    allocate(nfk(ns,norb))    
    allocate(nfq(ns,norb)) 
    ! allocate(xi(ns,norb**2,norb**2))
    ! allocate(tensor_product(norb, norb))
    
    
    ! Open the data files
    Open(12, File="wte2_band.dat", Status="Unknown") 
    Open(13, File="kmesh_file.dat", Status="Unknown")

    print*, "reading data files ..."
   
    Do i = 1, ns
        Read(13, *) kx(i), ky(i)
    !    write(*, *) kx(i), ky(i)
    End Do
    print*, "data file reading done!" 
    ! Close the files
    Close(13)

    do k = 1, ns

        !---------------------------------------------------------!
        !               Hamiltonian (kx,ky)                       !
        !---------------------------------------------------------!
        
        Bk(k,:,:) = hamiltonian(kx(k)/a,ky(k)/b)

        !-------------diagonalising the hamiltonian---------------!
        Call ZHEEV('V','U',Nm,Bk(k,:,:),LDA,Wk(k,:),WORK,LWORK,RWORK,INFO)

        If( INFO.Gt.0 ) Then
        Write(*,*)'The algorithm failed to compute eigenvalues.'
        Stop
        End If
        
        write(12,*) kx(k)/a,ky(k)/b,Wk(k,:)
    
    enddo
end program band

