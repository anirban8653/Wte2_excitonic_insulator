program linhard
    Use wte2_hamiltonian

    implicit none
    integer :: i, j, k, l, m, n, o, q, ii, jj
    integer, parameter :: kkm = 100
    integer, parameter :: ns = (kkm+1)**2
    integer, parameter :: qqm = 70
    integer, parameter :: nsq = (qqm+1)**2
    integer, parameter :: norb = 8
    real*8, dimension(:), allocatable :: kx, ky,qx,qy
    complex*16, parameter :: iota = (0.D0, 1.D0)
    real*8, parameter :: kT = 0.001
    real*8, dimension(:,:), allocatable :: nfk,nfq
    !complex*16, dimension(:,:,:), allocatable :: xi
    complex*16 :: xi 
    real*8, parameter :: omega = 0.06
    real*8, parameter :: eta = 0.005
    real, parameter :: a = 3.477
    real, parameter :: b = 6.249
    Real*8,Parameter::Pi=ACos(-1.D0) 

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
    allocate(qx(nsq))
    allocate(qy(nsq))
    allocate(nfk(ns,norb))    
    allocate(nfq(ns,norb)) 
    ! allocate(xi(ns,norb**2,norb**2))
    ! allocate(tensor_product(norb, norb))
    
    
    ! Open the data files
    Open(12, File="wte2_sus.dat", Status="Unknown") 
    Open(13, File="kmesh_file.dat", Status="Unknown")
    Open(14, File="qmesh_file.dat", Status="Unknown")

    print*, "reading data files ..."
   
    Do i = 1, ns
        Read(13, *) kx(i), ky(i)
    !    write(*, *) kx(i), ky(i)
    End Do
    Do i = 1, nsq
        Read(14, *) qx(i), qy(i)
    !    write(*, *) kx(i), ky(i)
    End Do
    print*, "data file reading done!" 
    ! Close the files
    Close(13)
    Close(14)

    do q = 1, nsq
        
        xi = (0.D0,0.D0)
        do k = 1, ns

            !---------------------------------------------------------!
            !               Hamiltonian (kx,ky)                       !
            !---------------------------------------------------------!
            
            Bk(k,:,:) = hamiltonian((kx(k))/a,(ky(k))/b)

            !-------------diagonalising the hamiltonian---------------!
            Call ZHEEV('V','U',Nm,Bk(k,:,:),LDA,Wk(k,:),WORK,LWORK,RWORK,INFO)

            If( INFO.Gt.0 ) Then
            Write(*,*)'The algorithm failed to compute eigenvalues.'
            Stop
            End If

            !------------calculating the fermi function---------------!

            nfk(k,:) = 1/(exp(Wk(k,:)/kT) + 1)

            !-------------------end part 1----------------------------!




            !---------------------------------------------------------!
            !               Hamiltonian (kxp,kyp)                     !
            !---------------------------------------------------------!
            
            Bq(k,:,:) = hamiltonian((kx(k)+qx(q))/a,(ky(k)+qy(q))/b)

            !-------------diagonalising the hamiltonian---------------!
            Call ZHEEV('V','U',Nm,Bq(k,:,:),LDA,Wq(k,:),WORK,LWORK,RWORK,INFO)

            If( INFO.Gt.0 ) Then
            Write(*,*)'The algorithm failed to compute eigenvalues.'
            Stop
            End If

            !------------calculating the fermi function---------------!
            nfq(k,:) = 1/(exp(Wq(k,:)/kT) + 1)
            
            !-------------------end part 1----------------------------!

            do i=1,norb
                do j=1,norb
                    ! do l = 1, norb
                    !     do m = 1, norb
                    !         Bks(k,(l-1)*norb+m) = Bk(k,i,l) * conjg(Bk(k,i,m))
                    !         Bqs(k,(l-1)*norb+m) = Bq(k,j,l) * conjg(Bq(k,j,m)) 
                    !     enddo 
                    ! enddo
                    ! do l = 1, norb**2
                    !     do m = 1, norb**2
                    !         ! xi(k,l,m) = xi(k,l,m)-(Bks(k,l)*Bqs(k,m)*(nfq(k,i) - &
                            ! nfk(k,j))/(omega + Wq(k,i) - Wk(k,j) + iota * eta))/ns
                    xi = xi-(nfq(k,i) - nfk(k,j))/ &
                    ((omega + Wq(k,i) - Wk(k,j) + iota * eta)*ns)   
                    !     enddo
                    ! enddo    
                enddo
            enddo
        enddo  

        write(12,*) qx(q)/a, qy(q)/b, real(xi), aimag(xi)
        print*, q    
    enddo
                 
    
end program linhard