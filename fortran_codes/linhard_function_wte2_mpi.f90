program linhard
    Use wte2_hamiltonian
    

    implicit none
    include 'mpif.h'

    integer :: i, j, k, l, m, n, o, q, ii, jj
    integer :: it, base, core, job ,mpi_err, mpi_size, num_procs, mpi_rank
    integer, parameter :: kkm = 100
    integer, parameter :: ns = (kkm+1)**2
    integer, parameter :: qqm = 69
    integer, parameter :: nsq = (qqm+1)**2
    integer, parameter :: norb = 8
    real*8, dimension(:), allocatable :: kx, ky,qx,qy
    complex*16, parameter :: iota = (0.D0, 1.D0)
    real*8, parameter :: kT = 0.001D0
    real*8, dimension(:,:), allocatable :: nfk,nfq
    !complex*16, dimension(:,:,:), allocatable :: xi
    complex*16 :: xi 
    real*8, parameter :: omega = 0.0
    real*8, parameter :: eta = 0.002
    real, parameter :: a = 3.477
    real, parameter :: b = 6.249
    integer, allocatable :: pardis(:,:)
    complex*16, dimension(:), allocatable :: res, res_red
    
   
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
    allocate(res(nsq))
    allocate(res_red(nsq))
    ! allocate(xi(ns,norb**2,norb**2))
    ! allocate(tensor_product(norb, norb))



    
    ! Open the data files
    Open(12, File="wte2_sus_mpi.dat", Status="Unknown") 
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

    !---------------------------starting mpi-----------------------------!
   CALL MPI_Init(mpi_err)
   CALL MPI_Comm_size(MPI_COMM_WORLD,mpi_size,mpi_err)
   CALL MPI_Comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

   mpi_rank=mpi_rank+1

   ALLOCATE(pardis(mpi_size,nsq))
   base=nsq/mpi_size

   core=1
   job=1
   do i=1,mpi_size
      do k=1,nsq
         pardis(i,k)=0
      enddo
   enddo



   do i=1,nsq
      pardis(core,job)=i
      job=job+1

      if (job.gt.base) then
      core=core+1
      job=1
      endif
   enddo

    it=1
    do while ((pardis(mpi_rank,it).gt.0).and.(it.le.nsq))   
        q = pardis(mpi_rank,it)
    
        
        xi = (0.D0,0.D0)
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
                    xi = xi-(nfq(k,i) - nfk(k,j))/ &
                    ((omega + Wq(k,i) - Wk(k,j) + iota * eta)*ns)   
                    !     enddo
                    ! enddo    
                enddo
            enddo
        enddo  
        
        res(q) = xi
        it = it+1
        ! write(12,*) kx(q)/a, ky(q)/b, real(xi), aimag(xi)
        ! print*, q    
    enddo

    !-------------------------------- Ending MPI---------------------------------!
    call mpi_allreduce(res,res_red, 2*nsq,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, mpi_err)
    CALL MPI_Barrier(MPI_COMM_WORLD,mpi_err)

    !write(*,*) "forking data from different cores to a sing file ..."
   

    if (mpi_rank == 1) then
        do i = 1,nsq
        write(12,*) qx(i)/a, qy(i)/b, real(res_red(i)), aimag(res_red(i))
        enddo
        write(*,*) "4. data generated successfully!"
    end if

    deallocate(pardis)
    call MPI_Finalize(mpi_err)
    close(12)

               
    
end program linhard