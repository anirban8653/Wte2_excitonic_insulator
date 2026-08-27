PROGRAM GenerateDataFile
   INTEGER, PARAMETER :: num_divisions = 300
   INTEGER, PARAMETER :: num_divisions2 = 99
   REAL :: pi, delta, x
   INTEGER :: i
   CHARACTER(20) :: file_name, file_name2
  

   ! Calculate the value of pi
   pi = ACOS(-1.0)

   ! Calculate the step size (delta) for each division
   delta = (2.0 * pi) / REAL(num_divisions)

   ! Open the file for writing
   WRITE(file_name, '(A)') 'kmesh_file.dat'
   OPEN(UNIT=10, FILE=file_name, STATUS='REPLACE')

   ! Generate and write the data to the file
   DO i = 1, num_divisions + 1
   DO j = 1, num_divisions + 1
  
      x = -pi + (REAL(i - 1) * delta)
      y = -pi + (REAL(j - 1) * delta)
     
      WRITE(10, *) x, y
!      WRITE(*, *) x, y
   ENDDO
   ENDDO
  

   ! Close the file
   CLOSE(10)

   WRITE(*, *) 'K Data file generated successfully!'


   ! Calculate the step size (delta) for each division
   delta2 = (2.0 * pi) / REAL(num_divisions2)

   ! Open the file for writing
   WRITE(file_name2, '(A)') 'qmesh_file.dat'
   OPEN(UNIT=12, FILE=file_name2, STATUS='REPLACE')

   ! Generate and write the data to the file
   DO i = 1, num_divisions2 + 1
   DO j = 1, num_divisions2 + 1
  
      x = -pi + (REAL(i - 1) * delta2)
      y = -pi + (REAL(j - 1) * delta2)
     
      WRITE(12, *) x, y
!      WRITE(*, *) x, y
   ENDDO
   ENDDO
  

   ! Close the file
   CLOSE(12)

   WRITE(*, *) 'q Data file generated successfully!'

END PROGRAM GenerateDataFile


