program duschinsky_read
    implicit none
    external sgetrf
    integer :: natom = 21, nmodes, arg_count, i, j, k, n, row, nsteps, test_IPIV(2), test_INFO
    character(len = 80) :: arg_value, duschinsky, displacement, singlet, triplet, filename
    character(len = 2000) :: line
    integer :: maxlines = 10000, counter
    real(8) :: time, atomic_time, adiabatic, damping, starting_point, end_point, temperature, c_cm, partition = 1.0, test_partition = 1.0, boltzmann = 3.166811563455607e-06, integrated_result, time_step, soc, total_time
    real(8), dimension(1) :: test_s1_freq, test_t1_freq, test_displacement
    real(8), dimension(1,1) :: test_duschinsky = 1.0
    real(8), dimension(:,:), allocatable :: duschinsky_matrix
    real(8), dimension(:), allocatable :: displacement_vector,s1_freq,t1_freq, tsteps
    double complex :: correlation, ehamiltonian, integrand, test_determinant, sqrt_previous = CMPLX(0, 0), new_sqrt
    double complex, dimension(:,:), allocatable :: t1_sin,t1_tan,s1_sin,s1_tan
    double complex, dimension(2,2) :: test_matrix
    logical :: duschinsky_exists, displacement_exists, singlet_exists, triplet_exists

    nmodes = 3 * natom - 6
    atomic_time = 2.418884326505e-17
    c_cm = 2.998e10
! Insert the adiabatic energy difference
    adiabatic = 0.001121619999935
! Insert the number of time steps you will take
    nsteps = 1000000
! Insert the temperature
    temperature = 300.0
! Insert the time step, in atomic units
    time_step = 0.1
    integrated_result = 0.5 * time_step
! Insert the SOC matrix element
    soc = SQRT(3.0) * 0.153 * c_cm * atomic_time * 2 * 3.141592
    end_point = time_step * nsteps
    starting_point = -end_point
    write(*,*) 'The integration will be performed from', starting_point * atomic_time * 1e15, 'femtoseconds to', end_point * atomic_time * 1e15, 'femtoseconds with time interval of', time_step * atomic_time * 1e15, 'femtoseconds'
    damping = 0.0 * c_cm * atomic_time * 2 * 3.141592
    allocate(duschinsky_matrix(nmodes,nmodes))
    allocate(displacement_vector(nmodes))
    allocate(s1_sin(nmodes,nmodes))
    allocate(s1_tan(nmodes,nmodes))
    allocate(t1_sin(nmodes,nmodes))
    allocate(t1_tan(nmodes,nmodes))
    allocate(tsteps(nsteps))
    allocate(s1_freq(nmodes))
    allocate(t1_freq(nmodes))
    !Get number of command line arguments
    arg_count= COMMAND_ARGUMENT_COUNT()

    if (arg_count<4) then
        write(*,*) 'ERROR: You did not give enough input'
        stop
    end if

    call GET_COMMAND_ARGUMENT(1, arg_value)
    duschinsky=arg_value
    call GET_COMMAND_ARGUMENT(2, arg_value)
    displacement=arg_value
    call GET_COMMAND_ARGUMENT(3, arg_value)
    singlet=arg_value
    call GET_COMMAND_ARGUMENT(4, arg_value)
    triplet=arg_value

    !Check if the files exist
    inquire(file=duschinsky, exist=duschinsky_exists)
    inquire(file=displacement, exist=displacement_exists)
    inquire(file=singlet, exist=singlet_exists)
    inquire(file=triplet, exist=triplet_exists)
    if(.not. duschinsky_exists) then
        write (*,*) "ERROR: Duschinsky matrix file does not exist"
        stop
    else if(.not. displacement_exists) then
        write (*,*) "ERROR: Displacement vector file does not exist"
        stop
    else if(.not. singlet_exists) then
        write (*,*) "ERROR: Singlet frequency file does not exist"
        stop
    else if(.not. triplet_exists) then
        write (*,*) "ERROR: Triplet frequency file does not exist"
        stop
    else
        write (*,*) "All files are present! Continuing with the calculation"
    end if

    !Open and read the files
    open(unit=10, file=duschinsky, status='old', action='read')
    counter = 1
    do i= 1, maxlines
        read(10, '(A)', iostat=j) line
        if (j/=0) exit
        row=counter
        read(line,'(57f25.20)') duschinsky_matrix(row,:)
        counter = counter+1
    end do
    close(unit=10)


    open(unit=10, file=displacement, status='old', action='read')
    counter=1
    do i=1,maxlines
        read(10,'(A)', iostat=j) line
        if (j/=0) exit
        row=counter
        read(line,'(57f25.20)') displacement_vector(row)
        counter=counter+1
    end do
    close(unit=10)


    open(unit=10, file=singlet, status='old', action='read')

    counter=1
    do i=1,maxlines
        read(10,'(A)', iostat=j) line
        if (j/=0) exit
        row=counter
        read(line,'(57f25.20)') s1_freq(row)
!        s1_freq(row)=1e-15/atomic_time*s1_freq(row)
        counter=counter+1
    end do
    close(unit=10)


    open(unit=10, file=triplet, status='old', action='read')
    counter=1
    do i=1,maxlines
        read(10,'(A)', iostat=j) line
        if (j/=0) exit
        row=counter
        read(line,'(57f25.20)') t1_freq(row)
!        t1_freq(row)=1e-15/atomic_time*t1_freq(row)
        counter=counter+1
    end do
    close(unit=10)

    test_s1_freq(1) = s1_freq(nmodes)
    test_t1_freq(1) = t1_freq(nmodes)
    test_displacement(1) = displacement_vector(nmodes)
    do i=1,nmodes
        partition = partition * EXP(-0.5 * s1_freq(i) / (boltzmann * temperature)) / (1 - EXP(-s1_freq(i) / (boltzmann * temperature)))
    end do
    write(*,*) "The partition is", partition

    test_partition = EXP(-0.5 * test_s1_freq(1) / (boltzmann * temperature)) / (1 - EXP(-test_s1_freq(1) / (boltzmann * temperature)))
    filename = 'correlation_long.txt'
    open(unit = 10, file = filename, status = 'replace', action = 'write')
    do i = 1, nsteps
        time = i * time_step
        call RHO_DSO(nmodes, t1_freq, s1_freq, duschinsky_matrix, displacement_vector, time, temperature, partition, adiabatic, damping, sqrt_previous, new_sqrt, correlation)
        integrated_result = integrated_result + time_step * REAL(correlation)
        sqrt_previous = new_sqrt
        write(10, *) REAL(correlation)
    end do
    close(10)
    write(*,*) "The integrated result is", 2 * integrated_result
    write(*,*) "The rate is", soc * soc * 2 * integrated_result / atomic_time, "s-1"
!    test_matrix(1,1) = CMPLX(1.55,0.5)
!    test_matrix(1,2) = CMPLX(0.3,-0.11)
!    test_matrix(2,1) = CMPLX(1.5,-0.7)
!    test_matrix(2,2) = CMPLX(0.9,-0.4)
!    call ZGETRF(2, 2, test_matrix, 2, test_IPIV, test_INFO)
!    test_determinant = CMPLX(1.0,0.0)
!    do i = 1,2
!        test_determinant = test_determinant * test_matrix(i,i)
!    end do
!    write(*,*) "The determinant of the test matrix", test_determinant

!############
!Subroutines
!############

    contains
        subroutine VECTORT1(nmodes,T1FREQ,time,sin_matrix,tan_matrix)
            real(8), intent(in) :: time
            integer :: i
            integer, intent(in) :: nmodes
            real(8), dimension(nmodes), intent(in) :: T1FREQ
            double complex, dimension(nmodes,nmodes),intent(out) :: sin_matrix, tan_matrix
            do i=1,nmodes
                sin_matrix(i,i)=T1FREQ(i)/SIN(T1FREQ(i)*time)
                tan_matrix(i,i)=T1FREQ(i)/TAN(T1FREQ(i)*time)
            end do
        end subroutine VECTORT1

        subroutine VECTORS1(nmodes,S1FREQ,time,temperature,sin_matrix,tan_matrix)
            integer :: i
            integer, intent(in) :: nmodes
            real(8), intent(in) :: time,temperature
            real(8) :: boltzmann = 3.166811563455607e-06, atomic_time = 2.418884326505e-17, unit_time = 2.418884326505e-17
            real(8), dimension(nmodes), intent(in) :: S1FREQ
            double complex, dimension(nmodes,nmodes),intent(out) :: sin_matrix, tan_matrix
            do i=1,nmodes
                sin_matrix(i,i)=S1FREQ(i)/SIN(S1FREQ(i)*CMPLX(-time,-1/(boltzmann*(unit_time/atomic_time)*temperature)))
                tan_matrix(i,i)=S1FREQ(i)/TAN(S1FREQ(i)*CMPLX(-time,-1/(boltzmann*(unit_time/atomic_time)*temperature)))
            end do
        end subroutine VECTORS1

        subroutine MATRIX_CREATION(nmodes,T1FREQ,S1FREQ,duschinsky_matrix,displacement_vector,time,temperature,t_sin_matrix,t_tan_matrix,s_sin_matrix,s_tan_matrix,E_matrix,A_matrix,B_matrix,K_matrix,f_vector)
            integer :: i,j,k
            integer, intent(in) :: nmodes
            real(8) :: time, temperature
            real(8), dimension(nmodes), intent(in) :: displacement_vector, T1FREQ, S1FREQ
            real(8), dimension(nmodes, nmodes), intent(in) :: duschinsky_matrix
            double complex, dimension(nmodes, nmodes), intent(out) :: t_sin_matrix,t_tan_matrix,s_sin_matrix,s_tan_matrix
            double complex, dimension(2 * nmodes), intent(out) :: f_vector
            double complex, dimension(nmodes, nmodes), intent(out) :: E_matrix,B_matrix,A_matrix
            double complex, dimension(2 * nmodes, 2 * nmodes), intent(out) :: K_matrix
            do i = 1, nmodes
                f_vector(i) = 0.0
                do j = 1, nmodes
                    A_matrix(i,j) = 0.0
                    B_matrix(i,j) = 0.0
                    K_matrix(2 * i, 2 * j) = 0.0
                    K_matrix(2 * i + 1, 2 * j) = 0.0
                    K_matrix(2 * i, 2 * j + 1) = 0.0
                    K_matrix(2 * i + 1, 2 * j + 1) = 0.0
                end do
            end do
            call VECTORT1(nmodes,T1FREQ,time,t_sin_matrix,t_tan_matrix)
            call VECTORS1(nmodes,S1FREQ,time,temperature,s_sin_matrix,s_tan_matrix)
            do i = 1, nmodes
                E_matrix(i,i) = s_tan_matrix(i,i) - s_sin_matrix(i,i)
                A_matrix(i,i) = t_sin_matrix(i,i)
                B_matrix(i,i) = t_tan_matrix(i,i)
                do j=1,nmodes
                    do k=1,nmodes
                        A_matrix(i,j) = A_matrix(i,j) + duschinsky_matrix(k,i) * s_sin_matrix(k,k) * duschinsky_matrix(k,j)
                        B_matrix(i,j) = B_matrix(i,j) + duschinsky_matrix(k,i) * s_tan_matrix(k,k) * duschinsky_matrix(k,j)
                    end do
                    K_matrix(i, j) = B_matrix(i, j)
                    K_matrix(nmodes + i, nmodes + j) = B_matrix(i, j)
                    K_matrix(nmodes + i, j) = -A_matrix(i, j)
                    K_matrix(i, nmodes + j) = -A_matrix(i, j)
                end do
            end do
            do i = 1, nmodes
                do j = 1, nmodes
                    f_vector(i) = f_vector(i) + displacement_vector(j) * E_matrix(j, j) * duschinsky_matrix(j, i)
                    f_vector(nmodes + i) = f_vector(nmodes + i) + displacement_vector(j) * E_matrix(j, j) * duschinsky_matrix(j, i)
                end do
            end do
        end subroutine MATRIX_CREATION


        subroutine KMATRIX_PROCESSING(N, A, determinant, A_inverse)
            external zgetrf, zgetri
            integer, intent(in) :: N
            double complex, intent(in) :: A(N,N)
            double complex, intent(out) :: determinant, A_inverse(N,N)
            integer :: IPIV(N), INFO, i
            double complex, allocatable :: WORK(:)
            allocate(WORK(2 * N))
            A_inverse=A
            ! Perform LU decomposition
            call ZGETRF(N, N, A_inverse, N, IPIV, INFO)
            if (INFO == 0) then
 !               write(*,*) "LU decomposition successful."
                determinant = CMPLX(1.0,0.0)
                do i = 1, N
                    determinant = determinant * A_inverse(i,i)
                    if (IPIV(i)/=i) then
                        determinant = determinant * (-1)
                    end if
                end do
                call zgetri(N, A_inverse, N, IPIV, WORK, 2 * N, INFO)
            else
                write(*,*) "Error in LU decomposition."
                return
            end if
            deallocate(WORK)
        end subroutine KMATRIX_PROCESSING

        subroutine diagonal_product(nmodes,A,answer)
            integer, intent(in) :: nmodes
            double complex, intent(inout) :: A(nmodes,nmodes)
            double complex, intent(out) :: answer
            integer :: i
            answer = 1.0
            do i = 1,nmodes
                answer = answer * A(i,i)
            end do
        end subroutine diagonal_product


        subroutine RHO_DSO(nmodes,T1FREQ,S1FREQ,duschinsky_matrix,displacement_vector,time,temperature,partition,adiabatic,damping,sqrt_previous,sqrt_updated,correlation)
            integer, intent(in) :: nmodes
            integer :: i,j
            real(8) :: sqrt_difference_abs, sqrt_sum_abs, damping
            real(8), intent(in) :: time, temperature, partition, adiabatic
            real(8), dimension(nmodes), intent(in) :: displacement_vector, S1FREQ, T1FREQ
            real(8), dimension(nmodes,nmodes), intent(in) :: duschinsky_matrix
            double complex :: k_determinant, prod_sin_sin, prod_tri_sin, fKf, dEd, new_sqrt, sqrt_difference, sqrt_sum
            double complex, dimension(nmodes) :: temp_d_vector
            double complex, dimension(nmodes,nmodes) :: s_sin_matrix,s_tan_matrix,t_sin_matrix,t_tan_matrix,E_matrix,A_matrix,B_matrix
            double complex, dimension(2 * nmodes) :: f_vector, temp_f_vector
            double complex, dimension(2 * nmodes, 2 * nmodes) :: K_matrix,K_inverse
            double complex, intent(in) :: sqrt_previous
            double complex, intent(out) :: correlation, sqrt_updated
            call MATRIX_CREATION(nmodes,T1FREQ,S1FREQ,duschinsky_matrix,displacement_vector,time,temperature,t_sin_matrix,t_tan_matrix,s_sin_matrix,s_tan_matrix,E_matrix,A_matrix,B_matrix,K_matrix,f_vector)
            call KMATRIX_PROCESSING(2 * nmodes,K_matrix,k_determinant,K_inverse)
            call diagonal_product(nmodes,t_sin_matrix,prod_tri_sin)
            call diagonal_product(nmodes,s_sin_matrix,prod_sin_sin)
            temp_f_vector = MATMUL(f_vector, K_inverse)
            temp_d_vector = MATMUL(displacement_vector, E_matrix)
            fKf = CMPLX(0, 0)
            dEd = CMPLX(0, 0)
            do i = 1, 2 * nmodes
                if (i < nmodes + 1) then
                    dEd = dEd + temp_d_vector(i) * displacement_vector(i)
                end if
                fKf = fKf + temp_f_vector(i) * f_vector(i)
            end do
            new_sqrt = SQRT(prod_tri_sin / partition / k_determinant * prod_sin_sin / partition)
            sqrt_difference = sqrt_previous - new_sqrt
            sqrt_sum = sqrt_previous + new_sqrt
            sqrt_difference_abs = REAL(sqrt_difference)**2 + AIMAG(sqrt_difference)**2
            sqrt_sum_abs = REAL(sqrt_sum)**2 + AIMAG(sqrt_sum)**2

            if (sqrt_previous .EQ. CMPLX(0, 0)) then
                sqrt_updated = new_sqrt
            ELSE
                if (sqrt_difference_abs > sqrt_sum_abs) then
                    sqrt_updated = - new_sqrt
                else
                    sqrt_updated = new_sqrt
                end if
            end if
            correlation = EXP(CMPLX(0, -1) * adiabatic * time) * sqrt_updated * EXP(CMPLX(0, -1) * (0.5 * fKf - dEd)) * EXP(- damping * time)
        end subroutine RHO_DSO

end program duschinsky_read
