program main
    use matrix_ops
    implicit none
    integer, parameter :: n = 5, ndim = n
    real :: A(ndim, n), A_inv(ndim, n), B(n-1), cond, work(n)
    integer :: ipvt(n), i, j
    real :: R(ndim, n), I(ndim, n), normR
    real, dimension(4) :: var_values = [1.5, 1.01, 1.001, 1.0001]
    integer :: k
    
    ! Задание вектора B
    B = [4.0, 3.0, 2.0, 0.0] ! Последний элемент заменяется на var позже
    
    do k = 1, 4
        B(4) = var_values(k) ! Устанавливаем текущее значение var
        print *, "------------------------------"
        print *, "For var =", B(4)
        
        ! Формирование матрицы A
        do i = 1, n
            do j = 1, n
                if (j < i) then
                    A(i, j) = 1.0
                else if (j > i) then
                    A(i, j) = B(j-1)
                else
                    A(i, j) = 1.0
                end if
            end do
        end do
        
        print *, "Matrix A:"
        do i = 1, n
            print *, A(i, :)
        end do
        
        ! Копируем A в A_inv для получения A^(-1)
        A_inv = A
        
        ! Вычисление A^(-1)
        call decomp(ndim, n, A_inv, cond, ipvt, work)
        do i = 1, n
            work = 0.0
            work(i) = 1.0
            call solve(ndim, n, A_inv, work, ipvt)
            A_inv(:, i) = work
        end do
        
        print *, "Inverse Matrix A_inv:"
        do i = 1, n
            print *, A_inv(i, :)
        end do
        
        ! Вычисление R = AA^(-1) - E
        R = matmul(A, A_inv)
        I = 0.0
        do i = 1, n
            I(i, i) = 1.0
        end do
        R = R - I
        
        print *, "Matrix R = AA^(-1) - I:"
        do i = 1, n
            print *, R(i, :)
        end do
        
        ! Вычисление нормы матрицы R
        normR = 0.0
        do i = 1, n
            do j = 1, n
                normR = normR + abs(R(i, j))
            end do
        end do
        
        ! Вывод результатов
        print *, "Norm of R:", normR
    end do
    
end program main
