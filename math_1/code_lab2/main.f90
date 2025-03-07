program matrix_norm
    use matrix_ops  ! Подключаем модуль с подпрограммами decomp и solve
    implicit none
    integer, parameter :: n = 5
    real(8) :: A(n,n), A_inv(n,n), B(n-1), R(n,n)
    real(8) :: a4_values(4) = [1.5d0, 1.01d0, 1.001d0, 1.0001d0]
    integer :: i, j, k, NPIV(n)
    real(8) :: D, normR, cond, work(n)

    do k = 1, 4
        ! Задаем вектор B
        B = [4.0d0, 3.0d0, 2.0d0, a4_values(k)]

        ! Формируем матрицу A
        do i = 1, n
            do j = 1, n
                if (j == 1) then
                    A(i,j) = 1.0d0
                elseif (i < j) then
                    A(i,j) = B(i)
                else
                    A(i,j) = 1.0d0
                end if
            end do
        end do

        ! Копируем A в A_inv перед разложением
        A_inv = A

        ! Выполняем LU-разложение
        call decomp(n, n, A_inv, cond, NPIV, work)

        ! Находим A^-1, решая AX = I
        do i = 1, n
            ! Единичная матрица в виде столбца
            R = 0.0d0
            R(i,i) = 1.0d0
            call solve(n, n, A_inv, R(:,i), NPIV)
        end do

        ! Вычисляем R = AA⁻¹ - I
        R = matmul(A, A_inv)
        do i = 1, n
            R(i,i) = R(i,i) - 1.0d0
        end do

        ! Вычисляем норму матрицы R (по максимальной сумме строк)
        normR = maxval(sum(abs(R), dim=2))

        print *, "Для a4 =", B(n-1), "норма R =", normR
    end do
end program matrix_norm