      SUBROUTINE DECOMP(NDIM,N,A,COND,IPVT,WORK)
C
      INTEGER NDIM,N
      REAL A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)
C
C     ПРОГРАММА ВЫЧИСЛЯЕТ РАЗЛОЖЕНИЕ ВЕЩЕСТВЕННОЙ МАТРИЦЫ
C     ПОСРЕДСТВОМ ГАУССОВА ИСКЛЮЧЕНИЯ И ОЦЕНИВАЕТ
C     ОБУСЛОВЛЕННОСТЬ  МАТРИЦЫ.
C
C     ОНА ИСПОЛЬЗУЕТСЯ ДЛЯ ВЫЧИСЛЕНИЯ РЕШЕНИЯ
C     ЛИНЕЙНЫХ СИСТЕМ.
C
C     ВХОДНАЯ ИНФОРМАЦИЯ.
C
C     NDIM - ДАВЛЕННАЯ СТРОЧНАЯ РАЗМЕРНОСТЬ МАССИВА,
C           СОДЕРЖАЩЕГО A.
C
C     N    - ПОРЯДОК МАТРИЦЫ.
C
C     A    - МАТРИЦА, КОТОРУЮ НУЖНО РАЗЛОЖИТЬ.
C
C     ВЫХОДНАЯ ИНФОРМАЦИЯ.
C
C     A     СОДЕРЖИТ ВЕРХНЮЮ ТРЕУГОЛЬНУЮ МАТРИЦУ U
C           И ЭЛЕМЕНТАРНЫЕ ПЕРЕСТАНОВКИ ВЕРСИИ
C           НИЖНЕЙ ТРЕУГОЛЬНОЙ МАТРИЦЫ L, ТАКИЕ,
C           ЧТО (МАТРИЦА  ПЕРЕСТАНОВОК) *A=L*U
C
C     COND - ОЦЕНКА ОБУСЛОВЛЕННОСТИ A.
C           ДЛЯ ЛИНЕЙНОЙ СИСТЕМЫ A*X=B ИЗМЕНЕНИЯ В A И B
C           МОГУТ ВЫЗЫВАТЬ  ИЗМЕНЕНИЯ В X, ПРОПОРЦИОНАЛЬНЫЕ COND.
C           ЕСЛИ COND+1.0.EQ.COND, ТО  A В ПРЕДЕЛАХ МАШИННОЙ
C           ТОЧНОСТИ ЯВЛЯЕТСЯ ВЫРОЖДЕННОЙ МАТРИЦЕЙ. COND
C           ПОЛАГАЕТСЯ РАВНЫМ 1.0E+32, ЕСЛИ ОНАНРЯДНА ТОЧНАЯ
C           ВЫРОЖДЕННОСТЬ.
C
C     IPVT - ВЕКТОР ВЕДУЩИХ ЭЛЕМЕНТОВ.
C           IPVT(K) = ИНДЕКС K-ОЙ ВЕДУЩЕЙ СТРОКИ
C           IPVT(N) = (-1)**(ЧИСЛО ПЕРЕСТАНОВОК)
C
C     РАБОЧЕЕ ПОЛЕ. ВЕКТОР WORK ДОЛЖЕН БЫТЬ ОПИСАН И
C             ВКЛЮЧЕН В ВЫЗОВ. ЕГО ВХОДНОЕ СОДЕРЖАНИЕ ОДНАКО
C             НЕ ДАЕТ ВАЖНОЙ ИНФОРМАЦИИ.
C
C     ОПРЕДЕЛИТЕЛЬ МАТРИЦЫ A МОЖЕТ БЫТЬ НАЙДЕН НА ВЫХОДЕ
C     ПО ФОРМУЛЕ:
C          DET(A)=IPVT(N)*A(1,1)*A(2,2)*...*A(N,N).
C
      REAL EK,T,ANORM,YNORM,ZNORM
      INTEGER NM1,I,J,K,KP1,KB,KM1,M
C
      IPVT(N)=1
      IF(N.EQ.1)GO TO 80
      NM1=N-1
C
C     ВЫЧИСЛИТЬ 1-НОРМУ МАТРИЦЫ A
C
      ANORM=0.0
      DO 10 J=1,N
        T=0.0
        DO 5 I=1,N
          T=T+ABS(A(I,J))
    5   CONTINUE
        IF(T.GT.ANORM) ANORM=T
   10 CONTINUE
C
C     ГАУССОВО ИСКЛЮЧЕНИЕ С ЧАСТИЧНЫМ ВЫБОРОМ
C     ВЕДУЩЕГО ЭЛЕМЕНТА
C
      DO 35 K=1,NM1
        KP1=K+1
C
C       НАЙТИ ВЕДУЩИЙ ЭЛЕМЕНТ
C
        M=K
        DO 15 I=KP1,N
          IF(ABS(A(I,K)).GT.ABS(A(M,K))) M=I
   15   CONTINUE
        IPVT(K)=M
        IF(M.NE.K)IPVT(N)=-IPVT(N)
        T=A(M,K)
        A(M,K)=A(K,K)
        A(K,K)=T
C
C       ПРОПУСТИТЬ ЭТОТ ШАГ, ЕСЛИ ВЕДУЩИЙ ЭЛЕМЕНТ РАВЕН НУЛЮ
C
        IF(T.EQ.0.0)GO TO 35
C
C       ВЫЧИСЛИТЬ МНОЖИТЕЛИ
C
        DO 20 I=KP1,N
          A(I,K)=-A(I,K)/T
   20   CONTINUE
C
C       ПЕРЕСТАВИТЬ И ИСКЛЮЧИТЬ ПО СТОЛБЦАМ
C
        DO 30 J=KP1,N
          T=A(M,J)
          A(M,J)=A(K,J)
          A(K,J)=T
          IF(T.EQ.0.0)GO TO 30
          DO 25 I=KP1,N
            A(I,J)=A(I,J)+A(I,K)*T
   25     CONTINUE
   30   CONTINUE
   35 CONTINUE
C
C     COND=(1-НОРМА МАТРИЦЫ A)*(ОЦЕНКА ДЛЯ 1-НОРМЫ МАТРИЦЫ,
C     ОБРАТНОЙ К A)
C     ОЦЕНКА ПОЛУЧАЕТСЯ ПОСРЕДСТВОМ ОДНОГО ЦИКЛА МЕТОДА
C     ОБРАТНЫХ ИТЕРАЦИЙ ДЛЯ НАИМЕНЬШЕГО СИНГУЛЯРНОГО
C     ВЕКТОРА. ЭТО ТРЕБУЕТ РЕШЕНИЯ ДВУХ СИСТЕМ УРАВНЕНИЙ,
C     (ТРАНСПОНИРОВАННАЯ ДЛЯ A) *Y=E И A*Z=Y, ГДЕ E-ВЕКТОР
C     С +1 И -1, ВЫБРАННЫЙ ТАК, ЧТОБЫ МАКСИМИЗИРОВАТЬ
C     ВЕЛИЧИНУ Y.
C     ESTIMATE=(1-НОРМА Z)/(1-НОРМА Y)
C
C     РЕШИТЬ СИСТЕМУ (ТРАНСПОНИРОВАННАЯ ДЛЯ A)*Y=E
C
      DO 50 K=1,N
        T=0.0
        IF(K.EQ.1)GO TO 45
        KM1=K-1
        DO 40 I=1,KM1
          T=T+A(I,K)*WORK(I)
   40   CONTINUE
   45   EK=1.0
        IF(T.LT.0.0)EK=-1.0
        IF(A(K,K).EQ.0.0)GO TO 90
        WORK(K)=-(EK+T)/A(K,K)
   50 CONTINUE
      DO 60 KB=1,NM1
        K=N-KB
        T=WORK(K)
        KP1=K+1
        DO 55 I=KP1,N
          T=T+A(I,K)*WORK(I)
   55   CONTINUE
        WORK(K)=T
        M=IPVT(K)
        IF(M.EQ.K)GO TO 60
        T=WORK(M)
        WORK(M)=WORK(K)
        WORK(K)=T
   60 CONTINUE
C
      YNORM=0.0
      DO 65 I=1,N
        YNORM=YNORM+ABS(WORK(I))
   65 CONTINUE

      CALL SOLVE(NDIM,N,A,WORK,IPVT)
      
      ZNORM=0.0
      DO 70 I=1,N
        ZNORM=ZNORM+ABS(WORK(I))
   70 CONTINUE

      COND=ANORM*ZNORM/YNORM
      IF(COND.LT.1.0)COND=1.0
      RETURN

   80 COND=1.0
      IF(A(1,1).NE.0.0)RETURN

   90 CONTINUE
      COND=1.0E+32
      RETURN
      END
