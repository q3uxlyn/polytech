module matrix_ops
    implicit none
contains

    subroutine decomp(ndim,n,a,cond,ipvt,work)

        integer ndim,n
        real a(ndim,n),cond,work(n)
        integer ipvt(n)
        real ek,t,anorm,ynorm,znorm
        integer nm1,i,j,k,kp1,kb,km1,m

        ipvt(n)=1
        if(n.eq.1)go to 80
        nm1=n-1

        anorm=0.0
        do 10 j=1,n
            t=0.0
            do 5 i=1,n
                t=t+abs(a(i,j))
        5   continue
            if(t.gt.anorm) anorm=t
    10 continue

        do 35 k=1,nm1
        kp1=k+1

        m=k
        do 15 i=kp1,n
            if(abs(a(i,k)).gt.abs(a(m,k))) m=i
    15  continue
        ipvt(k)=m
        if(m.ne.k)ipvt(n)=-ipvt(n)
        t=a(m,k)
        a(m,k)=a(k,k)
        a(k,k)=t

        if(t.eq.0.0)go to 35

        do 20 i=kp1,n
            a(i,k)=-a(i,k)/t
    20  continue

        do 30 j=kp1,n
            t=a(m,j)
            a(m,j)=a(k,j)
            a(k,j)=t
            if(t.eq.0.0)go to 30
            do 25 i=kp1,n
            a(i,j)=a(i,j)+a(i,k)*t
    25     continue
    30   continue
    35 continue

        do 50 k=1,n
        t=0.0
        if(k.eq.1)go to 45
        km1=k-1
        do 40 i=1,km1
            t=t+a(i,k)*work(i)
    40   continue
    45   ek=1.0
        if(t.lt.0.0)ek=-1.0
        if(a(k,k).eq.0.0)go to 90
        work(k)=-(ek+t)/a(k,k)
    50 continue
        do 60 kb=1,nm1
        k=n-kb
        t=work(k)
        kp1=k+1
        do 55 i=kp1,n
            t=t+a(i,k)*work(i)
    55   continue
        work(k)=t
        m=ipvt(k)
        if(m.eq.k)go to 60
        t=work(m)
        work(m)=work(k)
        work(k)=t
    60 continue

        ynorm=0.0
        do 65 i=1,n
        ynorm=ynorm+abs(work(i))
    65 continue

        call solve(ndim,n,a,work,ipvt)
        
        znorm=0.0
        do 70 i=1,n
        znorm=znorm+abs(work(i))
    70 continue

        cond=anorm*znorm/ynorm
        if(cond.lt.1.0)cond=1.0
        return

    80 cond=1.0
        if(a(1,1).ne.0.0)return

    90 continue
        cond=1.0e+32
        return
    end

    subroutine solve(ndim,n,a,b,ipvt)

        integer ndim,n,ipvt(n)
        real a(ndim,n),b(n)
        integer kb,km1,nm1,kp1,i,k,m
        real t

        if(n.eq.1) go to 50
        nm1=n-1
        do 20 k=1,nm1
        kp1=k+1
        m=ipvt(k)
        t=b(m)
        b(m)=b(k)
        b(k)=t
        do 10 i=kp1,n
            b(i)=b(i)+a(i,k)*t
    10   continue
    20 continue

        do 40 kb=1,nm1
        km1=n-kb
        k=km1+1
        b(k)=b(k)/a(k,k)
        t=-b(k)
        do 30 i=1,km1
            b(i)=b(i)+a(i,k)*t
    30   continue
    40 continue
    50 b(1)=b(1)/a(1,1)
        return
    end        
        
end module matrix_ops