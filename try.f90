module fft_mod
    integer,parameter::dp=selected_real_kind(15,300)
    real(kind=dp),parameter::pi=3.141592653589793238460_dp
contains
    recursive subroutine fft(x,state)!fft递归实现
        complex(kind=dp), dimension(:), intent(inout)::x
        complex(kind=dp)::t
        integer::N!信号点数
        integer::i!循环变量
        integer::state
        complex(kind=dp), dimension(:), allocatable::even,odd!偶序列/奇序列
        N=size(x)
        if(N .le. 1) return
        allocate(odd((N+1)/2))
        allocate(even(N/2))
        ! divide
        odd =x(1:N:2)
        even=x(2:N:2)
        ! conquer
        call fft(odd,state)
        call fft(even,state)
        ! combine
        ! 在这里实现并行
        if(state == 1)then
            x=conjg(x)
            end if
        !end if
        !$OMP PARALLEL
        !$OMP DO
        do i=1,N/2
            t=exp(cmplx(0.0_dp,-2.0_dp*pi*real(i-1,dp)/real(N,dp),kind=dp))*even(i)
            x(i)     = odd(i) + t
            x(i+N/2) = odd(i) - t
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        deallocate(odd)
        deallocate(even)
    end subroutine fft
end module fft_mod

program test
    use fft_mod
    complex(kind=dp),allocatable::data(:,:)
    complex(kind=dp)::temp
    integer::M,N,i,j
    temp=cmplx(-1.0_dp,0.0_dp)
    !读文件
    open(99,file="data.dat")
    read(99,*)M,N
    allocate(data(M,N))
    do j=1,N
        do i=1,M
            read(99,*)data(i,j)
        end do
    end do
    !!!!!!!!!!!!!!!!!!FFT
    !$OMP PARALLEL
    !$OMP DO
    do j=1,N
        call fft(data(:,j),0)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO
    do i=1,M
        call fft(data(i,:),0)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    !!!!!!!!!!!!!!!!!!!!Hilbert
    data=data*temp;
    !!!!!!!!!!!!!!!!!!!!!IFFT
    !$OMP PARALLEL
    !$OMP DO
    do j=1,N
        call fft(data(:,j),1)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO
    do i=1,M
        call fft(data(i,:),1)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    open(98,file="data1.dat")
    do j=1,N
        do i=1,M
            write(98,*)data(i,j)
        end do
    end do
end program test