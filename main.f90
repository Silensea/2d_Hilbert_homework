module par_1dfft(x)
    implicit none
    integer,parameter::dp=selected_real_kind(15,300)
    real(kind=dp),parameter::pi=3.141592653589793238460_dp
contains
    recursive subroutine fft(x)!fft递归实现
        complex(kind=dp), dimension(:), intent(inout)::x
        complex(kind=dp)::temp!
        integer::N!信号点数
        integer::i!循环变量
        complex(kind=dp), dimension(:), allocatable::even,odd!偶序列/奇序列
        N=size(x)
        if(N .le. 1) return!
        allocate(odd((N+1)/2))
        allocate(even(N/2))
        ! divide
        odd =x(1:N:2)
        even=x(2:N:2)
        ! conquer
        call fft(odd)
        call fft(even)
        ! combine
        ! 在这里实现并行
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
end module par_1dfft

program test
    use fft_mod
    include"mpif.h"
    implicit none
    complex(kind=dp), dimension(8)::data
    integer::i

    integer::node,np,ierr
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,node,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

    ionode=(node==0)

    if(Ionode)then
        !读文件
        data=(/1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0/)
        n=size(data)
        do i=1,np
            call MPI_SEND(n,1,MPI_INTEGER,i,99,MPI_COMM_WORLD,ierr)
        end do
    else
        call MPI_RECV(n,1,MPI_INTEGER,0,99,MPI_COMM_WORLD,ierr)
    end if

    call fft(data)

    if(Ionode) then
        do i=1,np
            call MPI_RECV()
        end do
    end if

    do i=1,8
        write(*,'("(", F20.15, ",", F20.15, "i )")') data(i)
    end do

end program test

!module par_1dfft
!    implicit none
!    integer,parameter::dp=selected_real_kind(15,300)
!    real(kind=dp),parameter::pi=3.141592653589793238460_dp
!contains
!    recursive subroutine fft(x,temp,wn,N)!fft递归实现
!        complex(kind=dp), dimension(:), intent(inout)::x
!        complex(kind=dp)::temp!
!        integer::N!信号点数
!        integer::i!循环变量
!        complex(kind=dp), dimension(:), allocatable::temp,wn!偶序列/奇序列
!        N=size(x)
!        if(N .le. 1) return!
!        allocate(temp(N))
!        allocate(wn(N))
!
!        if(N==1)
!            call fft_base(x);
!            return
!        end if
!
!        !$OMP PARALLEL
!        !$OMP DO
!        !loop 1
!        do i=1,N/2
!            temp[i]=x[i]+x[i+N/2]
!            temp[i+n/2]=(x[i]-x[i+N/2])*wn[i]
!        end do
!        !$OMP END DO
!        !$OMP END PARALEL
!        !$OMP PARALLEL
!        !$OMP DO
!        !loop 2
!        do k=1,2
!            call fft(temp[N/2*k],x[N/2*k],wn,N/2)
!        end do
!        !$OMP END DO
!        !$OMP END PARALLEL
!
!        !$OMP PARALLEL
!        !$OMP DO
!        !loop 3
!        do i=0,n/2
!            x[2*i+1]=temp[i]
!            x[2*i]=temp[i+n/2]
!        end do
!        !$OMP END DO
!        !$OMP END PARALLEL
!
!    end subroutine fft
!    subroutine fft_base
!end module par_1dfft