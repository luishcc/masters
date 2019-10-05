program main
implicit none

!Accuracy options:
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: qp = selected_real_kind(33, 4931)

!Declare var:
real(dp), allocatable::a(:,:)   
real(dp), allocatable :: x(:), b(:), stocktime(:)
integer::i, j, k, num, cont, numdavez
real(dp)::dt, t1, t2, sumdt, R
integer :: ndavez(100) = (/(i, i=100,10000, 100)/)

!Saving at:
open (10, file='linha-fortran.txt')

!Random number gen.:
R=(rand(0)*10.0)+0.0

!number of times done:
do k=1, 100
	!Update num:
	num = ndavez(k)

	!Vec alloc:
	allocate (a(1:num,1:num))
	allocate (x(1:num))
	allocate (b(1:num))

	!Average loop:
	sumdt=0
	do cont=1, 10
		t1=0
		t2=0
		dt=0
		do i=1, num
			b(i)=0.0
			x(i)=R
			do j=1, num
				a(i,j)=R
			enddo
		enddo

		!Main loop:
		call CPU_TIME( t1 )
		do i=1, num
			do j=1, num
				b(i) = b(i) + a(i,j)*x(j) 
			enddo
		enddo
		call CPU_TIME( t2 )

		dt = t2-t1
		sumdt = sumdt + dt
		write(*,*) cont	
	enddo
	!Ave value:
	sumdt = sumdt * 0.1

	!Writing results:
	write(*,*) sumdt, num
	write(10,206) sumdt, num
	deallocate (a, b, x)
	206 format (2x, E20.15,',',i6)
enddo

end
