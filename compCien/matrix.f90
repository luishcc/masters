program main
implicit none

real, allocatable::a(:,:)   
real, allocatable :: x(:), b(:)
integer::i, j, num, cont
real::dt, t1, t2, sumdt

open (10, file='result.txt')

num = 10  
 
allocate (a(1:num,1:num))
allocate (x(1:num))
allocate (b(1:num))


sumdt=0
do cont=1, 10

t1=0
t2=0
dt = 0

do i=1, num
  b(i)=0.0
  x(i)=1.0
  do j=1, num
    a(i,j)=1.0
  enddo
enddo

call CPU_TIME( t1 )

do i=1, num
  do j=1, num
    b(i) = b(i) + a(i,j)*x(i) 
    enddo
  enddo

call CPU_TIME( t2 )

dt = t2-t1

sumdt = sumdt + dt

enddo

sumdt = sumdt * 0.001

write(*,*) sumdt

write(10,206) sumdt, num

deallocate (a, b, x)

206 format (2x, f12.6,';',i6)

end
