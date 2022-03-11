c test eigen.f
	real*8 a(2,2)
	real*8 z(2,2),wr(2),wi(2)

	a(1,1) = 1.0d0
	a(1,2) = 3.0d0
	a(2,1) = 0.0d0
	a(2,2) = 2.0d 0
	N=2
	M=2
	call deigen(M,N,a,wr,wi,z)

	do 10 i=1,2
		write(6,*) i , wr(i),wi(i)
 10	continue

	stop
	end
