********************************************************************
*  蒙特卡洛法一维积分(One Dimention Monte Carlo Integration)
********************************************************************
	subroutine MCI_1D(a,b,f,s)
	double precision  a,b,f,r,x,s,k
	real nrnd1
	r=1.0
	s=0.0
	do k=1,10000
	  x=a+(b-a)*nrnd1(r)
	  s=s+f(x)
	enddo
	s=(b-a)*s/10000.0
	end subroutine