********************************************************
*  高振荡函数积分子程序测试程序。
********************************************************
	external ff
	dimension aa(4),bb(4)
	double precision a,b,aa,bb,s1,s2
	data aa/0.0,1.0,0.0,-3.0/
	data bb/6.2831852,1.0,-6.2831852,-3.0/
	a=0.0
	b=6.2831852
	m=30
	n=4
	call I_HOF(a,b,n,m,aa,bb,s1,s2)
	write(*,*) s1,s2
	end
	
	subroutine I_HOF(a,b,n,m,aa,bb,s1,s2)
	dimension aa(n),bb(n)
	double precision a,b,aa,bb,s1,s2
	s1=0.0
	s2=0.0
	mm=1.0
	do k=0,n-1
	  mm=mm*m
	  s1=s1+(bb(k+1)*ff(0,k,m*b)-aa(k+1)*ff(0,k,m*a))/mm
	  s2=s2+(aa(k+1)*ff(1,k,m*a)-bb(k+1)*ff(1,k,m*b))/mm
	enddo
	end subroutine

	double precision function ff(l,k,x)
	double precision x
	if(l.eq.0)then
	  select case(mod(k,4))
	  case(0)
	    ff=sin(x)
	  case(1)
		ff=cos(x)
	  case(2)
		ff=-sin(x)
	  case(3)
		ff=-cos(x)
	  end select
	else
	  select case(mod(k,4))
	  case(0)
	    ff=cos(x)
	  case(1)
		ff=-sin(x)
	  case(2)
		ff=-cos(x)
	  case(3)
		ff=sin(x)
	  end select
	endif
	end function


