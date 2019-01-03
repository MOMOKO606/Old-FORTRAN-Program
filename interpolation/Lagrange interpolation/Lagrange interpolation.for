*****************************************************************
*一维拉格朗日插值函数
*****************************************************************
	
	real function lagrange_interpolation(n,x,xx,f)
	real x(0:n),f(0:n)
	real xx,lxx
	lagrange_interpolation=0.0
	do k=0,n  
	  lxx=1.0
	  do i=0,n
	    if(i.ne.k)then
	    lxx=lxx*(xx-x(i))/(x(k)-x(i))
	    endif
	  enddo
	  lagrange_interpolation=lagrange_interpolation+f(k)*lxx
	enddo
	end function


*****************************************************************
*一维拉格朗日插值子程序
*****************************************************************

	subroutine lagrange_interpolation(n,x,xx,f,ff)
	real x(0:n),f(0:n)
	real xx,lxx,ff
	ff=0.0
	do k=0,n  
	  lxx=1.0
	  do i=0,n
	    if(i.ne.k)then
	    lxx=lxx*(xx-x(i))/(x(k)-x(i))
	    endif
	  enddo
	  ff=ff+f(k)*lxx
	enddo
	end subroutine


*****************************************************************
*二维拉格朗日插值子程序
*****************************************************************

	subroutine lagrange_2D(n,m,x,xx,f,ff)
	real x(0:n),f(0:n),xx(1:m),ff(1:m)
	real lxx
	do j=1,m
	  ff(j)=0.0
	  do k=0,n  
	    lxx=1.0
	    do i=0,n
	      if(i.ne.k)then
	      lxx=lxx*(xx(j)-x(i))/(x(k)-x(i))
	      endif
	    enddo
	    ff(j)=ff(j)+f(k)*lxx
	  enddo
	enddo
	end subroutine