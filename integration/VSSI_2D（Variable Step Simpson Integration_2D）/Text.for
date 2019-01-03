***************************************************************
*  变步长辛普森二重积分子程序测试程序
****************************************************************
	external f,fs
	double precision a,b,f,s
	a=0.0
	b=1.0
	eps=0.0001
	call VSSI_2D(a,b,eps,fs,f,s)
	write(*,30) s
30    format(1x,'s=',d13.6)
	end
	
	subroutine VSSI_2D(a,b,eps,fs,f,s)
	external f,fs
 	double precision a,b,f,s,h,ga,gb,t1,t2,x,g,s0
	n=1
	h=(b-a)/2
	call getgx(a,fs,f,eps,ga)
	call getgx(b,fs,f,eps,gb)
	t1=h*(ga+gb)
10	t2=0.0
      do i=1,n
        x=a+(2*i-1)*h
	  call getgx(x,fs,f,eps,g)
	  t2=t2+g
	enddo
	t2=t1*0.5+h*t2
	s=(4*t2-t1)/3.0
	if(abs(s-s0).ge.eps*(1+abs(s)))then
	  s0=s
	  t1=t2
	  h=h/2.0
	  n=n+n
	  goto 10
	endif
	end subroutine

	subroutine getgx(x,fs,f,eps,g)
	double precision x,f,g,hh,y2,y1,tt1,tt2,y,g0
	call fs(x,y1,y2)
	n=1
	hh=(y2-y1)/2.0
	tt1=hh*(f(x,y1)+f(x,y2))
20	tt2=0.0
	do i=1,n
	  y=y1+(2*i-1)*hh
	  tt2=tt2+f(x,y)
	enddo
	tt2=tt1*0.5+tt2*hh
	g=(4*tt2-tt1)/3.0
	if(abs(g-g0).ge.eps*(1+abs(g)))then
	  g0=g
	  tt1=tt2
	  hh=hh/2.0
	  n=n+n
	  goto 20
	endif
	end subroutine

	subroutine fs(x,y1,y2)
	double precision x,y1,y2
	y2=sqrt(1.0-x*x)
	y1=-sqrt(1.0-x*x)
	end subroutine

	function f(x,y)
	double precision x,y
	f=exp(x*x+y*y)
	end function
	


