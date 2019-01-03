*******************************************************************************
*  变步长辛普森二重积分子程序VSSI_2D（Variable Step Simpson Integration_2D）
*  a,b——外层积分的上下限。
*  fs——子程序名，计算内层积分的上下限。
*  f——双精度实型函数子程序名，输入参数，用于计算被积函数f(x)。
*  eps——为精度要求。
*  s——返回积分值。
*******************************************************************************
*  写法1：
*******************************************************************************
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
********************************************************************************
*  写法2：
********************************************************************************
	subroutine VSSI_2D(a,b,eps,fs,f,s)
	external fs,f
	double precision a,b,f,s,h,ga,gb,t1,t2,x,g,s0,c
	n=1
	h=(b-a)/2.0
	c=(b-a)*1.0e-06
	call getgx(a,fs,f,eps,ga)
	call getgx(b,fs,f,eps,gb)
	t1=h*(ga+gb)
10	t2=t1*0.5
      x=a-h
      do i=1,n
        x=x+2*h
	  call getgx(x,fs,f,eps,g)
	  t2=t2+h*g
	enddo
	s=(4*t2-t1)/3.0
	n=n+n
	if(n.ge.16)then
	  if(abs(s-s0).le.eps*(1.0+abs(s)))return
	endif
	s0=s
	t1=t2
	h=h/2.0
	if(abs(h).ge.abs(c))goto 10  
	return
	end subroutine

	subroutine getgx(x,fs,f,eps,g)
	double precision x,f,g,h,y2,y1,tt1,tt2,y,g0,c
	call fs(x,y1,y2)
	n=1
	h=(y2-y1)/2.0
	c=(y2-y1)*1.0e-06
	tt1=h*(f(x,y1)+f(x,y2))
20    tt2=0.5*tt1
	y=y1-h
	do i=1,n
	  y=y+2*h
	  tt2=tt2+f(x,y)
	enddo
	tt2=tt1*0.5+tt2*h
	g=(4*tt2-tt1)/3.0
	n=n+n
	if(n.ge.16)then
	  if(abs(g-g0).le.eps*(1.0+abs(g)))return
	endif
	g0=g
	tt1=tt2
	h=h/2.0
	if(abs(h).ge.abs(c))goto 20
	return
	end subroutine