*************************************************************
*  拉艾尔米特-高斯求积法子程序测试程序
*************************************************************	
	external f
	double precision a,f,r
	call HGI(f,r)
	write(*,*) r
	end
	
	subroutine HGI(f,r)
	dimension t(5),coef(5)
	double precision t,coef,f,r,x
	data t/-2.02018200,-0.95857190,0.0,0.95857190,2.02018200/
	data coef/1.181469599,0.9865791417,0.9453089237,0.9865791417,
    &          1.181469599/ 
	r=0.0
	do i=1,7
	  x=t(i)
	  r=r+coef(i)*f(x)
	enddo
	end subroutine

	function f(x)
	double precision x
	f=x*x*exp(-x*x)
	end function