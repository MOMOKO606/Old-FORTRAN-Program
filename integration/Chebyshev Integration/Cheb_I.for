*************************************************************
*  �б�ѩ�������(Chebyshev Integration)
*  b,a�����ֱ�Ϊ���������ޡ�
*  eps����Ϊ����Ҫ��
*  f����˫����ʵ�ͺ����ӳ�������������������ڼ��㱻������f(x)��
*  s�������ػ���ֵ��
*************************************************************
	subroutine Cheb_I(a,b,eps,f,s)	
	dimension t(5)
	double precision a,b,aa,bb,t,h,r,x,temp,f,s
	data t/-0.8324975,-0.3745414,0.0,0.3745414,0.8324975/
	n=1
10	h=(b-a)/n
	s=0.0
	do k=1,n
	  aa=a+(k-1)*h
	  bb=a+k*h
	  r=0.0
	  do i=1,5
		x=((bb-aa)*t(i)+(bb+aa))/2.0
		r=r+f(x)
	  enddo
	  r=r*(bb-aa)/5.0
	  s=s+r
	enddo
	if(abs(s-temp).ge.eps)then
	  temp=s
	  n=n+1
	  goto 10
	endif
	end subroutine
	  
