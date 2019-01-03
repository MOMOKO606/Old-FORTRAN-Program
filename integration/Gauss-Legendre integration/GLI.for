********************************************************************
*  ��˹���õ»���(Gauss-Legendre Integration)
*  b,a�����ֱ�Ϊ���������ޡ�
*  eps����Ϊ����Ҫ��
*  f����˫����ʵ�ͺ����ӳ�������������������ڼ��㱻������f(x)��
*  r�������ػ���ֵ��
********************************************************************
*д��1��
********************************************************************
	subroutine GLI(a,b,eps,f,r)
	dimension t(5),coef(5)
	double precision a,b,aa,bb,f,r,h,x,q,temp
	data t/-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459/
	data coef/0.2369268851,0.4786286705,0.5688888889,
     &          0.4786286705,0.2369268851/
	n=1
	temp=(b-a)*(f(a)+f(b))/2.0
10	h=(b-a)/n
	r=0.0
	do k=1,n
	  aa=a+(k-1)*h
	  bb=a+k*h
	  q=0.0
	  do i=1,5
	    x=((bb-aa)*t(i)+(bb+aa))/2.0
		q=q+f(x)*coef(i)
	  enddo
	  r=r+q
	enddo
	r=r*h/2.0
	if(abs(r-temp).ge.eps)then
	  n=n+1
	  temp=r
	  goto 10
	endif
	end subroutine
********************************************************************
*д��2��
********************************************************************
	subroutine GLI(a,b,eps,f,r)
	dimension t(5),coef(5)
	double precision a,b,aa,bb,f,r,h,x,q,temp,s,w
	data t/-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459/
	data coef/0.2369268851,0.4786286705,0.5688888889,
     &          0.4786286705,0.2369268851/
	m=1
	s=(b-a)*0.001
	temp=0.0
10	h=(b-a)/m
	r=0.0
	do k=1,m
	  aa=a+(k-1)*h
	  bb=a+k*h
	  q=0.0
	  do i=1,5
	    x=((bb-aa)*t(i)+(bb+aa))/2.0
		q=q+f(x)*coef(i)
	  enddo
	  r=r+q
	enddo
	r=r*h/2.0
	w=abs(r-temp)/(1.0+abs(r))
	if((w.ge.eps).and.(abs(h).gt.abs(s)))then
	  temp=r
	  m=m+1
	  goto 10
	endif
	end subroutine