*********************************************************************
*  ���񵴺��������(Integration of Highly-Oscillating-Function)	
*  b,a�����ֱ�Ϊ���������ޡ�
*  aa����Ϊ˫����ʵ�����飬���f(a)��k�׵�����k=0,1,2,��n-1��
*  bb����Ϊ˫����ʵ�����飬���f(b)��k�׵�����k=0,1,2,��n-1��
*  ff����˫����ʵ�ͺ����ӳ�����������sin(k��/2+x)��cos(k��/2+x)��
*  s1,s2�������ظ��񵴺�������S1(m)��S2(m)�Ļ���ֵ��
*  д��1��
**********************************************************************
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

****************************************************************
*  �������ܣ�����sin(k��/2+x)��cos(k��/2+x)��
*  ��l=0ʱ����sin(k��/2+x)����l=1ʱ����cos(k��/2+x)��
*  xΪ˫����ʵ�ͱ�����
****************************************************************
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
	
************************************************************************
*   д��2��
************************************************************************
	subroutine I_HOF(a,b,n,m,aa,bb,s1,s2)
	dimension aa(n),bb(n),sa(4),sb(4),ca(4),cb(4)
	double precision a,b,aa,bb,s1,s2,sa,sb,ca,cb
	sma=sin(m*a)
	smb=sin(m*b)
	cma=cos(m*a)
	cmb=cos(m*b)
	sa(1)=sma
	sa(2)=cma
	sa(3)=-sma
	sa(4)=-cma
	sb(1)=smb
	sb(2)=cmb
	sb(3)=-smb
	sb(4)=-cmb
	ca(1)=cma
	ca(2)=-sma
	ca(3)=-cma
	ca(4)=sma
	cb(1)=cmb
	cb(2)=-smb
	cb(3)=-cmb
	cb(4)=smb
	s1=0.0
	s2=0.0
	mm=1.0
	do k=1,n
	  l=k
20	  if(l.ge.5)then
	    l=l-4
	    goto 20
	  endif
	  mm=mm*m
	  s1=s1+(bb(k)*sb(l)-aa(k)*sa(l))/mm
	  s2=s2+(bb(k)*cb(l)-aa(k)*ca(l))/mm
	enddo
	s2=-s2
	end subroutine
	  