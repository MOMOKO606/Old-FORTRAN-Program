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

	  
