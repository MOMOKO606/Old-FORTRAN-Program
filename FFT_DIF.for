*****************************************************************************************
*	时间抽选法（DIF）：Decimation In Frequency
*     I:pr,pi 输入数据的实部、虚部
*     P: 1、计算WN
*        2、计算前k-1级蝶形运算
*        3、计算最后一级蝶形运算
*        4、倒位序转为自然序
*     O:fr,fi 输出数据的实部、虚部
*	
*     n：数据的总点数
*     k: n=2**k	
*     l:l等于0时做Fourier正变换，l不等于零时做Fourier反变换
******************************************************************************************
	subroutine FFT_DIF(pr,pi,n,k,fr,fi,l)
	dimension pr(n),pi(n),fr(n),fi(n)
	double precision pr,pi,fr,fi,p,q,s,poddr,poddi,vr,vi
*     real(8) pr(n),pi(n),fr(n),fi(n)
*	real(8) p,q,s,vr,vi,poddr,poddi
	fr(1)=1.0
	fi(1)=0.0
	fr(2)=cos(6.283185306/n)
	fi(2)=-sin(6.283185306/n)
	if(l.ne.0) fi(2)=-fi(2)
	do  i=3,n
	  p=fr(i-1)*fr(2)
	  q=fi(i-1)*fi(2)
	  s=(fr(i-1)+fi(i-1))*(fr(2)+fi(2))
	  fr(i)=p-q
	  fi(i)=s-p-q
	enddo
	nv=2*n
	do lo=0,k-2
	  m=2**lo
	  nv=nv/2
	  do it=0,(m-1)*nv,nv
	    do j=0,(nv/2)-1
	      poddr=pr(it+j+1)
		  poddi=pi(it+j+1)
		  pr(it+j+1)=poddr+pr(it+j+1+nv/2)
		  pi(it+j+1)=poddi+pi(it+j+1+nv/2)
		  pr(it+j+1+nv/2)=poddr-pr(it+j+1+nv/2)
		  pi(it+j+1+nv/2)=poddi-pi(it+j+1+nv/2)
		  p=pr(it+j+1+nv/2)*fr(m*j+1)
		  q=pi(it+j+1+nv/2)*fi(m*j+1)
		  s=(pr(it+j+1+nv/2)+pi(it+j+1+nv/2))
		  s=s*(fr(m*j+1)+fi(m*j+1))
		  pr(it+j+1+nv/2)=p-q
		  pi(it+j+1+nv/2)=s-p-q
		enddo
	  enddo
	enddo
	do it=0,n-2,2
	  vr=pr(it+1)
	  vi=pi(it+1)
	  pr(it+1)=vr+pr(it+2)
	  pi(it+1)=vi+pi(it+2)
	  pr(it+2)=vr-pr(it+2)
	  pi(it+2)=vi-pi(it+2)
	enddo
	do it=0,n-1
	  m=it
	  is=0
	  do i=0,k-1
	    j=m/2
	    is=is*2+(m-j*2)
*         is=is+(m-j*2)*2**(k-1-i)
		m=j
	  enddo
	fr(is+1)=pr(it+1)
	fi(is+1)=pi(it+1)
	enddo
	if(l.ne.0)then
	  do i=1,n
	    fr(i)=fr(i)/n
	    fi(i)=fi(i)/n
	  enddo
	endif
	end subroutine  


