*****************************************************************************
�������ֹ�ѩ��ʽ�Ķ��ڴ���ʵ�������Լ������ʡ�
��ʽһ��1����1�棬�������1�꣬һ���12�Σ�����2.25%��
��ʽ����1����1�棬�������3���£�һ���12�Σ�����0.428%��
��ʽ����3����1�棬�������3���£�һ���4�Σ�����0.428%��
����˵����
year���ܹ�Ԥ�ƴ�����ޡ�
pr��principal��������
pattern����ʽѡ��
ir��interest rate�����������ʡ�
lx1������Ϣ��ƴ������ĸ����ʵ������Ϣ��
lx2������Ϣ��ƴ������ĸ����ʵ������Ϣ��
rair��real annual interest rate��:ʵ���������ʡ�
sum1�����ڱ�Ϣ�ܺ͡�
sum2��������Ϣ�ܺ͡�
*****************************************************************************
	program fixed_deposit
	character *40 cmdfile,outputfile
	integer year,pr,pattern
	real ir,sum1,sum2
	real,allocatable::lx1(:),lx2(:),rair(:)
	call read_cmd(cmdfile,year,pr,pattern,ir,outputfile)
	allocate(lx1(1:year),lx2(1:year),rair(1:year))
	call calulate_result(pattern,year,pr,ir,lx1,lx2,rair,
     &sum1,sum2)
	call output_result()
	deallocate()
	end program

	subroutine read_cmd(cmdfile,year,pr,pattern,ir,outputfile)
	character*(*) cmdfile,outputfile
	integer year,pr,pattern
	real ir
	cmdfile='cmd.txt'
	open(11,file=cmdfile,status='old')
	read(11,*) pattern
	read(11,*) year
	read(11,*) pr
	if (pattern.eq.1) then
	  ir=0.0225
	else
	  ir=0.004275
	endif
	read(11,*) outputfile
	close(11)
	endsubroutine 	

	subroutine calulate_result(pattern,year,pr,ir,lx1,lx2,rair,
     &sum1,sum2)
	real temp
	real lxx(1:year,1:4)
	temp=0.0
	if (pattern.eq.1)then
	  do i=1,year,1
	    lx1(i)=(temp+pr*i)*ir*0.95
	    lx2(i)=lx1(i)*12
	    temp=temp+lx1(i)
	    rair(i)=temp/pr*i
	  enddo
	else	  
	  do i=1,year,1
	    do j=1,4,1
	        lxx(i,j)=temp+pr*((i-1)*4+j)*ir*0.95
		    temp=temp+lxx(i,j)
	    enddo
	    if(pattern.eq.2)then
		  lx1(i)=temp/(i*4)
	      lx2(i)=(temp*3)/i
	    else
	      lx1(i)=temp/(i*12)
	      lx2(i)=temp/i
		endif
	      rair(i)=temp/(i*4*pr)
	  enddo
	endif    
	if(pattern.ne.3)then
	  sum1=pr*year*12+temp
	else
	  sum1=pr*year*4+temp 
	endif  
	  sum2=temp
	endsubroutine  

	subroutine output_result(year,lx1,lx2,rair,sum1,sum2,outputfile)
	
	