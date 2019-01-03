*****************************************************************************
计算三种滚雪球方式的定期存款的实际收益以及收益率。
方式一：1个月1存，储蓄存期1年，一年存12次，利率2.25%。
方式二：1个月1存，储蓄存期3个月，一年存12次，利率0.428%。
方式三：3个月1存，储蓄存期3个月，一年存4次，利率0.428%。
参数说明：
year：总共预计存款年限。
pr（principal）：本金。
pattern：方式选择。
ir（interest rate）：银行利率。
lx1（“利息”拼音首字母）：实得月利息。
lx2（“利息”拼音首字母）：实得年利息。
rair（real annual interest rate）:实际年收益率。
sum1：到期本息总和。
sum2：到期利息总和。
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
	
	