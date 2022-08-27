//This code is used to simulate the 2D Ising model by Monte-Carlo algorithm 
//author: JYNi 2020/10/1
program exam1
parameter(L=40,neg=20000,ntal=500)  !改变体系的尺寸和模拟的MCS。
parameter(cm=ntal-neg)
dimension s(L,L)
integer*4::s,i,j,n,k
integer*4::n1,b
real*8::summ,t,Mav,Mas,y(7),sui,h
real*8::expde1(0:4),expde2(0:4),c2,M,x,ii
call random_seed()
c2=1.0*L*L
t=0.1            !模拟的初始温度
h=0.0       !是否加磁场，0表示不加磁场。
do while(t<=5.0)
  call inial(s,L)
  Mav=0.0;Mas=0.0	
  summ=0.0
  do n1=0,4,2
    expde1(n1)=exp(-22.0*(real(n1)+h)/t)
  enddo
  do n1=0,4,2
    expde2(n1)=exp(-22.0*(real(n1)-h)/t)
  enddo
    t=t+0.1
  do b=1,ntal,1
    do n=1,L*L
	  do k=1,7,1
	    call random_number(ii)
		y(k)=ii
	  enddo
	  i=int(1.0*L*y(1)+1.0)
	  j=int(1.0*L*y(4)+1.0)
	  sui=y(7)
	  call mcp(expde1,expde2,i,j,L,s,sui)
	enddo
	  if(b>neg) then
		call sum0(L,s,summ)
		Mav=Mav+abs(summ)/cm
		Mas=Mas+summ*summ/cm
	  endif
    enddo
  M=Mav/c2
  x=(Mas-Mav*Mav)/c2/t
  open(9,file="DATA.TXT",FORM="formAtted",access="sequential")		   
  write(9,100)t,M,x
  100 format(30F30.7)				
  print*,t,M,x
enddo
end
!对每个自旋初始化
subroutine inial(s,L)
  integer*4::s(L,L)
  integer L
  do I=1,L
    do J=1,L
	  s(I,J)=1
    enddo
  enddo									 
end
!每次翻转的函数
subroutine mcp(expde1,expde2,i,j,L,s,sui)
      integer*4::s(L,L),i,j,tempe
	  real*8::expde1(0:4),expde2(0:4),sui
	  tempe=s(i,j)*(s(mod(L+i-2,L)+1,j)+s(mod(i,L)+1,j)+s(i,mod(L+j-2,L)+1)+s(i,mod(j,L)+1))
	  if(s(i,j)==1) then
		if(tempe>=0) then
			if(sui<=expde1(tempe))  s(i,j)=-s(i,j)		
		else
			s(i,j)=-s(i,j)
		endif
	  else
		if(tempe>=0) then
			if(sui<=expde2(tempe))  s(i,j)=-s(i,j)
		else
			s(i,j)=-s(i,j)
		endif
	  endif
end
!计算每次翻转的自旋
subroutine sum0(L,s,summ)
  integer*4::s(L,L)
  real*8::summ
		summ=0.0
	    do i=1,L,1
	      do j=1,L,1
			summ=summ+1.0*s(i,j)
		  enddo
	    enddo
end
