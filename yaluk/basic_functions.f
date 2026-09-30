
!**************************************************************************************************
!*  BASIC FUNCTIONS and SUBROUTINES
!--------------------------------------------------------------------------------------------------
!*	FUNCTION								
!*  flipud(W,Z,nn,gg)									: Flip vector up and down
!*  mult (W,X,Z,nn,mm,zz)								: Matrix multiplication
!*  inv(B,n,np,A)										: Matrix inverse (Decomposition LU)
!*  ludcmp(a,n,np,indx,d)								: LU Decomposition
!*  lubksb(a,n,np,indx,b)								: Backsubstitution subroutinte
!*  INTEGRAL(FUN,xlow,xhigh,integ,tol,gn1)				: Modified Quadrature Gauss integral
!*  INTEGRAL2(FUN,xlow,xhigh,integ,tol,gn1,var1,var2)	: Modified Quadrature Gauss integral
!*  GAUSSQJ(FUN,xlow,xhigh,gn,int)						: Gauss coefficients 
!*  BESSELI0(K1)										: Bessel Integral 0
!*  BESSELI1(K1)										: Bessel Integral 1
!*  comparar(A,B,cond)									: Comparison Function
!*
!**************************************************************************************************


!*	----------------------------------------
!*	FLIPUD Flip matrix in up/down direction
!*	------------------------------------------	
	SUBROUTINE flipud(W,Z,nn,gg)
	INTEGER nn,gg,i,j
	DOUBLE PRECISION W(nn,gg),Z(nn,gg)
    DO 42 j=1,gg
	    DO 43 i=1,nn
			Z(i,j)=W(nn-i+1,j)
43	    CONTINUE
42	CONTINUE
	RETURN
	END



!Cm			----------------------------------------
!Cm				MATRIX MULTIPLICATION
!Cm			----------------------------------------
	SUBROUTINE mult(W,X,Z,nn,mm,zz)
	INTEGER nn,mm,zz,i,j,k
	DOUBLE PRECISION W(mm,nn),X(nn,zz),Z(mm,zz)
	
	DO 63 i=1,mm
		DO 62 j=1,zz
			Z(i,j)=0.D0
			DO 61 k=1,nn
				Z(i,j)=W(i,k)*X(k,j)+Z(i,j)
61			CONTINUE
62		CONTINUE
63	CONTINUE
	RETURN
	END



!Cm			-----------------------------------
!Cm					INVERSION MATRIX
!Cm			   	 USES ludcmp AND lubksb
!Cm			-----------------------------------
	SUBROUTINE inv(B,n,np,Ai)
	INTEGER n,np,indx(n),i,j
	DOUBLE PRECISION A(np,np),E,Ai(np,np),B(np,np)
	A(:,:)=B(:,:)
	
	DO 51 j=1,n
	    DO 50 i=1,n
			Ai(i,j)=0.D0
50		CONTINUE
		Ai(j,j)=1
51	CONTINUE
	CALL ludcmp(A,n,np,indx,E)
	
		do 52 j=1,n
			call lubksb(A,n,np,indx,Ai(1,j))
52		CONTINUE
	
	RETURN

	END



!Cm			---------------------------
!Cm			      LU DECOMPOSITION
!Cm			---------------------------
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (TINY=1.0e-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(n)
      d=1.
	  
      do 12 i=1,n
        aamax=0.D0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
		
        if (aamax.eq.0.) THEN 
			write(*,*) '******************************************************************'
			write(*,*) '*ERR* Singular Matrix in ludcmp Decomposition'
            write(*,*) '*ERR* It may be two conductors in the same position, or heigth equal to zero'
            write(*,*) '*PROGRAM WILL STOP*'
            write(*,*) '******************************************************************'
            STOP
        ENDIF
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.D0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END	



!Cm			---------------------------
!Cm			BACKSUBSTITUTION ROUTINE
!Cm			---------------------------

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END


!*																		*
!*																		*
!************************************************************************
	
!************************************************************************
!*	MODIFIED INTEGRATION SUBROUTINE										*
!*	*******************************										*
	SUBROUTINE INTEGRAL(FUN,xlow,xhigh,integ,tol,gn1)

	IMPLICIT NONE
	INTEGER gn, i, converge, nk, gn1
	Real(8) tol, intold, tol1, xlow,xhigh, integ, mult 
	EXTERNAL FUN

!*---------------------------------
!*	Definición Número de Puntos
	mult=1.7
	IF (gn1 .GT. 7) THEN
	mult=1.7
	gn=ceiling(gn1/mult)
   	ELSE
	gn=4
	mult=1.7
	END IF
	IF (gn1 .GT. 20) THEN
	mult=1.5
	gn=ceiling(gn1/mult)
	END IF
	IF (gn1 .GT. 200) THEN
	gn=20
	mult=1.5

	END IF

!*------------------------------------
!*Ciclo d Integración

	CALL GAUSSQJ(FUN, xlow, xhigh, gn, integ)
	
	!CALL trapzd(FUN,xlow,xhigh,integ,gn)
	intold=integ
	integ = 0
	converge = 0

	DO  i=1,10
		gn=int(gn*mult)
		CALL GAUSSQJ(FUN, xlow, xhigh, gn, integ)
	
		!CALL trapzd(FUN,xlow,xhigh,integ,gn)
!WRITE(*,*) integ
		tol1 = ABS(intold-integ) 
  
		IF (tol1 .GT. ABS(tol*integ)) THEN
			nk = 1
		ELSE
			converge = 1
	        gn1=gn
	        gn=0
			EXIT
		    
		END IF
	intold=integ
	integ=0
	gn1=gn
	END DO

	IF (converge .EQ. 0) THEN
	WRITE(*,*)'No convergencia de la integral',tol1/integ
	END IF

	RETURN

	END
	
	SUBROUTINE GAUSSQJ(FUN,xlow,xhigh,gn,int)

	IMPLICIT NONE
	INTEGER gn, i
	DOUBLE PRECISION xlow, xhigh, integ, jacob, x(gn), cb(gn),	&
     cw(gn), y(gn), int, FUN
	 EXTERNAL FUN
	int=0.D0
	integ=0.D0
	CALL QRULE(gn,cb,cw)
      
	jacob=(xhigh-xlow)/2

	x = ((cb+1) * jacob + xlow )

	!$OMP SIMD REDUCTION(+:int)	
	DO i=1,gn
		y(i) = FUN(x(i))
          integ=(cw(i)*y(i))*jacob
	    int = int + integ
    END DO		

	
	RETURN

	END

!***********************************************************************
! INTEGRAL GUNCTION FOR EVALUATION WITH TWO MORE ARGUMENTS VAR1, VAR2 
! MADE ESPETIALLY FOR INCLUDING THE INTEGRAL CURRENT AND TIME
!***********************************************************************
		SUBROUTINE INTEGRAL2(FUN,xlow,xhigh,integ,tol,gn1,var1,var2)

	IMPLICIT NONE
	INTEGER gn, i, converge, nk, gn1
	Real(8) tol, intold, tol1, xlow,xhigh, integ, mult,var1(*),var2(*)
	EXTERNAL FUN

!*---------------------------------
!*	Definición Número de Puntos
	mult=1.7
	IF (gn1 .GT. 7) THEN
	mult=1.7
	gn=ceiling(gn1/mult)
   	ELSE
	gn=4
	mult=1.7
	END IF
	IF (gn1 .GT. 20) THEN
	mult=1.5
	gn=ceiling(gn1/mult)
	END IF
	IF (gn1 .GT. 200) THEN
	gn=20
	mult=1.5

	END IF
!*------------------------------------
!*Ciclo d Integración

	CALL GAUSSQJ2(FUN, xlow, xhigh, gn, integ,var1,var2)
	
	!CALL trapzd2(FUN,xlow,xhigh,integ,gn,var1,var2)

	intold=integ
	integ = 0
	converge = 0

	DO  i=1,10
		gn=gn*mult  
		CALL GAUSSQJ2(FUN, xlow, xhigh, gn, integ,var1,var2)

	!CALL trapzd2(FUN,xlow,xhigh,integ,gn,var1,var2)

		tol1 = ABS(intold-integ) 
  
		IF (tol1 .GT. ABS(tol*integ)) THEN
			nk = 1
		ELSE
			converge = 1
	        gn1=gn
	        gn=0
			EXIT 
		    
		END IF
	intold=integ
	integ=0.D0
	gn1=gn
   END DO

	IF (converge .EQ. 0) THEN
	WRITE(*,*)'No convergencia de la integral2',tol1/integ
	END IF

	RETURN

	END
	
	SUBROUTINE GAUSSQJ2(FUN,xlow,xhigh,gn,int,var1,var2)

	IMPLICIT NONE
	INTEGER gn, i
	DOUBLE PRECISION xlow, xhigh, integ, jacob, x(gn), cb(gn),	&
     cw(gn), FUN, y(gn), int,var1(*),var2(*)
	int=0.D0
	integ=0.D0
	CALL QRULE(gn,cb,cw)
      
	jacob=(xhigh-xlow)/2

	x = ((cb+1) * jacob + xlow )

	!$OMP SIMD REDUCTION(+:int)
	DO i=1,gn
		y(i) = FUN(x(i),var1,var2)
          integ=(cw(i)*y(i))*jacob
	    int = int + integ
	END DO
	RETURN

	END


!*--------------------------------------
!*Subrutina para el cálculo de raices de integración
			
	SUBROUTINE QRULE(n,bp,wf)
	
	IMPLICIT NONE
	INTEGER m,n,iter,n1j,jj,i, k 
	DOUBLE PRECISION e1, pi, t((n+1)/2), xo((n+1)/2), bp(n), wf(n),		&
     nn,j, pkm1((n+1)/2), pk((n+1)/2), t1((n+1)/2), pkp1((n+1)/2),   	&
     den((n+1)/2),d1((n+1)/2),dpn((n+1)/2),d2pn((n+1)/2),d3pn((n+1)/2),	&
     d4pn((n+1)/2),u((n+1)/2),v((n+1)/2),Hcan((n+1)/2),p((n+1)/2),		&
     dp((n+1)/2),fx((n+1)/2)
	PARAMETER (pi = 3.141592653589793D0)

	iter=2
	e1=n*(n+1)
	m=int((n+1)/2)
			
	j=3D0
	DO 10, i=1,m
	 t(i)=(pi/(4D0*n+2D0))*j
	 j=j+4D0
10	CONTINUE		

	nn=1D0-(1D0-1D0/n)/(8D0*n*n)
	
	xo=nn*cos(t)

	DO 20, i=1,iter
		pkm1=1 
		pk=xo
			DO 30, k=2,n
				t1=xo*pk
				pkp1=t1-pkm1-(t1-pkm1)/k+t1
				pkm1=pk 
				pk=pkp1
30			CONTINUE
		den=1-xo*xo 
		d1=n*(pkm1-xo*pk)
		dpn=d1/den
		d2pn=(2D0*xo*dpn-e1*pk)/den
		d3pn=(4*xo*d2pn+(2D0-e1)*dpn)/den
		d4pn=(6*xo*d3pn+(6D0-e1)*d2pn)/den
		u=pk/dpn 
		v=d2pn/dpn
		Hcan=-u*(1D0+(.5*u)*(v+u*(v*v-u*d3pn/(3D0*dpn))))
		p=pk+Hcan*(dpn+(.5*Hcan)*(d2pn+(Hcan/3D0)*(d3pn+.25*Hcan*d4pn)))
		dp=dpn+Hcan*(d2pn+(.5*Hcan)*(d3pn+Hcan*d4pn/3D0))
		Hcan=Hcan-p/dp 
		xo=xo+Hcan
20	CONTINUE

	bp(1:(n+1)/2)=-xo-Hcan
	fx=d1-Hcan*e1*(pk+(Hcan/2)*(dpn+(Hcan/3)*(d2pn+(Hcan/4)*(d3pn+(.2*Hcan)*d4pn))))
	wf(1:(n+1)/2)=2*(1-bp**2)/(fx*fx)

	IF ((m+m).GT. n) THEN
	 bp(m)=0.D0
	END IF

	IF (.NOT.((m+m) == n)) THEN
	 m=m-1
	END IF
	
	DO 40, jj=1,m 
		n1j=(n+1-jj) 
		bp(n1j)=-bp(jj)
		wf(n1j)=wf(jj)
40	CONTINUE

	RETURN
	END

	DOUBLE PRECISION FUNCTION BESSELI0(K1)

	IMPLICIT NONE

	DOUBLE PRECISION K, i1, K1

	IF (K1 .LE. 3.75) THEN
		
		K=K1/3.75
		i1 = 1 + 3.5156229 * K**2 + 3.0899424 * K**4	&
     		+ 1.2067492 * K**6 + 0.2659732 * K**8		&
     		+ 0.0360768 * K**10 + 0.0045813 * K**12		
		BESSELI0 = i1 * DEXP(-K1)
   
	ELSE
		K=K1/3.75
		i1 = 0.39894228 + 0.01328592 * K**(-1)			&
     		+ 0.00225319 * K**(-2) - 0.00157565 * K**(-3)	&
     		+ 0.00916281 * K**(-4) - 0.02057706 * K**(-5)	&
     		+ 0.02635537 * K**(-6) - 0.01647633 * K**(-7)	&
     		+ 0.00392377 * K**(-8)
		BESSELI0 = i1/DSQRT(K1)

		END IF

	END

	DOUBLE PRECISION FUNCTION BESSELI1(K1)

	IMPLICIT NONE

	DOUBLE PRECISION K, i1, K1

	IF (K1 .LE. 3.75) THEN
		K=K1/3.75
		i1 = 0.5 + 0.87890594 * K**2 + 0.51498869 * K**4	&
     		+ 0.15084934 * K**6 + 0.02658733 * K**8			&
     		+ 0.00301532 * K**10 + 0.00032411 * K**12
		BESSELI1 = i1 * DEXP(-K1) * K1      
	ELSE
		K=K1/3.75
		i1 = 0.39894228 - 0.03988024 * K**(-1)				& 
     		- 0.00362018 * K**(-2) + 0.00163801 * K**(-3)		&
     		- 0.01031555 * K**(-4) + 0.02282967 * K**(-5)		&
     		- 0.02895312 * K**(-6) + 0.01787654 * K**(-7)		&
     		- 0.00420059 * K**(-8)
			BESSELI1 = i1/DSQRT(K1)
	END IF
	END
!* ----------------------------------------------
!Cm     INTERPOLATION FOR FUNCTIONS
!*----------------------------------------------
      SUBROUTINE interp(vt,vx,size,tr,xr)
      INTEGER size,  inds, indi
      DOUBLE PRECISION vt(size),vx(size),tr,xr,dt,ind1
      
	dt=vt(2)-vt(1)
	ind1=(tr+dt)/dt
	indi=int(ind1)
	inds=indi+1
	xr=vx(indi)+(vx(inds)-vx(indi))*(ind1-indi)
	RETURN
	END

!
!-------------------------------------------

	SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END

	SUBROUTINE trapzd2(func,a,b,s,n,var1,var2)
      INTEGER n
      DOUBLE PRECISION a,b,s,func,var1(*),var2(*)	
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a,var1,var2)+func(b,var1,var2))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x,var1,var2)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
	 write(*,*) 's=',s,'  n=',n
      return
      END


	
!******************************************************************************************
!										COMPARAR VECTOR
!******************************************************************************************
LOGICAL function comparar(A,B,cond)
	INTEGER i,cond
	DOUBLE PRECISION A(cond),B(cond)
	
	LOGICAL resultado
		DO i=1,cond
			IF (A(i).EQ.B(i)) THEN
			    comparar =.TRUE.
			ELSE
			    comparar =.FALSE.
			    EXIT	
			ENDIF
		END DO

end function comparar

