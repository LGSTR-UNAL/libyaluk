!***************************************************************************
!***************************************************************************
!*	ELECTROMAGNETIC FIELD COMPUTATION SUBROUTINE Ver.-YALUK-ATP			   *
!*
!*
!*
!*
!***************************************************************************
!***************************************************************************
	SUBROUTINE	EMFIELD(SIZE_T, EZ, EH,ang,Icorr,tcorr,i02,i12,kcampo,t0_min)
	IMPLICIT NONE
	
	INTEGER i,SIZE_T,j, SIZE_1,SIZE_2, gn1v,gn2v,gn3v,gn1h,gn2h,		&
      gn3h,gn4h,gn5h,gn6h,gn1b,gn2b, conductividad, 	&
      indi,inds

	DOUBLE PRECISION c,v, zlam, r, h,Ih1, Ih2, tao11, tao12, tao21, &
     tao22, n1, n2, eta1, eta2,t,ang,e0,dx,								&
     EIRZ(SIZE_T,1), EELZ(SIZE_T,1), EZ2(SIZE_2,1),EZ(SIZE_T,1),			&
     B(SIZE_T,1), BRAD(SIZE_T,1),							&
     EIRR1, EELR1,EIRR2, EELR2, EH(SIZE_T,1),EH2(SIZE_2,1),			&
     EIRR(SIZE_T,1), EELR(SIZE_T,1),EX2(SIZE_T,1),	&
     tolh,tolz,tolb,B2(SIZE_2,1),m2(SIZE_T,1),m22(SIZE_T,1),		&
     beta, Hcan, mu, c1, c3, dt,dt2, zini,zfin,zfin2,pi,RSLT, &
     tini,ind,er2,e,o,t2,c4,kcampo,						&
	 Icorr(*),tcorr(*),i02(*),i12(*),trapzd,midexp,chint,t0_min
	
	LOGICAL(4) triangular

	COMMON /COND2/ conductividad
	COMMON /COND/  o,e
	COMMON /LINE_YALUK/  SIZE_1,SIZE_2
	COMMON /ITYPE/ triangular
	!*	----------------------------------------------
	!Cm	COMMON DATA FOR LINE DESCRIPTION
	!Cm		r: Distance for Electromagnetic Field Calculation
	!Cm		h: Line Height
	!Cm		XR: Cordinate X	
	!Cm		YR: Cordinate Y
	!*	---------------------------------------------- 
	COMMON /POINT/ r,h

	!*	----------------------------------------------
	!Cm	COMMON DATA /DAT/ FOR CHANNEL
	!Cm		c: light velocity
	!Cm		v: return stroke velocity
	!Cm		zlam: Atenuation MTL model if zlam=0  TL model
	!Cm		Hcan: Height Channel
	!Cm		dt: Discretization Time
	!*	----------------------------------------------

	COMMON /DAT/   v,zlam,Hcan,dt, dt2

	!*	----------------------------------------------
	!Cm	COMMON /DAT2/ DATA FOR HEIDLER CURRENT WAVEFORM
	!*	-----------------------------------------------
	!Cm		Ih1: Amplitud Corriente (function 1)
	!Cm		Ih2: Amplitud Corriente (function 2)
	!Cm		tao11: Front Time Constant (function 1)
	!Cm		tao21: Decay Time Constant (function 1)
	!Cm		tao12: Front Time Constant (function 2)
	!Cm		tao22: Decay Time Constant (function 2)
	!Cm		n1: Amplitude Constant (function 1)
	!Cm		n2: Amplitude Constant (function 2)
	!*	-----------------------------------------------

	COMMON /DAT2/ Ih1, Ih2, tao11, tao12,tao21,tao22, n1, n2, eta1, eta2
	COMMON /TIME/ t


    PARAMETER (pi=3.141592653589793D0)
 	PARAMETER (mu=(4*pi)*1E-7)
	PARAMETER (c=2.99792458E8)
	PARAMETER (e0=1/(mu*c**2))


	

	EXTERNAL EZP, EELZP, ERP,BINDP, EELRP, IPULSE, DIPULSE,trapzd,midexp,chint
	  
    
!*	________________________________________________
!Cm	Constants Definition
!Cm	beta: relation between v and c
!Cm	mu: permeability vacuum constant
!*	________________________________________________

	
    beta = v / c
	C1 = -mu * r*c / (2D0 * pi)
	c3 = 1. / (4. * pi * e0)
	c4 = c3/c
	!t = -dt
	t = t0_min-dt
!* Field arriving time
	tini = DSQRT(r**2+h**2)/c
!*	----------------------------------------------
!Cm	CONSTANTS FOR CONTROL INTEGRALS
!*	----------------------------------------------
	gn1v=0
	gn2v=0
	gn3v=0
	gn1h=0
	gn2h=0
    gn3h=0
	gn4h=0
	gn5h=0
	gn6h=0
	gn1b=0
	gn2b=0
!*	----------------------------------------------
!Cm	INTEGRATION TOLERANCES
!Cm	tolh: Horizontal Electric Field tolerance
!Cm	tolv: Vertical Electric Field tolerance
!Cm	tolb: Magnetic Field tolerance
!*	----------------------------------------------

	tolh=1.D-3
	tolz=1.D-3
	tolb=1.D-3
!*	----------------------------------------------
!Cm	FIRST INTEGRATION WINDOW
!*	----------------------------------------------


!10 - start of cycle of EM calculation
	DO i = 1,SIZE_1  !10
	t=t+dt
	    IF (t .LE. tini) THEN
			EZ(i,1) = 0.D0
			EELZ(i,1) = 0.D0
			BRAD(i,1) = 0.D0
			B(i,1) = 0.D0
			EELR(i,1) = 0.D0
			EIRR(i,1) = 0.D0
	        EH(i,1)=0.D0
		ELSE
!*	-------------------------------------------
!Cm            INTEGRATION LIMITS CALCULATION
!Cm            Solving analytically the equation
!Cm            t=z/v+R/c
!Cm            z:Channel height
!Cm            v:return stroke velocity
!Cm            R:Distance Calcularion
!Cm            c:light velocity
!Cm
!Cm             x              Line
!Cm         *<------->======================
!Cm		 !  	+
!!Cm		 !	  +   
!Cm		y!	+    R
!Cm		 ! +
!Cm		 *
!*	-------------------------------------------	

           zini=(-c*t-v/c*h+((-v*t-h)**2+(1-(v/c)**2)*r**2)**.5)/(c/v*(1-(v/c)**2))
           zfin2=(c*t-v/c*h-((v*t-h)**2+(1-(v/c)**2)*r**2)**.5)/(c/v*(1-(v/c)**2))
!Cm		Limite de integracion con h=0
	       zfin=(c*t-((v*t)**2+(1-(v/c)**2)*r**2)**.5)/(c/v*(1-(v/c)**2))	

				IF( zfin2 .GT. Hcan ) THEN
       				zfin2 = Hcan
				END IF
				IF( zfin .GT. Hcan ) THEN
       				zfin = Hcan
				END IF
				IF( zini .LT. -Hcan ) THEN
       				zini = -Hcan
				END IF
			IF (kcampo.EQ.1) THEN
				CALL INTEGRAL(EZP,0.D0,zfin,RSLT,tolz,gn1v)
				EIRZ(i,1) =(2*c3) *RSLT

				CALL INTEGRAL2(EELZP,0.D0,zfin,RSLT, tolz,gn3v,Icorr,tcorr)

				EELZ(i,1) = 2*c3*RSLT
				EZ(i,1) = EIRZ(i,1)+EELZ(i,1)
			ELSE
				EZ(i,1) =	0
			END IF

	        IF (conductividad .EQ. 1) THEN
			  CALL INTEGRAL(BINDP,0.D0,zfin,RSLT, tolb,gn1b)
			  B(i,1) = C1*RSLT
!*			  CALL INTEGRAL(BRADP,0.D0,zfin,RSLT, tolb,gn2b)
!*			  BRAD(i,1) = (C1/c)*RSLT
            END IF
			IF (ang.EQ.0) THEN !evaluating if the angle is 0
				EH(i,1) = 0
			ELSE
				CALL INTEGRAL(ERP,zini,0.D0,RSLT,tolh,gn1h)
				EIRR1 = (r*c4)*RSLT
				CALL INTEGRAL2(EELRP,zini,0.D0,RSLT, tolh,gn3h,Icorr,tcorr)
				EELR1 = 3*r*c3*RSLT
				CALL INTEGRAL(ERP,0.D0,zfin2,RSLT,tolh,gn4h)
				EIRR2 = (r*c4)*RSLT
				CALL INTEGRAL2(EELRP,0.D0,zfin2,RSLT,tolh,gn6h,Icorr,tcorr)
				EELR2 = 3*r*c3*RSLT
			
				EIRR(i,1) = EIRR1 + EIRR2
				EELR(i,1) = EELR1 + EELR2
				EH(i,1) = (EELR(i,1)+EIRR(i,1)) !Without taking into account the incidence angle
			END IF		
		END IF
	    B2(1,1) =B(SIZE_1,1)
		EH2(1,1) =EH(SIZE_1,1)	
								  
        EZ2(1,1) =EZ(SIZE_1,1)
	END DO 

!finishing do 10

!*	-----------------------------------------------------------------
!Cm					Calculation Second Window
!Cm	_________________________________________________________________
!Cm	Note:
!Cm	Check "if" the second window is outside of the calculation window
!*	-----------------------------------------------------------------

	IF (SIZE_1 .LT. SIZE_T) THEN	
	  
  
	  DO  j = 2,SIZE_2

	    IF (t .LE. tini) THEN
			EZ(j,1) = 0
			EELZ(j,1) = 0
			BRAD(j,1) = 0
			B2(j,1) = 0
			EELR(j,1) = 0
			EIRR(j,1) = 0

		ELSE
		t = t+dt2
!*	Limites d Integración
			h = ABS(h)
!*			t = t + DSQRT(r**2+h**2)/c
              zini=(-c*t-v/c*h+((-v*t-h)**2+(1-(v/c)**2)*r*r)**.5)	&
     			/(c/v*(1-(v/c)**2))
              zfin2=(c*t-v/c*h-((v*t-h)**2+(1-(v/c)**2)*r*r)**.5)	&
     			/(c/v*(1-(v/c)**2))
!Cm			Limite de integracion con h=0
	        zfin=(c*t-((v*t)**2+(1-(v/c)**2)*r**2)**.5)				&
     			/(c/v*(1-(v/c)**2))

!*			H1 = (beta*(c*t-v/c*h-((beta*c*t-h)**2+(1-beta**2)*r**2)
!*     $		   **0.5))/(1-beta**2)

!*			H2 = (beta*(c*t-((beta*c*t)**2+(1-beta**2)*r**2)
!*     $		   **0.5))/(1-beta**2)
				IF( zfin2 .GT. Hcan ) THEN
       				zfin2 = Hcan
				END IF
		
				IF( zfin .GT. Hcan ) THEN
       				zfin = Hcan
				END IF
				IF( zini .LT. -Hcan ) THEN
       				zini = -Hcan
				END IF
			IF (kcampo.EQ.1) THEN ! Evaluating if it is neccesary to calculate the vertical component

				CALL INTEGRAL(EZP,0.D0,zfin,RSLT,tolz,gn1v)
				EIRZ(SIZE_1+j,1) =(2*c3) *RSLT

				CALL INTEGRAL2(EELZP,0.D0,zfin,RSLT, tolz,gn3v,Icorr,tcorr)
				EELZ(SIZE_1+j,1) = 2*c3*RSLT
				EZ2(j,1) = EIRZ(SIZE_1+j,1)+EELZ(SIZE_1+j,1)
			ELSE
				EZ2(j,1) = 0
			END IF

	        IF (conductividad .EQ. 1) THEN
			  CALL INTEGRAL(BINDP,0.D0,zfin,RSLT,tolb,gn1b)
			  B2(j,1) = C1*RSLT
!*			  CALL INTEGRAL(BRADP,0.D0,zfin,RSLT, tolb,gn2b)
!*			  BRAD(j,1) = (C1/c)*RSLT
            END IF

			IF (ang.EQ.0) THEN !evaluating if the angle is 0
				EH2(j,1) = 0
			ELSE
				CALL INTEGRAL(ERP,zini,0.D0,RSLT,tolh,gn1h)
				EIRR1 = (r*c4)*RSLT

				CALL INTEGRAL2(EELRP,zini,0.D0,RSLT, tolh,gn3h,Icorr,tcorr)
				EELR1 = 3*r*c3*RSLT


				CALL INTEGRAL(ERP,0.D0,zfin2,RSLT,tolh,gn4h)
				EIRR2 = (r*c4)*RSLT
				RSLT = 0

				CALL INTEGRAL2(EELRP,0.D0,zfin2,RSLT, tolh,gn6h,Icorr,tcorr)
				EELR2 = 3*r*c3*RSLT
		
				EIRR(SIZE_1+j,1) = EIRR1 + EIRR2
				EELR(SIZE_1+j,1) = EELR1 + EELR2
				EH2(j,1) = EELR(SIZE_1+j,1)+EIRR(SIZE_1+j,1)
			END IF
		END IF
			
	
	

    END DO
!finishing do 20

        !t=dt2
		t=t0_min+dt2
	  DO 25  i=SIZE_1+1,SIZE_T
	   	
			t=t+dt
			!ind=(t)/dt2
			ind=(t-t0_min)/dt2
			indi=floor(ind)
	   IF (indi .LT. SIZE_2) THEN
		  inds=indi+1
	 	  EH(i,1)=EH2(indi,1)+(EH2(inds,1)-EH2(indi,1))*(ind-indi)
		  IF (kcampo.EQ.1) THEN 
			EZ(i,1)=EZ2(indi,1)+(EZ2(inds,1)-EZ2(indi,1))*(ind-indi)
		  ELSE
			EZ(i,1)=0
		  END IF
	      IF (conductividad .EQ. 1) THEN
			  B(i,1)=B2(indi,1)+(B2(inds,1)-B2(indi,1))*(ind-indi)
		  END IF

	   ELSE
	        EH(i,1)=EH2(indi,1)
	    	EZ(i,1)=EZ2(indi,1)

	   END IF
25	  CONTINUE	

	END IF
!*	----------------------------------------------
!Cm	Calculation of Ground Conductivity Influence
!Cm	Convolution
!*	------------------------------------------------
 	
	    IF (conductividad .EQ. 1) THEN
			IF (ang.NE.0) THEN
				
			
				DO 26 i=1,SIZE_T-1
	        		m2(i,1)=(B(i+1,1)-B(i,1))/dt
26		        CONTINUE
					m2(SIZE_T,1)=m2(SIZE_T-1,1)
				  m22(SIZE_T,1)=0.D0
				DO 27 i=2,SIZE_T
					m22(i,1)=m2(i,1)-m2(i-1,1)
27		        CONTINUE
                
					EX2=0.D0
					t2=-dt
				DO 28 i=2,SIZE_T
				  t2=t2+dt
					EX2(i,1)=(m2(1,1)/e**0.5)*(t2)*(i02(i)+i12(i))
28			    CONTINUE
				  EH(1,1)=EH(1,1)+EX2(1,1)

				DO 65 j=2,SIZE_T
					 t2=-dt
				   DO 45 i=j,SIZE_T
						t2=t2+dt
						er2=(m22(j,1)/e**0.5)*(t2)*(i02(i-j+1)+i12(i-j+1))
					  EX2(i,1)=EX2(i,1)+er2
45				   CONTINUE
					  EH(j,1)=EH(j,1)+EX2(j,1)  
65				CONTINUE
				END IF
	       END IF

		EH(:,1)=EH(:,1)*ang !Taking into account the incidence angle
	END SUBROUTINE

!*																		*
!*																		*
!***************************************************************************
!***************************************************************************
!*	Radiated and Induced Electric Field Horizontal Component 			*
!*	********************************************************			*

	DOUBLE PRECISION FUNCTION EZP(z)
	!$OMP DECLARE SIMD (EZP)
	IMPLICIT NONE
	COMMON /POINT/ r, h
	COMMON /TIME/ t
	COMMON /DAT/   v,zlam,Hcan,dt, dt2

	DOUBLE PRECISION Rz, A, r, h, c , v, z, ic, t, zlam,	&
     IPULSE, Hcan, dt,ERADZP,EINDZP,dic, DIPULSE, dt2

	PARAMETER (c=2.99792458E8)




	EXTERNAL IPULSE, DIPULSE
	Rz = DSQRT(r*r + z*z)
	A = t - (Rz/c + ABS(z)/v)

	IF (A .LE. 0.) THEN
		ic = 0
	    dic = 0
	ELSE
   		ic = IPULSE(A)
	    dic= DIPULSE(A)
	END IF
	IF (zlam.EQ.0) THEN
		ERADZP = -dic/(Rz**3)
		EINDZP = ic*(2*(z*z) - r*r)/(Rz**4)
	ELSE
		ERADZP = -dic*DEXP(-abs(z)/zlam)/(Rz**3)
		EINDZP = ic*(2*(z*z) - r*r)*DEXP(-abs(z)/zlam)/(Rz**4)
	ENDIF
	EZP=(ERADZP*r*r/(c*c)+EINDZP/c)
	END
	
!*																		*
!*																		*
!***************************************************************************
!***************************************************************************
!*	Electrostatic Field Vertical Component 								*
!*	**************************************								*

	DOUBLE PRECISION FUNCTION EELZP(z,Icorr,tcorr)
	!$OMP DECLARE SIMD (EELZP)
	IMPLICIT NONE

	INTEGER nmaxi1,indi,inds
	DOUBLE PRECISION Rz, A, r, h, c , v, z, ic, t, zlam,					&
     Hcan, dt,dti,ind,dti1,  dt2,Icorr(*),tcorr(*)
	 
	
	COMMON /POINT/ r, h
	COMMON /TIME/ t
    COMMON /INTEG2/ nmaxi1
	COMMON /INTEG/ dti1
	COMMON /DAT/   v,zlam,Hcan,dt, dt2



	PARAMETER (c=2.99792458E8)


!*	EXTERNAL IPULSE

     	Rz = DSQRT(r*r + z*z)
	A = t - (Rz/c + ABS(z)/v)
    
	IF (A .LE. 0.) THEN
		ic = 0
	ELSE
!Cm		--------------------------------------------------------
!Cm		Ther is Used a Linear Interpolation for current Integral
!Cm		    _____________________________________________
!Cm		 	Note:
!Cm			For Icorr (Current Integral) It is necessary 
!Cm			to have a high accuracy 10-8 and timestep at 
!Cm			least 10 times lower than timestep simulation	
!Cm			_____________________________________________
!Cm		-----------------------------------------------------
		IF (A .LE. tcorr(nmaxi1)) THEN
			dti=tcorr(2)-tcorr(1)
			ind=((A)/dti)
			indi=floor(ind)
			inds=indi+1
			IF (indi .LE. 0) THEN
				ic=Icorr(1)*(ind-indi)
			ELSE
				ic=Icorr(indi)+(Icorr(inds)-Icorr(indi))*(ind-indi)
			ENDIF

		ELSE
			dti=tcorr(nmaxi1+1)-tcorr(nmaxi1)
			ind=nmaxi1+(A-tcorr(nmaxi1))/dti
			indi=floor(ind)
			inds=indi+1
			ic=Icorr(indi)+(Icorr(inds)-Icorr(indi))*(ind-indi)
		END IF	
	END IF
	IF (zlam.EQ.0) THEN
		EELZP = ic*(2*z*z - r*r)/(Rz**5)
	ELSE
		EELZP = ic*(2*z*z - r*r)*DEXP(-abs(z)/zlam)/(Rz**5)
	ENDIF
	END

!*																		*
!*																		*
!************************************************************************
!************************************************************************
!*	Induced and Radiated Magnetic Field Horizontal Component 						*
!*	*******************************************							*

	DOUBLE PRECISION FUNCTION BINDP(z)
	!$OMP DECLARE SIMD (BINDP)
	IMPLICIT NONE
	COMMON /POINT/ r, h
	COMMON /TIME/ t
	COMMON /DAT/   v,zlam,Hcan,dt, dt2

	DOUBLE PRECISION Rz, A, r, h, c , v, z, ic, t, zlam, &
     IPULSE, Hcan, dt,dic,BIND,BRAD,DIPULSE, dt2

		PARAMETER (c=2.99792458E8)

	EXTERNAL IPULSE,DIPULSE
     
     	Rz = DSQRT(r*r + z*z)
	A = t - (Rz/c + ABS(z)/v)

	IF (A .LE. 0) THEN
		ic = 0
	ELSE
  		ic=IPULSE(A)
  		dic=DIPULSE(A)
	END IF
	IF (zlam.EQ.0) THEN
		BIND = ic/(Rz**3)
		BRAD = dic/(c*Rz**2)
	ELSE
		BIND = ic*DEXP(-abs(z)/zlam)/(Rz**3)
		BRAD = dic*DEXP(-abs(z)/zlam)/(c*Rz**2)
	ENDIF
	BINDP=BRAD+BIND
	END
!*																		*
!*																		*
!************************************************************************
!************************************************************************
!*	Radiated and Induced Electric Field Horizontal Component 			*
!*	********************************************************			*

	DOUBLE PRECISION FUNCTION ERP(z)
	!$OMP DECLARE SIMD (ERP)
	IMPLICIT NONE
	COMMON /POINT/ r, h
	COMMON /TIME/ t
	COMMON /DAT/   v,zlam,Hcan,dt, dt2

	DOUBLE PRECISION Rz, A, r, h, c , v, z, ic, t, zlam,   &
     IPULSE, Hcan, dt,dic,DIPULSE,EINDR,ERADR, dt2

		PARAMETER (c=2.99792458E8)

	EXTERNAL IPULSE,DIPULSE
	Rz = DSQRT(r*r + (z - h)**2)
	A = t - (Rz/c + ABS(z)/v)

	IF (A .LE. 0) THEN
		ic = 0
		dic=0
	ELSE
   		ic = IPULSE(A)
	    dic=DIPULSE(A)
	END IF

	IF (zlam.EQ.0) THEN
		EINDR = ic*(h - z)/(Rz**4)
		ERADR = dic*(h - z)/(Rz**3)
	ELSE
		EINDR = ic*(h - z)*DEXP(-abs(z)/zlam)/(Rz**4)
		ERADR = dic*(h - z)*DEXP(-abs(z)/zlam)/(Rz**3)
	ENDIF

	ERP=(3*EINDR+ERADR/c)

	END
!*																		*
!*																		*
!************************************************************************

!************************************************************************
!*	Electrostatic Field Horizontal Component 							*
!*	************************************************					*

	DOUBLE PRECISION FUNCTION EELRP(z,Icorr,tcorr)
	!$OMP DECLARE SIMD (EELRP)
	IMPLICIT NONE
	COMMON /POINT/ r, h
	COMMON /TIME/ t
	COMMON /INTEG/ dti1
    COMMON /INTEG2/ nmaxi1
	COMMON /DAT/   v,zlam,Hcan,dt, dt2
	
	DOUBLE PRECISION Rz, A, r, h, c , v, z, ic, t, zlam,	&
      Hcan, dt,ind,dti1,dti, dt2, Icorr(*),tcorr(*)
	

	INTEGER indi,inds,nmaxi1

	PARAMETER (c=2.99792458E8)
     
	

    Rz = DSQRT(r*r + (z - h)**2)
	A = t - (Rz/c + ABS(z)/v)

	IF (A .LE. 0) THEN
		ic = 0
	ELSE

!Cm		--------------------------------------------------------
!Cm		Ther is Used a Linear Interpolation for current Integral
!Cm		    _____________________________________________
!Cm		 	Note:
!Cm			For Icorr (Current Integral) It is necessary 
!Cm			to have a high accuracy 10-8 and timestep at 
!Cm			least 10 times lower than timestep simulation	
!Cm			_____________________________________________
!Cm		-----------------------------------------------------
		IF (A .LE. tcorr(nmaxi1)) THEN
			dti=tcorr(2)-tcorr(1)
			ind=(A)/dti
			indi=floor(ind)
			inds=indi+1
			IF (indi .EQ. 0) THEN
			ic=Icorr(1)*(ind-indi)
			ELSE
			ic=Icorr(indi)+(Icorr(inds)-Icorr(indi))*(ind-indi)
			ENDIF
		ELSE
			dti=tcorr(nmaxi1+1)-tcorr(nmaxi1)
			ind=nmaxi1+(A-tcorr(nmaxi1))/dti
			indi=floor(ind)
			inds=indi+1
			ic=Icorr(indi)+(Icorr(inds)-Icorr(indi))*(ind-indi)
		END IF	
	END IF
	IF (zlam.EQ.0) THEN
		EELRP = ic*(h - z)/(Rz**5)
	ELSE
		EELRP = ic*(h - z)*DEXP(-abs(z)/zlam)/(Rz**5)
	ENDIF

	
	END
!*																		*
!*																		*
!************************************************************************

!************************************************************************
!*	Current Function with Heidler Expression							*
!*	************************************************					*

	DOUBLE PRECISION FUNCTION IPULSE(ta)
	
	IMPLICIT NONE

	COMMON /DAT2/ Ih1, Ih2, tao11, tao12,tao21,	&
     tao22, n1, n2, eta1, eta2
	COMMON /ITYPE/ triangular
	
	DOUBLE PRECISION i1, i2, Ih1, Ih2, tao11, tao12, &
     tao21, tao22, n1, n2, eta1, eta2,  ta
	LOGICAL(4) triangular


	IF (triangular .EQV. .FALSE.) THEN
		i1 = (Ih1/eta1)*(((ta/tao11)**n1)/(1+((ta/tao11)**n1)))*DEXP(-ta/tao21)
		IF (Ih2 .EQ.0 .OR. eta2 .EQ.0) THEN
			i2=0.D0
		ELSE
			i2 = (Ih2/eta2)*(((ta/tao12)**n2)/(1+((ta/tao12)**n2)))*DEXP(-ta/tao22)
		END IF
    
	 
		IPULSE = i1 + i2

	ELSE
		IF (ta .LE. tao11 .AND. ta .GT. 0.) THEN
			IPULSE=Ih1*ta/tao11
		ELSEIF(ta .GT. tao11 .OR. ta .LE. 2*tao21-tao11) THEN
			IPULSE=Ih1 ! -Ih1*(ta-tao11)/(2*tao21-tao11)
		ELSE
			IPULSE=0
		END IF
	ENDIF

	END FUNCTION IPULSE
!*																		*
!*																		*
!***************************************************************************

!***************************************************************************
!*	Current Derivate Function with Heidler Expression					*
!*	*************************************************					*

	DOUBLE PRECISION FUNCTION DIPULSE(ta)
	
	IMPLICIT NONE

	COMMON /DAT2/ Ih1, Ih2, tao11, tao12,tao21,tao22, n1, n2, eta1, eta2
	COMMON /ITYPE/ triangular

	DOUBLE PRECISION  Ih1, Ih2, tao11, tao12,	&
     tao21, tao22, n1, n2, eta1, eta2, ta, di11, di21,di12,di22,	&
		di1,di2,a,b
	LOGICAL(4) triangular
	
	IF (triangular .EQV. .FALSE.) THEN
		a=(ta/tao11)**n1
		b=(ta/tao12)**n2

		IF(Ih1 .EQ. 0 .OR. eta1 .EQ.0) THEN
			!di11 = 0
			!di21 = 0
			di1=0.D0
		ELSE
			!di11 = n1*DEXP(-ta/tao21)/tao11**n1/(1+(ta/tao11)**n1)**2*ta**(n1-1)
			!di21 = 1/tao21*DEXP(-ta/tao21)/tao11**n1/(1+(ta/tao11)**n1)*ta**n1
		
			di1=a*DEXP(-ta/tao21)/(1+a)*(n1/ta/(1+a)-1/tao21)
		      
	
		ENDIF
		IF(Ih2 .EQ. 0 .OR. eta2 .EQ.0) THEN
			!di12=0.D0
			!di22=0.D0
			di2=0.D0
		ELSE
			!di12 = n2*DEXP(-ta/tao22)/tao12**n2/(1+(ta/tao12)**n2)**2*ta**(n2-1)
			!di22 = 1/tao22*DEXP(-ta/tao22)/tao12**n2/(1+(ta/tao12)**n2)*ta**n2
			di2=b*DEXP(-ta/tao22)/(1+b)*(n2/ta/(1+b)-1/tao22)


		END IF
	
		!di1 = di11 - di21
		!di2 = di12 - di22

		IF(eta1 .EQ.0 .AND. eta2 .EQ.0) THEN
			DIPULSE=0.D0
		ELSEIF(eta1 .EQ.0) THEN
			DIPULSE = Ih2/eta2*di2
		ELSEIF(eta2 .EQ.0) THEN
			DIPULSE = Ih1/eta1*di1
		ELSE
			DIPULSE = Ih1/eta1*di1 + Ih2/eta2*di2
		ENDIF

	ELSE
		IF (ta .GT. 0. .AND. ta .LE. tao11 ) THEN
			DIPULSE=Ih1/tao11
		ELSEIF(ta .GT. tao11 .AND. ta .LE. 2*tao21-tao11) THEN
			DIPULSE=0 !-Ih1/(2*tao21-tao11)
		ELSE
			DIPULSE=0
		END IF
	ENDIF

	END FUNCTION DIPULSE
