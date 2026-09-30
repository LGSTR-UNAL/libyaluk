!*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!*                                                              %
!*                        YALUK-ATP CODE                        %
!*                                                              %
!*   - INIT CODE												%
!      Condiciones de Frontera se realiza con una linea tipo
!		sin p�rdidas con inducci�n electromagn�tica
!		se modifica la distancia de la l�nea para incluir un dx	
!		que es igual a v*dt a cada extremo de la l�nea
!		Grabar un solo archivo con todo los campos 
!		y mantenerlo en memoria sin leer archivo
!		Calcula el tiempo de inicio
!       Incluye perfil de la l�nea
!*   - Link with FOREIGN MODELS (Subroutine)                    %
!*    Incluye Resistencia en el cable                           %
!*    _ EV and EX equal for each conductor                                                          %
!*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!****************************************************************
!*                INPUT AND OUTPUT DATA DESCRIPTION             *
!*                                                              *
!*		Xdata(1):timestep   (Tiempo de Execuci�n)               *
!*		Xdata(2):timestop										*
!*		Xdata(3): n�mero de conductores                         *
!*		Xdata(4):Num_lin										*
!*		Xvar(1): kmax											*
!*		Xvar(2): SIZE_T											* 
!*		Xvar(3)=dx												*
!*		Xvar(4)=o												*
!*		Xvar(5)=e												*
!*		Xvar(6)=conductividad									*
!*		Xvar(7):t0												*
!*		Xvar(8)=dx2
!*      Xvar(9)=nlin											*
!*      Xvar(10):maximum number of division in x				*
!*      Xvar(11):maximum number of conductors                   *
!*      Xvar(12): status_perfil (1 .True. 0 .False.)
!*		Xvar(13): imprimir_inf                  				*
!*      Xvar(20)=ind_ini1 !initial variable for saving h values
!*		Xin(1):  t    (Tiempo de Execuci�n)                     *
!*		Xin(2):  Vant(1,1)                                      *
!*		Xin(3):  Vant(1,kmax)                                   *
!*		Xin(4):  Iant(1,1)                                      *
!*		Xin(5):  Iant(1,kmax)									* 
!*		Xin(6):  Vant(2,1)                                      *
!*		Xin(7):  Vant(2,kmax)                                   *
!*		Xin(8): Iant(2,1)                                       *
!*		Xin(9): Iant(2,kmax) 									* 
!*		Xin(10): Vant(3,1)                                      *
!*		Xin(11): Vant(3,kmax)                                   *
!*		Xin(12): Iant(3,1)                                      *
!*		Xin(13): Iant(3,kmax)									* 
!*		Xin(14): Vant(4,1)                                      *
!*		Xin(15): Vant(4,kmax)                                   *
!*		Xin(16): Iant(4,1)                                      *
!*		Xin(17): Iant(4,kmax)									* 
!*		Xout(1): Vi(1,1)                                        *
!*		Xout(2): Vi(1,kmax)                                     *
!*		Xout(3): Ii(1,1)                                        *
!*		Xout(4): Ii(1,kmax)                                     *
!*		Xout(5): -Ev(1,1)*hm                                    *
!*		Xout(6): -Ev(1,kmax)*hm                                 *
!*		Xout(7): Vi(2,1)                                        *
!*		Xout(8): Vi(2,kmax)                                     *
!*		Xout(9): Ii(2,1)                                        *
!*		Xout(10):Ii(2,kmax)                                     *
!*		Xout(11):Vi(3,1)                                        *
!*		Xout(12):Vi(3,kmax)                                     *
!*		Xout(13):Ii(3,1)                                        *
!*		Xout(14):Ii(3,kmax)                                     *
!*		Xout(15):V(4,1)                                         *
!*		Xout(16):Vi(4,kmax)                                     *
!*		Xout(17):Ii(4,1)                                        *
!*		Xout(18):Ii(4,kmax)                                     *
!****************************************************************

SUBROUTINE yaluk_init(xdata,xin,xout,xvar)
!DEC$ ATTRIBUTES STDCALL, DLLEXPORT::yaluk_init

!! For Compaq Fortran USE DFLIB  !Module which manage information on directories 
!USE DFLIB !visual fortran 6.0
!USE DFPORT
!USE IFPORT ! visual fortran 9.
!USE IFCORE
!USE IFQWIN
	
IMPLICIT NONE

INTEGER SIZE_T, i, j, k,cond, kmax, kmaxt, cond_max,SIZE_1,SIZE_2,g,n,ind_ini1,ind_ini,&
     nmaxi1,GNI,conductividad,num_lin,ALLOC_ERR,ERRNUM,ERRNUM2,ERRNUM_2, p_kmax
DOUBLE PRECISION xin(*), xout(*), xvar(*),xdata(*)
DOUBLE PRECISION r, rant,dt, tmax,XA,XB, 				&
     a,c, v, t, zlam, Hcan,tmax1,dt2,	dx2,dx,						&
     o, e,ang,t0,dti,ti,dti1,RSLT,h,angant,		&
     Ih1, Ih2, tao11, tao12,tao21,tao22, n1, n2, eta1, eta2,YR,L,	&
     BESSELI0,BESSELI1,t2,p,p12,pi, E_TIME,t02,				&
        s,rc,e0,mu,vl,kcampo,  	X0,Y0,XA0,XB0,YA0,YB0,t0_min,dist_camp
		
!	INTEGER(2) ihr, imin, isec, i100th,ihr2, imin2, isec2, i100th2	

DOUBLE PRECISION v_s,zlam_s,Hcan_s,o_s,e_s,dx_s,Ih1_s,Ih2_s, &
				tao11_s,tao21_s,n1_s,n2_s,X0_s,Y0_s,rc_s,YA0_s,YB0_s,XA0_s,XB0_s, &
				dt_s,tmax_s 
INTEGER Conductividad_s,nlin
INTEGER status_int 

CHARACTER (LEN=5) ncase_s
CHARACTER (LEN=3) indarchiv,kmax_s
CHARACTER (LEN=15) FORMATO
CHARACTER (LEN=30) archivo, casename,evalue,prueba_archiv
CHARACTER(512) dir
INTEGER(4) length,ncase
LOGICAL(4) status_perfil,status, status_bas,status_ef,status_hm,triangular,Imprimir_campo,campo_distante,Imprimir_inf,Read_campo
LOGICAL comparar
INTEGER, SAVE :: ALLOC_ERR_I

PARAMETER (pi=3.141592653589793D0)
PARAMETER (mu=(4*pi)*1E-7)
PARAMETER (c=2.99792458E8)
PARAMETER (e0=1/(mu*c**2))


DOUBLE PRECISION, ALLOCATABLE :: Rm(:,:),h_perfil(:,:),hm(:),Xi(:),hm_s(:),Xi_s(:),Cm(:,:),Lm(:,:),Li2(:,:),	&
		Ci(:,:),Li(:,:),LCi(:,:),LC(:,:),D1(:,:),D3i(:,:,:),D2(:,:),D6(:,:),D7(:,:),Mv1(:,:,:),Mv2(:,:,:),Mv3(:,:,:),		&
		D3(:,:),D5(:,:,:),Ro(:),Rf(:),Mi1(:,:,:),Mi2(:,:,:),Mi3(:,:,:),Mi4(:,:,:),				&
		Evini(:),Evfin(:),Ev(:), EZ(:),EH(:),Ex(:,:,:),					&
		Zc(:,:),Zc_o(:,:,:),Zci(:,:),Zci_o(:,:,:)

DOUBLE PRECISION, ALLOCATABLE :: Icorr(:),tcorr(:),i02(:),i12(:)

! Evini: 3 dimensiones (primera n�mero de conductores, segunda div, en x, y tercera N�mero de l�neas
! Evfin: 3 dimensiones (primera n�mero de conductores, segunda div, en x, y tercera N�mero de l�neas
! Ex: 4 dimensiones (primera n�mero de conductores, segunda div, en t, terc. div en x y cuarta N�mero de l�neas
! kmaxt; n�mero m�ximo de divisiones en las l�neas
! max_lin: numero m�xmio de l�neas
!		Ci(4,4),Li(4,4),Lm(4,4),Cm(4,4),Ro(4),Rf(4),Rt,					&
!    	LCi(4,4),LC(4,4),D1(4,4),D3i(4,4),D2(4,4),D4i(4,4),				&
!       D3(4,4),
     	
!*	----------------------------------------------
!Cm	COMMON DATA FOR LINE DESCRIPTION
!Cm		r: Distance for Electromagnetic Field Calculation
!Cm		hm: Line Height
!Cm		XR: Cordinate X
!Cm		YR: Cordinate Y
!*	----------------------------------------------
!*	----------------------------------------------
!Cm	COMMON DATA /DAT/ FOR CHANNEL
!Cm		c: light velocity
!Cm		v: return stroke velocity
!Cm		zlam: Atenuation MTL model if zlam=0  TL model
!Cm		Hcan: Height Channel
!Cm		dt: Discretization Time
!Cm		t0: Initial Calculation Time
!*	----------------------------------------------
!*	----------------------------------------------
!Cm	COMMON /DAT2/ DATA FOR HEIDLER CURRENT WAVEFORM
!*	-----------------------------------------------
!Cm		Ih1: Current Amplitude(function 1)
!Cm		Ih2: Curren Amplitude (function 2)
!Cm		tao11: Front Time Constant (function 1)
!Cm		tao21: Decay Time Constant (function 1)
!Cm		tao12: Front Time Constant (function 2)
!Cm		tao22: Decay Time Constant (function 2)
!Cm		n1: Amplitude Constant (function 1)
!Cm		n2: Amplitude Constant (function 2)
!*	-----------------------------------------------
!*	----------------------------------------------
!Cm	COMMON /MAGN/ DATA FOR MATRIX DIMMENSIONS
!*	-----------------------------------------------
!Cm		cond: Number of Line Conductors
!Cm		kmax: Maximum number of divisions on the Line
!Cm		nmax: Maximum number of divisions on Time
!*	------------------------------------------------

COMMON /COND2/ conductividad
COMMON /COND/  o,e
COMMON /LINE_YALUK/  SIZE_1,SIZE_2
COMMON /POINT/ r,h
COMMON /DAT/   v,zlam,Hcan,dt, dt2
COMMON /DAT2/ Ih1, Ih2, tao11, tao12,tao21,tao22, n1, n2, eta1, eta2
COMMON /ITYPE/ triangular
 
COMMON /INTEG2/ nmaxi1
COMMON /INTEG/ dti1
COMMON /LPARAM/ L
!	SAVE Cm,Lm,LC,LCi,Li,Ci,D1,D2,D3i,D4i,/LPARAM/ 
!	SAVE Ex,Evini,Evfin,Icorr, tcorr,i02,i12

EXTERNAL IPULSE,BESSELI0,BESSELI1,comparar


!kmaxt=100	! kmaxt; n�mero m�ximo de divisiones en las l�neas
!max_lin=20	! max_lin: numero m�xmio de l�neas
	  

archivo='yaluk.ini'
ERRNUM=0
Imprimir_campo=.FALSE.
!****************************
!* STATUS_FILE.YLK
!* FILE READ for identify simulation parameters between lines
!**********************************
OPEN (UNIT = 11, FILE = 'status_file.ylk', FORM='UNFORMATTED', STATUS = 'OLD', ERR=1011,IOSTAT=ERRNUM2)
READ (UNIT=11) dt_s,tmax_s,ncase_s,casename,t0_min,kmaxt,cond_max,imprimir_inf 
CLOSE (UNIT=11)
    IF (Xvar(12) .EQ. 1) THEN
        status_perfil=.TRUE.
    ELSE
        status_perfil=.FALSE.
    ENDIF



1011 SELECT CASE(ERRNUM2>0)


	CASE (.TRUE.)
	!write(*,*)'ERRNUM=', ERRNUM2
	!  Get current directory
	!dir = FILE$CURDRIVE
	!length = GETDRIVEDIRQQ(dir)
	!IF (length .EQ. 0) THEN
	!  WRITE (*,*) 'Failed to get current directory'
	!END IF
	
	CALL GETCWD(dir, status_int)	
	IF (status_int .NE. 0) THEN	
	  WRITE (*,*) 'Failed to get current directory. Error:  ',status_int	
	END IF	
	
	
	
	
1012 SELECT CASE(ERRNUM>0)
	CASE (.TRUE.)
			WRITE(*,*) 'Configuration file not founded:',TRIM(dir),'/',archivo
			write(*,*) 'please enter name file:'
			read(*,*) archivo
	END SELECT 

		CALL OUTSIX( "*****************************************************", 80)
		write(*,*) "! NATIONAL UNIVERSITY OF COLOMBIA"
		write(*,*) "! Runing YALUK_Linux version with R in conductor 2020 - Line Profile"
		CALL GETENV ("CPU", evalue) 
		write(*,*) '! Using CPU=',evalue
		write(*,*) '! Authors: Ernesto Perez - Edison Soto'
		write(*,*) '! Contact: eperezg@unal.edu.co - easotor@unal.edu.co - jpnorenam@unal.edu.co'
		write(*,*) "******************************************************"


	OPEN (UNIT = 12, FILE = TRIM(archivo), STATUS = 'OLD', ERR=1012,IOSTAT=ERRNUM)
		READ(12,*) casename !Case Name
		READ(12,*) ncase    !Case Number
		READ(12,*) nlin ! Maximum number of Lines
		READ(12,*) cond_max     ! Maximum number of conductors
		IF (cond_max>8) THEN
	    		    cond_max=8
	    ELSEIF(cond_max .EQ. 0)THEN
	    	cond_max=4
	    	write(*,*) '*WARN* maximum number of conductors was set equal to zero and it is not logical'
   	    	write(*,*) '*WARN* Correct the problem modifying file Yaluk.ini'
            write(*,*) '*INF* maximum number of conductors will be set automatically in 4'
	    ENDIF

	CLOSE (UNIT=12)
    !****************************
    ! CONFIGURATION MISC. FILE
    !******************************
    WRITE(*,*) 'directorio actual: ', casename, ncase, dir
	! WRITE(*,*) 'please enter case files folder name:'
	! READ(*,*) casename !Case Name

    OPEN (UNIT = 21, FILE = 'yaluk_status.ini', FORM='FORMATTED', STATUS = 'OLD', ERR=9011,IOSTAT=ERRNUM2)
        READ   (21, *) Imprimir_campo  !Identify if Electromagnetic field will be printed
        READ   (21, *) Imprimir_inf    !!Identify if information of the process will be printed
		status_perfil=.FALSE.
        READ   (21, *, end=995) Read_campo      !Identify if Electromagnetic field will extarnally read
	    
		READ   (21, *) status_perfil
995 CONTINUE
		
    CLOSE (21)
    9011 SELECT CASE(ERRNUM2>0)
	        CASE (.TRUE.)
            Imprimir_campo=.FALSE.
            Imprimir_inf=.FALSE.
            Read_campo=.FALSE.
            status_perfil=.FALSE.
          END SELECT 


	!status = CHANGEDIRQQ(dir(1:length)//'\'//casename(1:LEN_TRIM(casename)))
	write(*,*) '*INF* Trying to open: ',TRIM(dir)//'/'//casename(1:LEN_TRIM(casename))	
	CALL CHDIR(TRIM(dir)//'/'//casename(1:LEN_TRIM(casename)), status_int)
		IF (status_int .EQ. 0) THEN
				!dir = FILE$CURDRIVE
				!length = GETDRIVEDIRQQ(dir)
			CALL GETCWD(dir, status_int)
			write(*,*) '*INF* Current Directory for cases:  ',TRIM(dir)
		ELSE
			!status = MAKEDIRQQ(dir(1:length)//'\'//casename(1:LEN_TRIM(casename)))
			write(*,*) 'Trying to create: ',TRIM(dir)//'/'//casename(1:LEN_TRIM(casename))	
			CALL SYSTEM( 'mkdir ' //TRIM(dir)//'/'//casename(1:LEN_TRIM(casename)), status_int)			
			
			IF (status_int .EQ. 0) THEN
				write(*,*) ' Directory for cases was created ',TRIM(dir)//'/'//casename(1:LEN_TRIM(casename))
				write(*,*) '**************************************************************'
				write(*,*) '*  PLEASE PUT THE LINE AND MISCELANEO FILES IN THIS DIRECTORY*'
				write(*,*) '**************************************************************'
				!pause

			ELSE
				write(*,*) ' Error Creating Directory: ',TRIM(dir)//'/'//casename(1:LEN_TRIM(casename))
				write(*,*) ' Please check permissions'
			ENDIF
		END IF
		!****************************************************************************
		!OPENING LINE DATA FOR CALCULATING t0
		! CALCULATING MAXIMUM KMAXT
		!****************************************************************************
		cond=1 ! int(Xdata(3))
		ALLOCATE (Rm(cond,cond),hm(1:cond),Xi(1:cond), STAT=ALLOC_ERR)
		t0_min=1.D200
		kmaxt=0 !n�mero m�ximo de divisiones
		IF (nlin.GT.0)		THEN
		121		FORMAT(I5.5)
		WRITE(ncase_s,FMT=121) ncase
			CALL OPEN_CORR_FILE (Ih1,Ih2,tao11,tao21,tao12,tao22,n1,n2,X0,Y0,ncase_s)
			CALL OPEN_MISC_FILE (v,zlam,Hcan,tmax1,SIZE_2,conductividad,o,e,dx,triangular,dist_camp)
			DO i=1,nlin
				120 FORMAT(I3.3) 
				WRITE(indarchiv,FMT=120) int(i)
				CALL OPEN_LINE_FILE(rc,XA0,YA0,XB0,YB0,hm,Rm,Xi,cond,indarchiv)
				
				IF (r .NE.-1) THEN
					a=ATAN2((YA0-YB0),(XA0-XB0))
					XA=(XA0-X0)*cos(a)+(YA0-Y0)*sin(a)
					XB=(XB0-X0)*cos(a)+(YB0-Y0)*sin(a)
					YR=(YA0-Y0)*cos(a)-(XA0-X0)*sin(a)
					L=XA-XB
                    
                   if (nint(L/dx)+4 .GT. kmaxt) THEN
                     kmaxt=nint(L/dx)+4
                     
                    END IF
					IF (XA .LT. 0 .AND. XA .GT. XB) THEN
						t0=(XA**2+Yr**2)**0.5/c;
					ELSEIF (XB .LT. 0 .AND. XB .GT. XA) THEN
						t0=(XB**2+Yr**2)**0.5/c;
					ELSEIF (XA .GT. 0 .AND. XA .LT. XB) THEN
						t0=(XA**2+Yr**2)**0.5/c;
					ELSEIF (XB .GT. 0 .AND. XB .LT. XA) THEN
						t0=(XB**2+Yr**2)**0.5/c;
					ELSE
						t0=abs(Yr)/c;
					ENDIF
					
					IF (t0.LE.t0_min) THEN
						t0_min=t0
					END IF
				END IF

				

			END DO
			t0_min=t0_min-dt
		ELSE
			t0_min=0
			kmaxt=200
		END IF
		
    !****************************************************************************
    !****************************************************************************
    IF (Imprimir_inf .EQV. .TRUE.) THEN
			write(*,*)'*INF* Initial time t0_min=',t0_min !WRITES THE INITIAL TIME
	ENDIF	
		DEALLOCATE (Rm,hm,Xi,STAT=ALLOC_ERR)
		!____________________________________________________________________________
		!****************************************************************************

		dt_s=Xdata(1)			!data for comparing if the file was opened before
		tmax_s=Xdata(2)
		101		FORMAT(I5.5)
		!File for knowing if any line is already calculated
		WRITE(ncase_s,FMT=101) ncase
		OPEN (UNIT = 13, FILE = 'status_file.ylk', FORM='UNFORMATTED', STATUS = 'UNKNOWN', ERR=1013,IOSTAT=ERRNUM)
		WRITE (UNIT=13) dt_s,tmax_s,ncase_s,casename,t0_min,kmaxt,cond_max,imprimir_inf
		CLOSE(UNIT=13)
		write(*,*)'CASE NO:',ncase_s
1013	SELECT CASE(ERRNUM>0)
			CASE (.TRUE.)
			WRITE(*,*) '*ERROR* saving status_file.ylk', ERRNUM
			END SELECT 

		END SELECT
ERRNUM=0
	IF (Imprimir_inf) THEN
	write(*,*) '*INF* Initializing Lines',indarchiv
      ENDIF




	
!*		-----------------------
!*		Setting Input Variables
!*		-----------------------
!		tmax = 20.D-6	- maximum simulation time given by ATP
!		dt=1D-8			- timestep given by ATP
!*		--------------------------------------------------------


dt=Xdata(1)
tmax=Xdata(2)
cond=int(Xdata(3))
num_lin=Xdata(4)
!t0_min=Xvar(9)
100 FORMAT(I3.3) 
WRITE(indarchiv,FMT=100) int(Xdata(4))
IF (imprimir_campo .EQV. .TRUE.) THEN
    open(UNIT = 99, FILE = 'cond_'//indarchiv//'.dat', FORM='UNFORMATTED', STATUS = 'UNKNOWN')
    write(UNIT = 99) cond
    close(UNIT = 99)
ENDIF



	  ALLOCATE (Rm(cond,cond),hm(1:cond),Xi(1:cond),hm_s(1:cond),Xi_s(1:cond),Cm(cond,cond),Lm(cond,cond),Li2(cond,cond),					&
		Ci(cond,cond),Li(cond,cond),LCi(cond,cond),LC(cond,cond),		&
		D1(cond,cond),D3i(cond,cond,2),D2(cond,cond),D5(cond,cond,2),D6(cond,cond),D7(cond,cond),Zc(cond,cond),Zc_o(cond,cond,2),Zci_o(cond,cond,2),		&
		Zci(cond,cond),D3(cond,cond),Ro(cond),Rf(cond), STAT=ALLOC_ERR)


IF (dt .EQ. dt_s .AND. tmax .EQ. tmax_s) THEN
	status_bas=.TRUE.
	open(UNIT = 14, FILE = 'bas_'//ncase_s//'_'//casename(1:LEN_TRIM(casename))//'_'//indarchiv//'.dat', FORM='UNFORMATTED', STATUS = 'OLD',ERR=1014,IOSTAT=ERRNUM)
	read(UNIT = 14, ERR=1021,IOSTAT=ERRNUM_2) v_s,zlam_s,Hcan_s,Conductividad_s,o_s,e_s,dx_s,Ih1_s,Ih2_s,tao11_s,tao21_s,n1_s,n2_s,X0_s,Y0_s,rc_s,YA0_s,YB0_s,XA0_s,XB0_s,hm_s,Xi_s,kmax,t0
	IF (Imprimir_inf) THEN
	write(*,*) '*INF* Loading database of previous case: ',indarchiv
      ENDIF
	close (UNIT = 14)
	
1014 SELECT CASE(ERRNUM>0)
		CASE (.TRUE.)
			status_bas=.FALSE.
		END SELECT 
		
1021 SELECT CASE(ERRNUM_2>0)
		CASE (.TRUE.)
			status_bas=.FALSE.
		END SELECT 

ELSE
	status_bas=.FALSE.
ENDIF




!indarchiv:  Es el indicador del n�mero de l�nea el cual se utiliza para los datos de l�nea
!			y para los archivos internos que se usan en Yaluk_exec
t=Xin(1)
!*		-----------------------------------------------------------------
!*         READING MISCELANEOUS PARAMETERS
!*		_________________________________________________________________
IF (Imprimir_inf ) THEN
	    write(*,*) '*Opening Miscelaneous file' 
	    
    ENDIF
CALL OPEN_MISC_FILE (v,zlam,Hcan,tmax1,SIZE_2,conductividad,o,e,dx,triangular,dist_camp)
IF (Imprimir_inf ) THEN
	    write(*,*) '*File succesfully read' 
	    
    ENDIF
!		------------------------------------
!        Reading Lightning Current Parameters
!		------------------------------------
IF (Imprimir_inf ) THEN
	    write(*,*) '*Opening Current file' 
	    
    ENDIF
CALL OPEN_CORR_FILE (Ih1,Ih2,tao11,tao21,tao12,tao22,n1,n2,X0,Y0,ncase_s)
IF (Imprimir_inf ) THEN
	    write(*,*) '*File succesfully read' 
	    
    ENDIF
!	------------------------------------------------------
!        Reading Line Position and Stroke Parameters
!	------------------------------------------------------
!	---------------------------
!*	Line Coordinates
!*	   <---- Xr-----> ==================
!*	   !             +      Line
!*	   !           +
!*	   !         +   
!*	YR !       +   
!*	   !     +  R
!*	   !   +
!*	   ! +
!*	   * (Lightning strike)   
!*	---------------------------
!      a:=atan2(YA0-YB0,XA0-XB0)
!      XA:=(XA0-X0)*cos(a)+(YA0-Y0)*sin(a)
!      XB:=(XB0-X0)*cos(a)+(YB0-Y0)*sin(a)
!      YR :=(YA0-Y0)*cos(a)-(XA0-X0)*sin(a)
!	L=XA-XB
!inicializando valores
IF (Imprimir_inf ) THEN
	    write(*,*) '*Opening Line file' 
	    
    ENDIF
 	CALL OPEN_LINE_FILE(rc,XA0,YA0,XB0,YB0,hm,Rm,Xi,cond,indarchiv)
IF (Imprimir_inf ) THEN
	    write(*,*) '*File succesfully read' 
	    
    ENDIF
    

!configuring maximum first window from front time
IF (.NOT. read_campo) THEN
tmax1=1.8*(tao11**n1*tao21*n1)**(1/(1+n1))
ELSE
tmax1=tmax-dt;
ENDIF

SIZE_T = nint(tmax / dt)		!number of time division for the whole window 
SIZE_1= nint(tmax1/dt)			!number of time division for first window
dti=.05D-8
dti1=dti
nmaxi1=nint(dt*SIZE_1/dti)
dt2= (tmax-tmax1)/SIZE_2 !time step for second window
SIZE_2=SIZE_2+1			!number of time division for first window

a=ATAN2((YA0-YB0),(XA0-XB0))
XA=(XA0-X0)*cos(a)+(YA0-Y0)*sin(a)
XB=(XB0-X0)*cos(a)+(YB0-Y0)*sin(a)
YR=(YA0-Y0)*cos(a)-(XA0-X0)*sin(a)
L=XA-XB

!*************************************************************************
!*      Courant Criteria
!*		-------------------------------------------------------
!*		dx=10			- line division size (given by equation)
!*		--------------------------------------------------------
!*************************************************************************

IF (dx.LE.dt*(1.1*c) .OR. dx .EQ.0) THEN
	dx=dt*(2*c)         !

    IF (Imprimir_inf ) THEN
	    write(*,*) '*WARNING*  dx is too small or equal to 0. unconvergence is very likely' 
	    write(*,*) '*WARNING* (500)dx is automatically computed dx=',dx 
    ENDIF
ENDIF
		

dx2=c*dt

kmax = nint((L-2*dx2)/dx+2) !taking into account the two bergeron lines for coupling with ATP
if (kmax.LE.3) THEN
    		write(*,*) '*INF* Number of divisions on the line is too small - reduce dt value'
    		STOP 
ENDIF
dx=(L-2*dx2)/real(kmax-2)
IF ((L-(kmax-2)*dx-2*dx2).GE.(0.5*dx)) THEN
	kmax=kmax+1
END IF

IF (Imprimir_inf .EQV. .TRUE.) THEN
    write(*,*)'*INF* (516)longitude=',(kmax-2)*dx+2*dx2
    write(*,*)'*INF*     ,dx2=',dx2, '  dx=',dx
END IF


IF (status_perfil .EQV. .TRUE.)THEN
    p_kmax=kmax !Define the maximum number of points in x for line when is introduced height non_linearities
    IF (Imprimir_inf .EQV. .TRUE.) THEN
        write(*,*)'*INF* Using line profile kmax=',kmax
    END IF
ELSE
    p_kmax=1
    IF (Imprimir_inf .EQV. .TRUE.) THEN
        write(*,*)'*INF* Homogenuous line kmax=1'
    END IF

ENDIF

!**********************
!* ALOCATING Mvs and Mis and hperfil variables
!*********************+

ALLOCATE (Mv1(cond,cond,p_kmax),Mv2(cond,cond,p_kmax),Mv3(cond,cond,p_kmax),Mi1(cond,cond,p_kmax),Mi2(cond,cond,p_kmax),&
        Mi3(cond,cond,p_kmax),Mi4(cond,cond,p_kmax),h_perfil(1:p_kmax,1:cond),STAT=ALLOC_ERR)
 
   
    IF (ALLOC_ERR .NE. 0) THEN
        IF (Imprimir_inf .EQV. .TRUE.) THEN
	        write(*,*) '*ERROR* Allocating Mvs (536_Yaluk.init) Err=',ALLOC_ERR 
        ENDIF
    ENDIF
    
    IF (status_perfil .EQV. .TRUE.) THEN
        CALL OPEN_LINE_PROFILE(h_perfil,p_kmax,cond,indarchiv)
    ELSE
        h_perfil(1,:)=hm(:)
    IF (Imprimir_inf .EQV. .TRUE.) THEN
        write(*,*)'*INF* Height ',h_perfil(1,:)
    END IF

    ENDIF

!**********************

IF( .NOT. ALLOCATED(Evini)) THEN
    ALLOCATE (Evini(SIZE_T+2),Evfin(SIZE_T+2),Ev(SIZE_T+2),EZ(SIZE_T+2),EH(SIZE_T+2),			&
		  Ex(cond,kmax+1,SIZE_T+2),STAT=ALLOC_ERR)
    IF (ALLOC_ERR .NE. 0) THEN
        IF (Imprimir_inf .EQV. .TRUE.) THEN
	        write(*,*) '*ERROR* Allocating Evini (542)Yaluk.init) Err=',ALLOC_ERR 
        ENDIF
    ENDIF

ENDIF

IF( .NOT. ALLOCATED(Icorr)) THEN
    ALLOCATE (Icorr(nmaxi1+11*SIZE_2),tcorr(nmaxi1+11*SIZE_2),i02(SIZE_T+2),i12(SIZE_T+2),STAT=ALLOC_ERR_I) 
    IF (ALLOC_ERR .NE. 0) THEN
        IF (Imprimir_inf .EQV. .TRUE.) THEN
	        write(*,*) '*ERROR* allocating Icorr (552)Yaluk.init) Icorr=',ALLOC_ERR_I 
        ENDIF
    ENDIF
ENDIF
!-------------------------------------------------------
!ALLOC_ERR_I - INDICA SI LA CORRIENTE YA FUE CALCULADA
!-------------------------------------------------------


IF (status_bas .AND. comparar(hm,hm_s,cond)) THEN
	status_bas=.TRUE.
	status_hm=.TRUE.
!	write(*,*) 'comparando h(verdadero)'
ELSE
	status_bas=.TRUE.
	status_hm=.FALSE.
!	write(*,*) 'comparando h(falso)'

ENDIF

IF (status_bas .AND. v_s .EQ. v .AND. zlam_s .EQ. zlam .AND. Hcan .EQ. Hcan .AND. Conductividad_s .EQ. conductividad		&
	.AND. o_s .EQ. o .AND. e_s .EQ. e .AND. dx_s .EQ. dx .AND. Ih1_s .EQ. Ih1 .AND. Ih2_s .EQ. Ih2 .AND. tao11_s .EQ. tao11	&
	.AND. tao21_s .EQ. tao21 .AND. n1_s .EQ. n1 .AND. n2_s .EQ. n2 .AND. X0_s .EQ. X0 .AND. Y0_s .EQ. Y0.AND. rc_s .EQ. rc	&
	.AND. YA0_s .EQ. YA0 .AND. YB0_s .EQ. YB0 .AND. XA0_s .EQ. XA0.AND. XB0_s .EQ. XB0 ) THEN

	    open(UNIT = 15, FILE = ncase_s//'_'//casename(1:LEN_TRIM(casename))//'_'//indarchiv//'.dat', FORM='UNFORMATTED', STATUS = 'OLD',ERR=1015,IOSTAT=ERRNUM)
		read(UNIT = 15) Evini,Evfin,Ex,D3i,D5,Zc_o,Zci_o,Mv1,Mv2,Mv3,Mi1,Mi2,Mi3,Mi4

IF (Imprimir_inf .EQV. .TRUE.) THEN
	write(*,*) '*INF* Loading EM data from calculated file: line'//indarchiv
ENDIF

	status_ef=.FALSE.
	close (UNIT = 15)
1015 SELECT CASE(ERRNUM>0)
		CASE (.TRUE.)
			status_ef=.TRUE.
		END SELECT 
ELSE
	status_ef=.TRUE.
ENDIF

 IF (.NOT. status_bas .OR. status_hm  .OR. status_ef ) THEN
!******************************************************************************************
!				CALCULATING LINE PARAMETERS
!******************************************************************************************



!*--------------------------------------------------------------
!*
!       MATRICES CONSTANTES PARA EL CALCULO DE DIFERENCIAS 
!                        FINITAS
!* D1:  Constante para calculo de frontera inicio l�nea (dx/dt*C*Rom-unos)
!* D3i: Constante inicio l�nea  inv(dx/dt*C*Rom+unos)
!* D2:  Constante para calculo de frontera final l�nea (dx/dt*C*Rfm-unos)
!* D4i  Constante para calculo de frontera final l�nea  inv(dx/dt*C*Rfm+unos)
!* D5   Constante para calculo de frontera final l�nea  (dx/dt*C*Zc)
!* Zci: Inversa Impedancia Caracter�stica
!*
!*--------------------------------------------------------------
!	READ(12) D3i_o,D5_o,Zc_o,Zci,Mv1_o,Mv2_o,Mv3_o,Mi1_o,Mi2_o,Mi3_o,Mi4_o



    DO k=1,p_kmax
        hm(:)=h_perfil(k,:)
	    CALL LINEPARAM(cond,Xi,hm,rc,Lm,Cm,Zc)
		CALL inv(Zc,cond,cond,Zci)
		CALL mult(Lm,Cm,LC,cond,cond,cond)
		CALL inv(Cm,cond,cond,Ci)
		CALL inv(Lm,cond,cond,Li)
		CALL inv(LC,cond,cond,LCi)
!*	--------------------------------
!*	matrix Construction M
!*	--------------------------------
    Mv1(:,:,k)=.5*dt/dx*Ci
    Mv2(:,:,k)=.25*dt*dt/dx*LCi
    Mv3(:,:,k)=.25*dt*dt/dx*matmul(Rm,LCi)
    Mi1(:,:,k)=dt*Li
    Mi2(:,:,k)=.5*dt*dt/(dx*dx)*LCi
    Mi3(:,:,k)=.5*dt*dt*matmul(Rm,Li)
!********************************************
!Considering loss in the conductor
    D6=dt*Rm
    D6=matmul(D6,Li)
    D7=0.5*dt*dt*matmul(Rm,Rm)
 !Li2=matmul(Li,Li)
    D7=matmul(D7,Li)
    Mi4(:,:,k)=-D6-D7
!********************************************	

!*	--------------------------------
!*	Construcci�n de D1, D2, D3i, D4i
!*	--------------------------------
	D2=D1
	!D1=D5
!	D2=D5
!	D4=D5
	



    DO 21, i=1,cond
!        D1(i,i)=D1(i,i)-1
!        D2(i,i)=D2(i,i)-1
        D3(i,i)=D3(i,i)+1
!        D4(i,i)=D4(i,i)+1
        Mi4(i,i,k)=Mi4(i,i,k)+1 !Considering loss in the conductor
21	CONTINUE
!	CALL inv(D4,cond,cond,D4i)
    
        if (k.EQ.1) THEN
            Zc_o(:,:,1)=Zc
            Zci_o(:,:,1)=Zci
        	D1=dx/dt*Cm
        	D5(:,:,1)=matmul(D1,Zc)
        	D3=D5(:,:,1)
      	    DO i=1,cond
                D3(i,i)=D3(i,i)+1
            ENDDO
		IF (Imprimir_inf .EQV. .TRUE.) THEN
        	write(*,*) '*INF* D3=',D3
        END IF

            CALL inv(D3,cond,cond,D3i(:,:,1))
            IF(p_kmax.EQ.1)THEN
                Zc_o(:,:,2)=Zc
                Zci_o(:,:,2)=Zci
                D5(:,:,2)=D5(:,:,1)
                D3i(:,:,2)=D3i(:,:,1)
            ENDIF

        ELSEIF (k.EQ.kmax)THEN
            Zc_o(:,:,2)=Zc
            Zci_o(:,:,2)=Zci
        	D1=dx/dt*Cm
        	D5(:,:,2)=matmul(D1,Zc)
        	D3=D5(:,:,2)
      	    DO i=1,cond
                D3(i,i)=D3(i,i)+1
            ENDDO
            CALL inv(D3,cond,cond,D3i(:,:,2))

        ENDIF

    ENDDO
    
!********************FINISH CALCULATING M************************
ENDIF

IF (Imprimir_inf) THEN
	write(*,*) '*INF* Surge Impedance Zc=',Zc
	write(*,*) '*INF* Surge Impedance Zci=',Zci
	write(*,*) '*INF* Inductance Lm=',Lm
	write(*,*) '*INF* Capacitance Cm=',Cm
	write(*,*) '*INF* Capacitance D2=',D2
	write(*,*) '*INF* Capacitance Mi4=',Mi4
	write(*,*) '*INF* Capacitance hm=',hm
	write(*,*) '*INF* Capacitance rc=',rc
ENDIF



IF (status_ef) THEN


	
!*************************************************
! Calculating constants for heidler function
!----------------------------------------------------
	IF (tao11 .EQ.0 .OR. tao21 .EQ.0 ) THEN
		eta1=0.D0
	ELSE
		eta1 = DEXP(-(tao11/tao21)*(n1*tao21/tao11)**(1/n1)) !calculating eta1 for heidler current
	END IF
	IF (tao12 .EQ.0 .OR. tao22 .EQ. 0) THEN
		eta2 =0.D0
	ELSE
		eta2 = DEXP(-(tao12/tao22)*(n2*tao22/tao12)**(1/n2)) !calculating eta2 for heidler current
	END IF


!*	--------------------------------------------------------------
!                  Computing Initial time
!*	--------------------------------------------------------------
	IF (XA .LT. 0 .AND. XA .GT. XB) THEN
		t0=(XA**2+Yr**2)**0.5/c;
	ELSEIF (XB .LT. 0 .AND. XB .GT. XA) THEN
		t0=(XB**2+Yr**2)**0.5/c;
	ELSEIF (XA .GT. 0 .AND. XA .LT. XB) THEN
		t0=(XA**2+Yr**2)**0.5/c;
	ELSEIF (XB .GT. 0 .AND. XB .LT. XA) THEN
		t0=(XB**2+Yr**2)**0.5/c;
	ELSE
	    t0=abs(Yr)/c;
	ENDIF






!***********************************************************************
! EVALUATING CALCULATION OF ELECTROMAGNETIC FIELD  
!
! Here is evaluated if the EM Field will be calculated depending on
!
! 1. lightning is too far from the line (defined on variable dist_camp)
! 2. Current is zero (Ih1=0) used when just the line model will be used
! 3. 
!***********************************************************************

	IF (Ih1 .EQ. 0) THEN
	        IF (Imprimir_inf) THEN
		        write(*,*) '*WARN* Current Amplitude is equal to zero, then, EM field is zero'
		        write(*,*) '*INF* Routine will be use only the line model'
   		        write(*,*) '*INF* too far away - for not disregarding EM field, change miscelaneo.txt file'
            ENDIF

	 campo_distante=.TRUE.
	ELSE    
	IF (t0 .GT. dist_camp/c) THEN
	
	        IF (Imprimir_inf .EQV. .TRUE.) THEN
		        write(*,*) '*INF* Disregarding the EM field calculation for line_',indarchiv
		        write(*,*) '*INF* line too far away - dist_camp=',dist_camp
   		        write(*,*) '*INF* too far away - for not disregarding EM field, change miscelaneo.txt file'
            ENDIF
		    campo_distante=.TRUE.
	    ELSE
		    campo_distante=.FALSE.
	    ENDIF
	ENDIF
!***************************************************************
!***************************************************************

!*	******************************************************************
!*                  Current Integral Computation
!*	--------------------------------------------------------------

!	Calculating First Window
!	--------------------------
	ti=0.D0

IF (.NOT. campo_distante) THEN
    IF (triangular .EQV. .TRUE.) THEN
	    IF (Imprimir_inf .EQV. .TRUE.) THEN
	        write(*,*) '*INF* Using triangular waveform'
        ENDIF
    ELSE
	    IF (Imprimir_inf .EQV. .TRUE.) THEN
	        write(*,*) '*INF* Using HEIDLER waveform'
        ENDIF
    ENDIF

	IF (ALLOC_ERR_I .EQ. 0 ) THEN 
		IF (Imprimir_inf .EQV. .TRUE.) THEN
		    write(*,*) '*INF* Calculating Icorr'
        ENDIF 
      DO 1, i=1,nmaxi1
      
		IF (ti .EQ. 0) THEN
			Icorr(i)=0.D0
			tcorr(i)=0.D0
		ELSE
			IF (triangular .EQV. .FALSE.) THEN

			CALL INTEGRAL(IPULSE,0.D0,ti,RSLT,1.D-7,GNI)
			Icorr(i)=RSLT
			tcorr(i)=ti
		
			
			ELSE
			!*****************************
			! Integral for double ramp function
			!*******************************

				IF (ti .LE. tao11 ) THEN
					Icorr(i)=.5*Ih1*ti**2/tao11
					tcorr(i)=ti
				ELSEIF(ti .GT. tao11 .OR. ti .LT. 2*tao21-tao11) THEN
					Icorr(i)=.5*Ih1*tao11 !+ Ih1*(ti-tao11)-Ih1*(.5*ti*ti-ti*tao11)/(2*tao21-tao11)
					tcorr(i)=ti
				ELSE
					Icorr(i)=.5*Ih1*tao11 !+ Ih1*(tao21-tao11)-Ih1*(.5*tao21*tao21-tao21*tao11)/(2*tao21-tao11)
					tcorr(i)=ti
				END IF
			ENDIF
		
		ENDIF

		ti=ti+dti 

1     CONTINUE
		ti=ti-dti  !resting the extra dti no calculated
        dti=(dt*SIZE_T-ti)/(10*SIZE_2)

!*******************************************************
!	--------------------------
!	Calculating Second Window
!	--------------------------
!*******************************************************
      DO 2, j=1,11*SIZE_2
	    ti=ti+dti
		IF (.NOT. triangular) THEN
			!*****************************
			! Using Heidler Function
			!*****************************
			CALL INTEGRAL(IPULSE,0.D0,ti,RSLT,1.D-7,GNI)
			Icorr(nmaxi1+j)=RSLT
			tcorr(nmaxi1+j)=ti
		ELSE
			!*****************************
			! Integral for double ramp function . if triangular =.TRUE.
			!*******************************
			IF (ti .LE. tao11 ) THEN
				Icorr(nmaxi1+j)=.5*Ih1*ti**2/tao11
				tcorr(nmaxi1+j)=ti
			ELSEIF(ti .GT. tao11 .OR. ti .LT. 2*tao21-tao11) THEN
				Icorr(nmaxi1+j)=.5*Ih1*tao11 !+ Ih1*(ti-tao11)-Ih1*(.5*ti*ti-ti*tao11)/(2*tao21-tao11)
				tcorr(nmaxi1+j)=ti
			ELSE
				Icorr(nmaxi1+j)=.5*Ih1*tao11 !+ Ih1*(tao21-tao11)-Ih1*(.5*tao21*tao21-tao21*tao11)/(2*tao21-tao11)
				tcorr(nmaxi1+j)=ti
			END IF
		ENDIF

2     CONTINUE
!---------------------------------------------------------------------
! END CALTULATING Lightning Current
!*	******************************************************************

!************************************************+
! Script para grabar la corriente calculada
!************************************************

!	OPEN(UNIT=3, FILE='integrali.txt', ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
 !    WRITE(UNIT=3,FMT='(1X,E12.6E2,E12.6E2)')(tcorr(j),Icorr(j), j=1,(nmaxi1+10*SIZE_2+2))
!	CLOSE(3)

!*	--------------------------------------------------------------
!                  BESSEL Integrals Computation
!*	--------------------------------------------------------------
		IF (conductividad .EQ. 1) THEN
			p=o/(2*e0*e)
			t2=-dt
			DO 3 i=1,SIZE_T+2
				t2=t2+dt
				p12=p*(t2)
				i02(i)=BESSELI0(p12)
				i12(i)=BESSELI1(p12)
3			CONTINUE
		ENDIF
	ENDIF

!	OPEN(UNIT=4, FILE='bessel.txt', ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
!     WRITE(UNIT=4,FMT='(1X,E12.3E2,E12.3E2)')(i02(j),i12(j), j=1,(SIZE_T+2))
!	CLOSE(4)


!*************************
! READING EM FIELD FROM EXTERNAL FILE
!*************************************
IF (read_campo) THEN
    write(kmax_s,'(I3)') kmax+1
    FORMATO='(1X,'//kmax_s//'E16.7E4)'
    OPEN(UNIT=2, FILE='vertical.txt', STATUS='UNKNOWN')
    OPEN(UNIT=3, FILE='horizontal.txt', STATUS='UNKNOWN')
    IF (Imprimir_inf .EQV. .TRUE.) THEN
        write(*,*) '*INF* Total line division read=',kmax+1
        write(*,*) '*INF* Total time steps read =',SIZE_T+2
    ENDIF
    READ(UNIT=2,FMT='(1X,2E16.7E4)')(Evini(j),Evfin(j), j=1,SIZE_T+2)
    READ(UNIT=3,FMT=FORMATO)((Ex(cond,i,j),i=1,kmax+1),j=1,SIZE_T+2)

    CLOSE (2)
    CLOSE (3)
ELSE
!asignando valor de h para cada c�lculo de la integral importante para Componente horizontal


   DO 51,k = 0,kmax
	  DO 50,g=1,cond
	    
	    IF (status_perfil) THEN
	        h=h_perfil(k,g)
	    ELSE
		    h=hm(g)
		ENDIF
		
	  IF (XA.LT.XB) THEN
		IF (k.EQ.0) THEN
			r = DSQRT( ( XA ) ** 2 + (YR+Xi(g)) * (YR+Xi(g)) )
			ang = ( (XA ) )/r

		ELSEIF (k.EQ.kmax) THEN

			!        2*vl*dt represents the two small lines at ends
		
			r = DSQRT( ( XA + 2*dx2 + (kmax-2) * dx ) ** 2 + (YR+Xi(g)) * (YR+Xi(g))  )
			ang = ( (XA + 2*dx2 + (kmax-2) * dx) )/r

		ELSE
			!vl*dt represents the small line at the begining
			r = DSQRT( ( XA + dx2 + (k-1) * dx ) ** 2 + (YR+Xi(g))  * (YR+Xi(g))  )
			ang = ( (XA + dx2 + (k-1) * dx) )/r

		ENDIF	  
	   ELSE
		IF (k.EQ.0) THEN
			r = DSQRT( ( XA ) ** 2 + YR * YR )
			ang = ( (-XA ) )/r
	

		ELSEIF(k.EQ.kmax) THEN
			!2*vl*dt represents the two small lines at ends
			r = DSQRT( ( XA - 2*dx2 - (kmax-2) * dx ) ** 2 + (YR+Xi(g))  * (YR+Xi(g))  )
			ang = ( (- XA + 2*dx2 + (kmax-2) * dx) )/r
		
		
		ELSE

			r = DSQRT( ( XA - dx2 - (k-1) * dx ) ** 2 + (YR+Xi(g))  * (YR+Xi(g))  )
			ang = ( (-XA + dx2 + (k-1) * dx) )/r	  
		END IF
	   END IF
			t02=r/c; !initial time of the signal

		!SIZE_1= (tmax1+t02-t0)/dt			!considering t0 for initiation
		SIZE_1= (tmax1+t02-t0_min)/dt
		!SIZE_2=(tmax-(tmax1+t02-t0))/dt2	!considering t0 for initiation
		SIZE_2=(tmax-(tmax1+t02-t0_min))/dt2
		SIZE_2=SIZE_2+2
		IF (SIZE_2 .LT. 2) THEN
			SIZE_2=2
		END IF
		!IF ((SIZE_1*dt+t02-t0) .GT. tmax) THEN
		IF ((SIZE_1*dt+t02-t0_min) .GT. tmax) THEN
			SIZE_1=tmax/dt
			SIZE_2=2
		END IF
		  IF (1.005 .GT. abs(rant/r) .AND. abs(rant/r) .GT. .995 .AND. h .EQ. hm(1)) THEN
			n=nint(abs(r-rant)/c/dt)+1 
			Ex(g,k+1,n:SIZE_T+2)=Ex(g-1,k+1,1:SIZE_T+3-n)*abs(rant/r)/angant*ang
			!IF (k.EQ.0 ) THEN
			!  Evini(n:SIZE_T+2)=Evini(1:SIZE_T+3-n)*abs(rant/r)/angant*ang
			!ELSEIF (k .EQ. kmax) THEN
			!  Evfin(g,n:SIZE_T+2)=Evfin(g-1,1:SIZE_T+3-n)*abs(rant/r)/angant*ang
			!ENDIF
		  ELSE
			IF (k.EQ.0) THEN
				kcampo=1	!Identify if Ev will be calculated or not
				CALL EMFIELD(SIZE_T+2, EZ, EH,ang,Icorr,tcorr,i02,i12,kcampo,t0_min)
				Evini(:)=EZ(:)				
			ELSEIF (k.EQ.kmax) THEN
				kcampo=1
				CALL EMFIELD(SIZE_T+2, EZ, EH,ang,Icorr,tcorr,i02,i12,kcampo,t0_min)
				Evfin(:)=EZ(:)
			ELSE
				kcampo=0
				CALL EMFIELD(SIZE_T+2, Ev,EH,ang,Icorr,tcorr,i02,i12,kcampo,t0_min)
				
			ENDIF
			
			Ex(g,k+1,:)=EH(:)
		  ENDIF
			rant=r
			angant=ang
	
 IF (Imprimir_inf) THEN
    WRITE(*,*) '*INF* Point ',k+1, 'of ', kmax+1
	write(*,*) '*INF* r=',r,' ang=', ang
ENDIF

50    CONTINUE
		rant=0.D0
	!IF (k .EQ. int(.25*kmax) .OR. k .EQ. int(.5*kmax) 	&
	!	.OR. k .EQ. int(.75*kmax) .OR. k .EQ. kmax ) THEN
	!	write(*,*)'Percentage calculated for line',indarchiv,' : ',int(real(k)/real(kmax)*100),'%'
	!END IF

51 CONTINUE

    ENDIF

    IF (Imprimir_inf) THEN
        write(*,*)'*INF* EM Field calculated for line_',indarchiv
    ELSE
	    write(*,'(A)', advance="no")'*'
    ENDIF
!****************************************
!Routina para imprimir campo
!***************************************
ELSE
	Evini=0
	Evfin=0
	Ex=0
ENDIF
!READ EXTERNAL FILE WITH E_M FIELD


IF (Imprimir_campo) THEN
	IF (k.EQ.1 .OR. k .EQ. (kmax+1)) THEN
		OPEN(UNIT=3, FILE='campo'//indarchiv//'.txt', ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
		WRITE(UNIT=3,FMT='(1X,E14.6E2,E14.6E2,E14.6E2,E14.6E2)')(Evini(j),Evfin(j),Ex(cond,1,j),Ex(cond,kmax+1,j), j=1,SIZE_T)
!		OPEN(UNIT=3, FILE='campo'//indarchiv//'.txt', ACCESS='SEQUENTIAL',STATUS='UNKNOWN')

!		WRITE(UNIT=3)((Ex(cond,i,j),i=1,kmax),j=1,SIZE_T+2)
		CLOSE(3)
	ENDIF
ELSE
	write(*,'(A)', advance="no")'e'
ENDIF



t0=t0-t0_min

open(UNIT = 16, FILE = ncase_s//'_'//casename(1:LEN_TRIM(casename))//'_'//indarchiv//'.dat', FORM='UNFORMATTED', STATUS = 'UNKNOWN',ERR=1016,IOSTAT=ERRNUM)
write(UNIT = 16) Evini, Evfin, Ex, D3i,  D5,  Zc_o,Zci_o,Mv1  ,Mv2  ,Mv3  ,Mi1  ,Mi2  ,Mi3  ,Mi4

IF (Imprimir_inf) THEN
    write(*,*) '*INF* Electromanetic field was sucsesfully saved'
    write(*,*) '*INF* File: ',ncase_s//'_'//casename(1:LEN_TRIM(casename))//'_'//indarchiv//'.dat'
ENDIF
close (UNIT = 16)

open(UNIT = 17, FILE = 'bas_'//ncase_s//'_'//casename(1:LEN_TRIM(casename))//'_'//indarchiv//'.dat', FORM='UNFORMATTED', STATUS = 'UNKNOWN',ERR=1017,IOSTAT=ERRNUM)
write(UNIT = 17) v,zlam,Hcan,conductividad,o,e,dx,Ih1,Ih2,tao11,tao21,n1,n2,X0,Y0,rc,YA0,YB0,XA0,XB0,hm,Xi,kmax,t0
close (UNIT = 17)
				 
1016 SELECT CASE(ERRNUM>0)
		CASE (.TRUE.)
			write(*,*)'*ERROR* saving',casename(1:LEN_TRIM(casename))//indarchiv//'.dat'
		END SELECT 
1017 SELECT CASE(ERRNUM>0)
		CASE (.TRUE.)
			write(*,*)'*ERROR* saving ','bas_'//casename(1:LEN_TRIM(casename))//indarchiv//'.dat'
		END SELECT 



ELSE
IF (Imprimir_inf) THEN
    write(*,*)'*INF* Reading from calculated EMField_',indarchiv
ELSE
	write(*,'(A)', advance="no")'.'
ENDIF
if(status_hm) THEN


open(UNIT = 19, FILE = ncase_s//'_'//casename(1:LEN_TRIM(casename))//'_'//indarchiv//'.dat', FORM='UNFORMATTED', STATUS = 'UNKNOWN',ERR=1019,IOSTAT=ERRNUM)
write(UNIT = 19) Evini,Evfin,Ex,D3i,D5,Zc_o,Zci_o,Mv1,Mv2,Mv3,Mi1,Mi2,Mi3,Mi4
close (UNIT = 19)

open(UNIT = 18, FILE = 'bas_'//ncase_s//'_'//casename(1:LEN_TRIM(casename))//'_'//indarchiv//'.dat', FORM='UNFORMATTED', STATUS = 'UNKNOWN',ERR=1018,IOSTAT=ERRNUM)
write(UNIT = 18) v,zlam,Hcan,conductividad,o,e,dx,Ih1,Ih2,tao11,tao21,n1,n2,X0,Y0,rc,YA0,YB0,XA0,XB0,hm,Xi,kmax,t0
close (UNIT = 18)
1019 SELECT CASE(ERRNUM>0)
		CASE (.TRUE.)
			write(*,*)'*ERR* Error saving ',casename(1:LEN_TRIM(casename))//indarchiv//'.dat'
		END SELECT 
1018 SELECT CASE(ERRNUM>0)
		CASE (.TRUE.)
			write(*,*)'*ERR* Error saving','bas_'//casename(1:LEN_TRIM(casename))//indarchiv//'.dat'
		END SELECT 

ENDIF



ENDIF
IF (Imprimir_inf) THEN
    write(*,*) '*INF* Maximum number of lines: ',nlin
    
ENDIF
!********************************
! Definiendo Variables del Model
!--------------------------------
Xvar(1)= kmax											
Xvar(2)= SIZE_T 
Xvar(3)=dx											
Xvar(4)=o
Xvar(5)=e
Xvar(6)=conductividad
Xvar(7)=t0     !calculation of initial time
Xvar(8)=dx2
Xvar(9)=nlin !enviando n�mero m�ximo de l�neas
Xvar(10)=kmaxt
Xvar(11)=cond_max
IF (status_perfil.EQV. .TRUE.) THEN
    Xvar(12)=1
ELSE
    Xvar(12)=0
ENDIF

!Definiendo Variables de Salida
ind_ini1=21
Xvar(20)=ind_ini1
Xvar(ind_ini1:ind_ini1-1+cond)=hm(:)

		ind_ini=ind_ini1+cond
Xvar(50)=merge(1.0,0.,campo_distante);  !Variable of the distant lightning

!***********************************************************
!	DO i=1,cond
!		Xvar(ind_ini             +(i-1)*cond:ind_ini-1            +i*cond)=D3i(i,:) 
!		Xvar(ind_ini +  cond*cond+(i-1)*cond:ind_ini-1+  cond*cond+i*cond)=D5(i,:)  
!		Xvar(ind_ini +2*cond*cond+(i-1)*cond:ind_ini-1+2*cond*cond+i*cond)=Zc(i,:)  
!		Xvar(ind_ini +3*cond*cond+(i-1)*cond:ind_ini-1+3*cond*cond+i*cond)=Li(i,:)  
!		Xvar(ind_ini +4*cond*cond+(i-1)*cond:ind_ini-1+4*cond*cond+i*cond)=Ci(i,:)  
!		Xvar(ind_ini +5*cond*cond+(i-1)*cond:ind_ini-1+5*cond*cond+i*cond)=LCi(i,:)
!	ENDDO
!***********************************************************




!	ind_ini=ind_ini+6*cond*cond+6*cond


!	DO i=1,cond
!		Xvar(ind_ini +3*cond*(kmax+1)+2*cond*(SIZE_T+2)+(i-1)*(SIZE_T+2)     :(ind_ini-1)+3*cond*(kmax+1)+2*cond*(SIZE_T+2)+i*(SIZE_T+2))=Evini(i,:)
!		Xvar(ind_ini +3*cond*(kmax+1)+3*cond*(SIZE_T+2)+(i-1)*(SIZE_T+2)     :(ind_ini-1)+3*cond*(kmax+1)+3*cond*(SIZE_T+2)+i*(SIZE_T+2))=Evfin(i,:)
!		write(*,*)'resta=',ind_ini +3*cond*(kmax+1)+4*cond*(SIZE_T+2)+(kmax+1)*(i-1)*(SIZE_T+2)+(k-1)*(SIZE_T+2)-((ind_ini-1)+3*cond*(kmax+1)+4*cond*(SIZE_T+2)+(kmax+1)*(i-1)*(SIZE_T+2)+k*(SIZE_T+2))
!		Do k=1,kmax+1
!		Xvar(ind_ini +3*cond*(kmax+1)+4*cond*(SIZE_T+2)+(kmax+1)*(i-1)*(SIZE_T+2)+(k-1)*(SIZE_T+2):(ind_ini-1)+3*cond*(kmax+1)+4*cond*(SIZE_T+2)+(kmax+1)*(i-1)*(SIZE_T+2)+k*(SIZE_T+2))=Ex(i,k,:)
		
!		End DO
!	ENDDO
DEALLOCATE (Zc,Zci,Ex,EH,EZ,Ev,Evfin,Evini,Cm,Lm,Ci,Li,LCi,LC,D1,D3i,D2, &
	D5,Ro,Rf,hm,Xi,hm_s,Xi_s,Rm,D7,D6,D3,STAT=ALLOC_ERR)

END SUBROUTINE

!******************************************************************************************
!				CALCULATING LINE PARAMETERS
!******************************************************************************************
SUBROUTINE LINEPARAM(cond,Xi,hm,rc,Lm,Cm,Zc)
	INTEGER cond,i,j
    DOUBLE PRECISION Lm(cond,cond),Cm(cond,cond),Zc(cond,cond),Xi(cond),hm(cond),s,e0,pi,mu,c,rc

	
	PARAMETER (pi=3.141592653589793D0)
	PARAMETER (mu=(4*pi)*1E-7)
	PARAMETER (c=2.99792458E8)
	PARAMETER (e0=1/(mu*c**2))


	DO 2, i=1,cond

      DO 1, j=1,cond

            s=((hm(i)-hm(j))**2+(Xi(i)-Xi(j))**2)**.5 !computing the space between conductors
            IF (i.EQ.j) THEN
               Lm(i,i)=mu/(2*pi)*log(2*hm(i)/rc)
            ELSE
               Lm(i,j)=mu/(4*pi)*log(1+4*hm(i)*hm(j)/(s*s))
            END IF

1         ENDDO

2     ENDDO

	  
      CALL inv(Lm,cond,cond,Cm)

	Cm=Cm*mu*e0
	Zc=Lm/(mu*e0)**0.5 !Calculating surge impedance matrix
	
	RETURN      
END SUBROUTINE



SUBROUTINE OPEN_LINE_FILE(rc,XA0,YA0,XB0,YB0,hm,Rm,Xi,cond,indarchiv)
!DEC$ ATTRIBUTES, DLLEXPORT::OPEN_LINE_FILE

INTEGER cond, ERRNUM,i
DOUBLE PRECISION rc,XA0,YA0,XB0,YB0,Rm(cond,cond),hm(cond),Xi(cond)

CHARACTER (LEN=3) indarchiv
    Rm(:,:)=0
    hm(:)=0
    Xi(:)=0
    
  OPEN (UNIT = 10, FILE = 'linea_'//indarchiv//'.txt', STATUS = 'OLD',ERR=1010,IOSTAT=ERRNUM)
	  READ   (10, *) rc
	  READ   (10, *) XA0,YA0
      READ   (10, *) XB0,YB0
      READ   (10, *) hm(1),Rm(1,1)
	  READ   (10, *) Xi(1)
        Rm(1,1)=Rm(1,1)/1000
!	  READ   (10, *) Rm(i,i) !Conductor Resistance ohm/km
!		Rm(1,1)=0.000
	  IF (cond .GT. 1) THEN
		DO i=2,cond
			READ   (10, *) hm(i),Rm(i,i)  !Heigh for multiple conductors
			READ   (10, *) Xi(i)  !Horizontal Relative distance for multiple conductors
            Rm(i,i)=Rm(i,i)/1000 !Considering Resistance per meter it is read as resistance per km
		END DO

	  ENDIF
1010 SELECT CASE(ERRNUM>0)
	    CASE (.TRUE.)
			write(*,*) '******************************************************************'
			write(*,*) '*ERR* Reading File: linea_'//indarchiv//'.txt'
            write(*,*) '*ERR* ERRNUM=',ERRNUM
            write(*,*) '*PROGRAM WILL STOP'
            write(*,*) '******************************************************************'
            STOP
		ENDSELECT 

  CLOSE (UNIT=10)
END SUBROUTINE



SUBROUTINE OPEN_MISC_FILE (v,zlam,Hcan,tmax1,SIZE_2,conductividad,o,e,dx,triangular,dist_camp)
!DEC$ ATTRIBUTES, DLLEXPORT::OPEN_MISC_FILE

DOUBLE PRECISION v,zlam,Hcan,tmax1,o,e,dx,dist_camp
INTEGER i, ERRNUM, SIZE_2,conductividad
LOGICAL triangular

      OPEN (UNIT = 8, FILE = 'miscelaneo.txt', STATUS = 'OLD', ERR=1008,IOSTAT=ERRNUM)
      READ   (8, *) v
	  READ   (8, *) zlam
      READ   (8, *) Hcan
      READ   (8, *) tmax1
      READ   (8, *) SIZE_2
      READ   (8, *) conductividad
      READ   (8, *) o
      READ   (8, *) e
	  READ   (8, *) dx
	  i=0
	  dist_camp=5000
	  READ   (8, *,end=990) i
	  READ   (8, *,end=990) dist_camp 

990 CONTINUE
	  
	  

	  IF (i .EQ. 1) THEN 
		triangular=.TRUE.
	  ELSE
		triangular=.FALSE.
	  ENDIF

		
1008 SELECT CASE(ERRNUM>0)
	    CASE (.TRUE.)
			write(*,*) '******************************************************************'
			write(*,*) '*WARN* Error opening file miscelaneo.txt. Error_number:',ERRNUM
			write(*,*) '*INF* Loading default values '
			write(*,*) '******************************************************************'
			write(*,*) '*********************'
			write(*,*) '* v = 1.3D8	        *'	
			write(*,*) '* zlam = 2000	    *'	
			write(*,*) '* Hcan = 8000	    *'	
			write(*,*) '* tmax1= 3.d-6      *'	
			write(*,*) '* SIZE_2=30	        *'	
			write(*,*) '* conductividad=1   *'	
			write(*,*) '* o=1.D-3		    *'
			write(*,*) '* e=10		        *'	
			write(*,*) '* dx=10		        *'
			write(*,*) '*********************'
			
    		v = 1.3D8	!	- return stroke velocity
	    	zlam = 2000		!- Atenuation for MTL model use 0 for TL model
		    Hcan = 8000		!- Channel Height
    		tmax1= 3.d-6	!- simulation time for first window
	    	SIZE_2=30		!- Number of samples for the second window
		    conductividad=1	!- Use 1(yes) to consider Conductivity
    		o=1.D-3			!- Ground Conductivity
	    	e=10			!- Relative ground permitivity
		    dx=10			!- space division if doesn't fulfill Courant criteria it is atumatically computed
!		----------------------------------------------------------------
 
		END SELECT 

CLOSE (UNIT=8)

END SUBROUTINE

SUBROUTINE OPEN_CORR_FILE(Ih1,Ih2,tao11,tao21,tao12,tao22,n1,n2,X0,Y0,ncase_s)
!DEC$ ATTRIBUTES, DLLEXPORT::OPEN_CORR_FILE


DOUBLE PRECISION Ih1,Ih2,tao11,tao21,tao12,tao22,n1,n2,X0,Y0
CHARACTER (LEN=5) ncase_s
INTEGER ERRNUM

 OPEN (UNIT = 9, FILE = 'corr_'//ncase_s//'.txt', STATUS = 'OLD', ERR=1009,IOSTAT=ERRNUM)

      READ   (9, *) Ih1
      READ   (9, *) Ih2
      READ   (9, *) tao11
      READ   (9, *) tao21
      READ   (9, *) tao12
      READ   (9, *) tao22
      READ   (9, *) n1
      READ   (9, *) n2
      READ   (9, *) X0
      READ   (9, *) Y0      

1009 SELECT CASE(ERRNUM>0)
	    CASE (.TRUE.)
			write(*,*) '******************************************************************'
			write(*,*) '*WARN* Error reading corr_'//ncase_s//'.txt file, Error_number=',ERRNUM
			write(*,*) '*INF* All values will set in zero, EM field will not be calculated'
			write(*,*) '*INF* EM field will be zero, then, just line model is used'
			write(*,*) '******************************************************************'
            Ih1=0
            Ih2=0
            tao11=0
            tao21=0
            tao12=0
            tao22=0
            n1=0
            n2=0
            X0=0
            Y0=0    
		END SELECT 

  CLOSE (UNIT=9)

  END SUBROUTINE
  
 !***********************************************
 ! opening and reading Profile file for each line
 !***********************************************
 
  SUBROUTINE OPEN_LINE_PROFILE(h_perfil,kmax,cond,indarchiv)
  INTEGER kmax,cond,ERRNUM,i
  DOUBLE PRECISION h_perfil(kmax,cond)
    CHARACTER (LEN=3) indarchiv

    !DEC$ ATTRIBUTES, DLLEXPORT::OPEN_LINE_PROFILE
    OPEN (UNIT = 101, FILE = 'line_profile_'//indarchiv//'.txt', STATUS = 'OLD', ERR=1009,IOSTAT=ERRNUM)
    DO i=1,kmax
        READ   (101, *) h_perfil(i,:)
    ENDDO
    write(*,*) '*INF* line profile, read succesfully',h_perfil(kmax,:),' kmax=',kmax
    
	
1009 SELECT CASE(ERRNUM>0)
	    CASE (.TRUE.)
			write(*,*) '******************************************************************'
			write(*,*) '*WARN* Error reading line_profile_'//indarchiv//'.txt'
            write(*,*) '*INF* height will be set in default values h=10'
			write(*,*) '******************************************************************'
			h_perfil(i,:)=10
		END SELECT 

  END SUBROUTINE