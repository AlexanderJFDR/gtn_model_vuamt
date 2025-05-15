!求解器设置双精度
!**********************************************************************************************************
!
      module NumKind
!
!**********************************************************************************************************
        implicit none
        integer (kind(1)), parameter :: ikind = kind(1), 
     &                                  rkind = kind(0.D0), 
     &                                  lkind = kind(.true.)
!       
      end module Numkind
      

!**********************************************************************************************************
!
      module ModelParam
!
!**********************************************************************************************************
        use NumKind
        implicit none

        ! Constants

        ! Flag of initilization
        logical(lkind) :: bInitialized = .false.

		! Max iteraition
		real(ikind), parameter :: Maxiter = 100

        ! Tolerance
        real(rkind), parameter :: TOL = 1.0d-8 
        
        ! pi parameter
        real(rkind), parameter :: pi = 3.1415926535897932384626433832d0
        
        ! 
        real(rkind) :: E, nu, sigma_0, Q_0, b_0, q1, q2, q3, f_0 
	 	real(rkind) :: f_c, f_F, S_N, f_N, epsilon_N
	 
        real(rkind) :: G, K, lamda, mu       

        !
        contains

        !===================================
          subroutine Initialize(props, nprops)
        !===================================

            integer (ikind), intent (in) :: nprops
            real    (rkind), intent (in) :: props(nprops)

            ! material properties
            E         =  props(1) ! props(1) -- Young's modulus
            nu        =  props(2) ! props(2) -- Poisson's ratio
            sigma_0   =  props(3) ! props(3) -- Initial yield stress
			Q_0       =  props(4) ! props(4) -- Harden parameter
			b_0       =  props(5) ! props(5) -- Harden exponent
            q1        =  props(6) ! props(6) -- GTN model 1st constant: q_1
			q2        =  props(7) ! props(7) -- GTN model 2st constant: q_2
			q3        =  props(8) ! props(8) -- GTN model 3st constant: q_3
            f_0       =  props(9) ! props(9) -- Initial void content
			f_c       =  props(10) ! props(10) -- Critical void content
			f_F       =  props(11) ! props(11) -- Failure void content
			S_N       =  props(12) ! props(12) -- Standard deviation in nucleation
			f_N       =  props(13) ! props(13) -- void content in nucleation
			epsilon_N =  props(14) ! props(14) -- mean plastic strain in nucleation
            
			! elastic parameters
            G       = E / (2.d0*(1.d0 + nu))
            K       = E / (3.d0*(1.d0 - 2.d0*nu))
            lamda   = E*nu / ((1.d0 + nu) * (1 - 2.d0*nu))
			mu      = E / (2.d0*(1.d0 + nu))
			
            bInitialized = .true.
            
            return
          end subroutine Initialize
      !========================================================================= 
      end module ModelParam
	  
!**********************************************************************************************************
	  !Return mapping algorithm
      subroutine rma(p_trial, sigma_q_trial, epsilon_m_p_trial, f_trial, 
     &      sigma_y_trial, delta_dp, delta_dq, flag)
!**********************************************************************************************************
		
		use NumKind
        use ModelParam
        implicit none
		
!**********************************************************************************************************		
		!define varibles
		real(rkind):: p_trial, sigma_q_trial, epsilon_m_p_trial, f_trial, 
     &      sigma_y_trial, delta_dp, delta_dq 
	 
		real(rkind):: p, sigma_q, epsilon_m_p, f, d_delta_dp, d_delta_dq,
     &      sigma_y, delta_epsilon_m_p, delta_f, A, f_star, f_F_bar, phi
		
		real(rkind):: A11, A12, A21, A22, B1, B2, C11, C12,
     &      C21, C22, D1, D2, E1, E2, det, det1, det2, res	
		
		real(rkind):: D_phi_sigmaq, D_phi_p, D_phi_epsilonmp,
     &      D_phi_f, D2_phi_p_epsilonmp, D2_phi_p_f, D2_phi_sigmaq_epsilonmp,
     &      D2_phi_sigmaq_f, D2_phi_p, D2_phi_sigmaq
	 
		real(rkind):: D_epsilonmp_deltadp, D_f_deltadp, D_epsilonmp_deltadq, D_f_deltadq
	 
		real(rkind):: D_deltaepsilonmp_epsionmp, D_deltaepsilonmp_f,
     &      D_deltaf_epsionmp, D_deltaf_f, D_deltaepsilonmp_deltadp,
     &      D_deltaf_deltadp, D_deltaepsilonmp_deltadq, D_deltaf_deltadq,
     &      D_deltaepsilonmp_p, D_deltaf_p, D_deltaepsilonmp_sigmaq,
     &      D_deltaf_sigmaq, D_sigmay_epsilonmp
	 
		integer(ikind):: iter, flag
	 
		!prepare varibles
		flag = 0
		p = p_trial
		sigma_q = sigma_q_trial
		epsilon_m_p = epsilon_m_p_trial
		f = f_trial
		
		if (f <= f_c) then
		  f_star = f
		else if ( f > f_c .and. f < f_F ) then
		  f_F_bar = (q1+sqrt(q1**2.d0+q3))/q3
          f_star = f_c + (f_F_bar-f_c)/(f_F-f_c)*(f-f_c)
		else 
		  f_star = f_F_bar
		end if  
		
		sigma_y = sigma_y_trial
		!Compute Yield Function
		phi = (sigma_q_trial**2.d0/sigma_y**2.d0)
     &		+ (2.d0*q1*f_star*cosh(-3.d0*q2*p_trial/2.d0/sigma_y))
     &    	- (1.d0+q3*f_star**2.d0)
	 
		delta_dp = 0.d0
		delta_dq = 0.d0
		iter = 0
		
		!solve delta_dp, delta_dq
		do iter = 1,Maxiter
		  delta_epsilon_m_p = (sigma_q*delta_dq-p*delta_dp)/(1.d0-f)/sigma_y
		  A = f_N/S_N/sqrt(2.d0*pi)*exp(-((epsilon_m_p-epsilon_N)/S_N)**2.d0/2.d0)
		  if (p >= 0.d0) then
			delta_f = A*delta_epsilon_m_p + (1.d0-f)*delta_dp
		  else
			delta_f = (1.d0-f)*delta_dp
		  end if
		  call Voce(sigma_y, D_sigmay_epsilonmp, epsilon_m_p)
		  
		!solve D_epsilonmp_deltadp, D_f_deltadp, D_epsilonmp_deltadq, D_f_deltadp
		  D_deltaepsilonmp_epsionmp = -delta_epsilon_m_p/sigma_y*D_sigmay_epsilonmp
		  D_deltaepsilonmp_f = delta_epsilon_m_p/(1.d0-f)
		  if (p >= 0.d0) then
			D_deltaf_epsionmp = A*(epsilon_m_p-epsilon_N)/S_N**2.d0*delta_epsilon_m_p
		  else
			D_deltaf_epsionmp = 0.d0
		  end if
		  D_deltaf_f = -delta_dp
		  D_deltaepsilonmp_deltadp = -p/(1.d0-f)/sigma_y
		  D_deltaf_deltadp = 1.d0-f
		  D_deltaepsilonmp_deltadq = sigma_q/(1.d0-f)/sigma_y
		  D_deltaf_deltadq = 0.d0
		  D_deltaepsilonmp_p = -delta_dp/(1.d0-f)/sigma_y
		  D_deltaf_p = 0.d0
		  D_deltaepsilonmp_sigmaq = delta_dp/(1.d0-f)/sigma_y
		  D_deltaf_sigmaq = 0.d0
		
		  C11 = 1.d0-D_deltaepsilonmp_epsionmp
		  C12 = -D_deltaepsilonmp_f
		  C21 = -D_deltaf_epsionmp
		  C22 = 1.d0-D_deltaf_f
		  
		  D1 = D_deltaepsilonmp_deltadp + K*D_deltaepsilonmp_p
		  D2 = D_deltaf_deltadp + K*D_deltaf_p
		  
		  E1 = D_epsilonmp_deltadq - 3.d0*G*D_deltaepsilonmp_sigmaq
		  E2 = D_deltaf_deltadq - 3.d0*G*D_deltaf_sigmaq
		
		  det = C11*C22 - C12*C21
		  det1 = D1*C22 - D2*C12
		  det2 = C11*D2 - C21*D1
		  D_epsilonmp_deltadp = det1/det
		  D_f_deltadp = det2/det
		  
		  det1 = E1*C22 - E2*C12
		  det2 = C11*E2 - C21*E1
		  D_epsilonmp_deltadq = det1/det
		  D_f_deltadq = det2/det
		  
		!solve d_delta_dp, d_delta_dq
		  D_phi_sigmaq = 2*sigma_q/sigma_y**2.d0
		  D_phi_p = -3.d0*q1*q2*f_star/sigma_y*sinh(-3.d0*q2*p/2.d0/sigma_y)
		  D2_phi_sigmaq = 2.d0/sigma_y**2.d0
		  D2_phi_p = 9.d0*q1*q2**2.d0*f_star/2.d0/sigma_y**2.d0 
     &   	  * cosh(-3.d0*q2*p/2.d0/sigma_y)
		  D_phi_epsilonmp = (-2.d0*sigma_q**2.d0/sigma_y**3.d0 
     &    	  + 3.d0*q1*q2*p/sigma_y**2.d0*f_star*sinh(-3.d0*q2*p/2.d0/sigma_y))
     &        * D_sigmay_epsilonmp
	 
		  if (f <= f_c) then
	        D_phi_f = 2.d0*q1*cosh(-3.d0*q2*p/2.d0/sigma_y)-2.d0*q3*f_star
		  else if ( f > f_c .and. f < f_F ) then
		    D_phi_f = (2.d0*q1*cosh(-3.d0*q2*p/2.d0/sigma_y)-2.d0*q3*f_star)
     &        * (f_F_bar-f_c)/(f_F-f_c)	 
		  else 
		    D_phi_f = 0.d0
		  end if
		  
		  D2_phi_p_epsilonmp = (3.d0*q1*q2/sigma_y**2*f_star*sinh(-3.d0*q2*p/2.d0/sigma_y)
     &      - 9*q1*q2**2.d0*p/2.d0/sigma_y**3.d0*f_star*cosh(-3.d0*q2*p/2.d0/sigma_y))*D_sigmay_epsilonmp	  
		  
		  if (f <= f_c) then
	        D2_phi_p_f = -3.d0*q1*q2/sigma_y*sinh(-3.d0*q2*p/2.d0/sigma_y)
		  else if ( f > f_c .and. f < f_F ) then
		    D_phi_f = (-3.d0*q1*q2/sigma_y*sinh(-3.d0*q2*p/2.d0/sigma_y))
     &        * (f_F_bar-f_c)/(f_F-f_c)	 
		  else 
		    D_phi_f = 0.d0
		  end if
		  
		  D2_phi_sigmaq_epsilonmp = -4.d0*sigma_q/sigma_y**3.d0*D_sigmay_epsilonmp
		  D2_phi_sigmaq_f = 0.d0
		  
		  A11 = D_phi_sigmaq + delta_dq*(K*D2_phi_p+D2_phi_p_epsilonmp 
     &   	  * D_epsilonmp_deltadp+D2_phi_p_f*D_f_deltadp) + delta_dp
     &        * (D2_phi_sigmaq_epsilonmp*D_epsilonmp_deltadp+D2_phi_sigmaq_f*D_f_deltadp)
	      A12 = D_phi_p + delta_dq*(D2_phi_p_epsilonmp*D_epsilonmp_deltadq
     &        + D2_phi_p_f*D_f_deltadq) + delta_dp*(-3.d0*G*D2_phi_sigmaq 
     &        + D2_phi_sigmaq_epsilonmp*D_epsilonmp_deltadq + D2_phi_sigmaq_f 
     &        * D_f_deltadq)
		  A21 = K*D_phi_p + D_phi_epsilonmp*D_epsilonmp_deltadp 
     &        + D_phi_f*D_f_deltadp
		  A22 = -3.d0*G*D_phi_sigmaq + D_phi_epsilonmp*D_epsilonmp_deltadq
     &        + D_phi_f*D_f_deltadq
		  
		  B1 = -delta_dq*D_phi_p - delta_dp*D_phi_sigmaq
		  B2 = -phi
		  
		  det = A11*A22 - A12*A21
		  det1 = B1*A22 - B2*A12
		  det2 = A11*B2 - A21*B1
		  d_delta_dp = det1/det
		  d_delta_dq = det2/det
		  
		!Update Values  
		  delta_dp = delta_dp + d_delta_dp
		  delta_dq = delta_dq + d_delta_dq
		  
		  p = p_trial+K*delta_dp
		  sigma_q = sigma_q_trial-3.d0*G*delta_dq
		  
		  delta_epsilon_m_p = (sigma_q*delta_dq-p*delta_dp)/(1.d0-f)/sigma_y
		  epsilon_m_p = epsilon_m_p_trial + delta_epsilon_m_p
		  
		  A = f_N/S_N/sqrt(2.d0*pi)*exp(-((epsilon_m_p-epsilon_N)/S_N)**2.d0/2.d0)
		  if (p >= 0) then
			delta_f = A*delta_epsilon_m_p+(1.d0-f)*delta_dp
		  else
			delta_f = (1.d0-f)*delta_dp
		  end if
		  f = f_trial + delta_f

		  if (f <= f_c) then
		  f_star = f
		  else if ( f > f_c .and. f < f_F ) then
		    f_F_bar = (q1+sqrt(q1**2.d0+q3))/q3
            f_star = f_c + (f_F_bar-f_c)/(f_F-f_c)*(f-f_c)
		  else 
		    f_star = f_F_bar
		  end if 

		  call Voce(sigma_y, D_sigmay_epsilonmp, epsilon_m_p)
		  
		  phi = (sigma_q**2.d0/sigma_y**2.d0)
     &		  + (2.d0*q1*f_star*cosh(-3.d0*q2*p/2.d0/sigma_y))
     &    	  - (1.d0+q3*f_star**2.d0)
	 
		!Compute convergence criterion
          res = abs(phi)
		  
		!Check for convergence
          if(res < TOL) then
		    flag = 1
			exit
		  end if
		end do
	  
		return
	 
	  end subroutine rma
	  
!**********************************************************************************************************
	  !Voce model
      subroutine Voce(sigma_y, D_sigmay_epsilonmp, epsilon_m_p)	  
	  
	    use NumKind
        use ModelParam
        implicit none
		
		real(rkind):: sigma_y, D_sigmay_epsilonmp, epsilon_m_p
		
		sigma_y = sigma_0 + Q_0*(1.d0-exp(-b_0*epsilon_m_p))
		D_sigmay_epsilonmp = Q_0*b_0*exp(-b_0*epsilon_m_p)
		
		return 
	  
	  end subroutine Voce
	  
!**********************************************************************************************************
      subroutine vumat(nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     &     			   stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     &                 props, density, strainInc, relSpinInc,
     &                 tempOld, stretchOld, defgradOld, fieldOld,
     &                 stressOld, stateOld, enerInternOld, enerInelasOld,
     &                 tempNew, stretchNew, defgradNew, fieldNew,
     &                 stressNew, stateNew, enerInternNew, enerInelasNew)
!**********************************************************************************************************

        use NumKind
        use ModelParam
        implicit none

!**********************************************************************************************************
      ! interface of vumat, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! variables passed in
		
        integer(ikind), intent (in    ) :: nblock, ndir, nshr, nstatev,
     &    nfieldv, nprops
     
        integer(ikind), intent (in    ) :: jInfoArray(*)
     
        real(rkind), intent (in    ) :: stepTime, totalTime
		
		
		real(rkind), intent (in    ) :: dtArray(*), coordMp(nblock),
     &    	charLength(nblock), props(nprops), density(nblock), strainInc(nblock, ndir+nshr),
     &      relSpinInc(nblock, nshr), tempOld(nblock), stretchOld(nblock, ndir+nshr),
     &      defgradOld(nblock, ndir+2*nshr), fieldOld(nblock,nfieldv), stressOld(nblock, ndir+nshr),
     &      enerInternOld(nblock), enerInelasOld(nblock), tempNew(nblock),
     &      stretchNew(nblock, ndir+nshr), defgradNew(nblock, ndir+2*nshr),
     &      fieldNew(nblock, nfieldv)
		
		character(len=8) :: cmname
  
        ! variables to be updated (the update of energy(8) is optional)
		
        real(rkind), intent (in out) :: stressNew(nblock, ndir+nshr), stateOld(nblock, nstatev), 
     &      stateNew(nblock, nstatev), enerInternNew(nblock), enerInelasNew(nblock)
	 
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !********************************************************************************************************
      !                   
      ! user defined varibles
        real(rkind):: traceInc, p_trial, sigma_q_trial, f_trial, 
     &      epsilon_m_p_trial, sigma_y_trial, delta_dp, delta_dq
	 
		real(rkind):: p, sigma_q, f, f_star, epsilon_m_p, delta_epsilon_m_p, 
     &	sigma_y, delta_f
	 
        real(rkind):: phi, f_F_bar, A, D_sigmay_epsilonmp, stressPower, stressPower_plastic
		
		real(rkind):: stress_trial(ndir+nshr), S_trial(ndir+nshr), delta_epsilon_p(ndir+nshr)
		
		integer(ikind):: i, flag

        ! initialize parameters, etc.                     
        if (.not. bInitialized) call Initialize(props, nprops)
        
		!********************************************************************************************************
        ! Calculate the preprocessor step
		if ((totalTime == 0.0).and.(stepTime == 0.0)) then
		  do i = 1, nblock
		! update stressNew
            traceInc = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
            stressNew(i,1) = stressOld(i,1)
     &          + 2.d0 * mu * strainInc(i,1) + lamda * traceInc
            stressNew(i,2) = stressOld(i,2)
     &          + 2.d0 * mu * strainInc(i,2) + lamda * traceInc
            stressNew(i,3) = stressOld(i,3)
     &          + 2.d0 * mu * strainInc(i,3) + lamda * traceInc
            stressNew(i,4) = stressOld(i,4) + 2.d0 * mu * strainInc(i,4)
            if (nshr > 1) then
              stressNew(i,5) = stressOld(i,5) + 2.d0 * mu * strainInc(i,5)
              stressNew(i,6) = stressOld(i,6) + 2.d0 * mu * strainInc(i,6)
            end if
		! update stateNew	
			stateOld(i,1) = f_0  !f
			stateOld(i,2) = 0.d0 !epsilon_m_p
			stateOld(i,3) = 1.d0 !is failure 0-yes, 1-no
          end do
		  
		!********************************************************************************************************
        else
		  do i = 1, nblock
		! compute stress_trial
			traceInc = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
            stress_trial(1) = stressOld(i,1)
     &        	+ 2.d0 * mu * strainInc(i,1) + lamda * traceInc
			stress_trial(2) = stressOld(i,2)
     &        	 + 2.d0 * mu * strainInc(i,2) + lamda * traceInc
			stress_trial(3) = stressOld(i,3)
     &        	+ 2.d0 * mu * strainInc(i,3) + lamda * traceInc
			stress_trial(4) = stressOld(i,4) + 2.d0 * mu * strainInc(i,4)
            if (nshr > 1) then
			  stress_trial(5) = stressOld(i,5) + 2.d0 * mu * strainInc(i,5)
			  stress_trial(6) = stressOld(i,6) + 2.d0 * mu * strainInc(i,6)
            end if
			
		!load parameters used in phi
		  !Calculate f_star
			f_trial = stateOld(i, 1)
			if (f_trial <= f_c) then
		    f_star = f_trial
			else if ( f_trial > f_c .and. f_trial < f_F ) then
			  f_F_bar = (q1+sqrt(q1**2.d0+q3))/q3
              f_star = f_c + (f_F_bar-f_c)/(f_F-f_c)*(f_trial-f_c)
			else 
			  f_star = f_F_bar
			end if         
		  !Calculate p_trial
			p_trial = -(stress_trial(1)+stress_trial(2)+stress_trial(3))/3.d0
		  !Calculate sigma_q_trial
			S_trial(1) = stress_trial(1) + p_trial
			S_trial(2) = stress_trial(2) + p_trial
			S_trial(3) = stress_trial(3) + p_trial
			S_trial(4) = stress_trial(4)
			if (nshr > 1) then
			  S_trial(5) = stress_trial(5)
			  S_trial(6) = stress_trial(6)
			end if      
			if (nshr == 1) then
			  sigma_q_trial = sqrt((3.0/2.0)*(S_trial(1)*S_trial(1)
     &            +S_trial(2)*S_trial(2)
     &            +S_trial(3)*S_trial(3))
     &            +3.d0*(S_trial(4)*S_trial(4)))
			else 
			  sigma_q_trial = sqrt((3.0/2.0)*(S_trial(1)*S_trial(1)
     &            +S_trial(2)*S_trial(2)
     &            +S_trial(3)*S_trial(3))
     &            +3.d0*(S_trial(4)*S_trial(4)
     &            +S_trial(5)*S_trial(5)
     &            +S_trial(6)*S_trial(6)))        
			end if
			
		  !Calculate epsilon_m_p_trial
			epsilon_m_p_trial = stateOld(i,2)
		  !Calculate sigma_y_trial
			call Voce(sigma_y_trial, D_sigmay_epsilonmp, epsilon_m_p_trial)
			
		!Compute Yield Function
			phi = (sigma_q_trial**2.d0/sigma_y_trial**2.d0)
     &			+(2.d0*q1*f_star*cosh(-3.d0*q2*p_trial/2.d0/sigma_y_trial))
     &    		-(1.d0+q3*f_star**2.d0)
		
		!********************************************************************************************************
		!Check for plasticity
			if(phi > 0.d0) then ! Plastic flow
		  !Initialize Variables
			  delta_dp = 0.d0
			  delta_dq = 0.d0
			  
			  call rma(p_trial, sigma_q_trial, epsilon_m_p_trial, f_trial, 
     &		      sigma_y_trial, delta_dp, delta_dq, flag)
			  
		  !Check for convergence
			  if (flag == 0) then
			    print*, "RMAP has not converged!"
			    stop
			  else
		      !Update Values  
				p = p_trial+K*delta_dp
				sigma_q = sigma_q_trial-3.d0*G*delta_dq
				  
				delta_epsilon_m_p = (sigma_q*delta_dq-p*delta_dp)/(1.d0-f_trial)/sigma_y_trial
				epsilon_m_p = epsilon_m_p_trial + delta_epsilon_m_p
				
				delta_epsilon_p(1) = delta_dq*(3.d0*S_trial(1)/2.d0/sigma_q_trial) + 1.d0/3.d0*delta_dp
                delta_epsilon_p(2) = delta_dq*(3.d0*S_trial(2)/2.d0/sigma_q_trial) + 1.d0/3.d0*delta_dp
                delta_epsilon_p(3) = delta_dq*(3.d0*S_trial(3)/2.d0/sigma_q_trial) + 1.d0/3.d0*delta_dp
                delta_epsilon_p(4) = delta_dq*(3.d0*S_trial(4)/2.d0/sigma_q_trial)
                if (nshr > 1) then
                  delta_epsilon_p(5) = delta_dq*(3.d0*S_trial(5)/2.d0/sigma_q_trial)
                  delta_epsilon_p(6) = delta_dq*(3.d0*S_trial(6)/2.d0/sigma_q_trial)
                end if
				  
				A = f_N/S_N/sqrt(2.d0*pi)*exp(-((epsilon_m_p-epsilon_N)/S_N)**2.d0/2.d0)
				if (p >= 0) then
			      delta_f = A*delta_epsilon_m_p+(1.d0-f)*delta_dp
				else
				  delta_f = (1.d0-f)*delta_dp
				end if
				f = f_trial + delta_f
				
			  !Update the stress tensor
                stressNew(i,1) = (stress_trial(1)+p_trial)*(sigma_q/sigma_q_trial) - p
                stressNew(i,2) = (stress_trial(2)+p_trial)*(sigma_q/sigma_q_trial) - p
                stressNew(i,3) = (stress_trial(3)+p_trial)*(sigma_q/sigma_q_trial) - p
                stressNew(i,4) = stress_trial(4)*(sigma_q/sigma_q_trial)
                if (nshr > 1) then
                  stressNew(i,5) = stress_trial(5)*(sigma_q/sigma_q_trial)
                  stressNew(i,6) = stress_trial(6)*(sigma_q/sigma_q_trial)
                end if
			  
			  !Update the history variables
                stateNew(i,1) = f             !f
                stateNew(i,2) = epsilon_m_p   !epsilon_m_p
				
				if (f > f_F) then
                  stateNew(i,3) = 0.d0         !is failure
				else
				  stateNew(i,3) = 1.d0   
				end if
			  end if
			  
		  !********************************************************************************************************  
		  !elastic not plastic
			else
			!Update the stress tensor
			  stressNew(i,1)   = stress_trial(1)
              stressNew(i,2)   = stress_trial(2)
              stressNew(i,3)   = stress_trial(3)
              stressNew(i,4)   = stress_trial(4)
              if (nshr > 1) then
                stressNew(i,5) = stress_trial(5)
                stressNew(i,6) = stress_trial(6)
              end if
			!Update the history variables
              stateNew(i,1) = f_trial            !f
              stateNew(i,2) = epsilon_m_p_trial  !epsilon_m_p
              if (f_trial > f_F) then
                stateNew(i,3) = 0.d0             !is failure 0-yes, 1-no
			  else
				stateNew(i,3) = 1.d0   
		      end if
			end if
			
		!********************************************************************************************************  	
		!Store the specific internal energy -
			if (nshr == 1) then
			  stressPower = dot_product(0.5d0*(stressOld(i,1:4)+stressNew(i,1:4)),
     &            strainInc(i,1:4))
			  stressPower_plastic = dot_product(0.5d0*(stressOld(i,1:4)+stressNew(i,1:4)),
     &            delta_epsilon_p(1:4))
			else
			  stressPower = dot_product(0.5d0*(stressOld(i,1:6)+stressNew(i,1:6)),
     &            strainInc(i,1:6))
			  stressPower_plastic = dot_product(0.5d0*(stressOld(i,1:6)+stressNew(i,1:6)),
     &            delta_epsilon_p(1:6))
			end if
			enerInternNew(i) = enerInternOld(i) + stressPower / density(i)
			enerInelasNew(i) = enerInelasOld(i) + stressPower_plastic / density(i)

		  end do
		end if
		
		write(*,*) ' stateOld(1,1)=', stateOld(1,1), ' stateNew(1,1)=', stateNew(1,1)
		
        return

      end subroutine vumat