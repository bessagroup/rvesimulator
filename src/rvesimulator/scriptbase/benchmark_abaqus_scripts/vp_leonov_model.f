      subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     +                RPL,DDSDDT,DRPLDE,DRPLDT,
     +                STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     +                CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     +                DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     +                KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      character*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     +          DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     +          TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     +          DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)	  

      REAL(8), dimension(3,3) :: SOld ! identity matrix
      REAL(8), dimension(3,3) :: eTriStrain ! trial strain
	    REAL(8), dimension(3,3) :: devStrain ! deviotoric strain 
      REAL(8), dimension(3,3) :: eStrain, eStrainO ! elastic strain matrix
      REAL(8), dimension(3,3) :: eDevStrain ! elastic deviotoric strain
      REAL(8), dimension(3,3) :: eTrialDevStrain ! trial deviotoric strain
      REAL(8), dimension(3,3) :: pDevStrainO, pDevStrain! plastic deviotoric strain 
      REAL(8), dimension(3,3) :: kTotalStress ! total stress (kirchoff driving stress + hardening stress)
      REAL(8), dimension(3,3) :: kDevStress ! deviotoric kirchoff driving stress
      REAL(8), dimension(3,3) :: kTrialDevStress ! trial kirchoff stress 
      REAL(8), dimension(3,3) :: kHStress ! hardening stress
      REAL(8), dimension(3,3) :: CauchyStress
      REAL(8), dimension(3,3) :: kStress ! kirchoff driving stress tensor
c
      REAL(8), dimension(3,3,3,3) :: FO1, FO2, FO3, FO4, FO5
      REAL(8), dimension(3,3,3,3) :: FOld, FTran
      REAL(8), dimension(3,3,3,3) :: FOAux1
      REAL(8), dimension(3,3,3,3) :: FOAux2
      REAL(8), dimension(3,3,3,3) :: FODiagTrace
      REAL(8), dimension(3,3,3,3) :: FOSym
      REAL(8), dimension(3,3,3,3) :: FODevProjSym
	    REAL(8), dimension(3,3,3,3) :: DEVPTensor	 
c    
      ! some Vogit notations of variables  
      REAL(8), dimension(NTENS) :: pDevStrainOVgt, pDevStrainVgt ! plastic deviotoric strain in voigt notation
      REAL(8), dimension(NTENS) :: eStrainOVgt, eStrainVgt ! elastic strain in voigt notation
      REAL(8), dimension(NTENS,NTENS) :: DMat ! the DDSDDE 
c
      REAL(8) :: eTrialVolStrain ! trial volumetric strain
      REAL(8) :: kTrialDevStressNorm ! trial deviotoric kirchoff stress norm
      REAL(8) :: visco, viscoO
      REAL(8) :: residual, residualDer! residual 
      REAL(8) :: C1, C2, C3, C4 
      REAL(8) :: C1Der, C2Der, C3Der, C4Der
      REAL(8) :: F1, F2, F3 
      REAL(8) :: trialHydPress ! hydrostatic pressure
      REAL(8) :: accPlastStrain, accPlastStrainO ! equivalent plastic strain
      REAL(8) :: convrgNorm
      REAL(8) :: Je ! jacobian
C 
      ! Arguments of this subroutine, which are some important parameters 
      ! to govern the Mirkhalaf model 
      REAL(8) :: E ! youngs modulus
      REAL(8) :: niu ! poissons ratio
      REAL(8) :: T ! temparature
      REAL(8) :: h_soft ! parameter influencing the softening behaviour
      REAL(8) :: D_inf ! saturation volume
      REAL(8) :: mue ! ratio of pressure activation volume and shear activation volume
      REAL(8) :: supPress ! atmospheric pressure
      REAL(8) :: H ! hardening modulus
      REAL(8) :: delta_H ! activation energy
      REAL(8) :: tau_star ! charecteristic stress
      REAL(8) :: v_star ! shear volume 
      REAL(8) :: A_star ! factor associated with the vibration
      REAL(8) ::  R_gas ! gas constant
      REAL(8) :: G ! shear modulus
      REAL(8) :: K ! bulk modulus

      INTEGER iter,i,j
      REAL*8 R0, R1, R2, R3,small,tolerance,maxiter
      parameter(R0=0.0D0, R1 =1.0D0, R2=2.0D0, R3=3.0D0)	 
      parameter(small = 1.0D-8,tolerance =1.0D-6,maxiter = 100)		
      parameter(R_gas = 8314.4598484848D0)
      E = props(1)           ! youngs modulus (MPa)
      niu = props(2)         ! poissons ratio (None)
      delta_H = props(3)     ! enthalpy (mJ/mol) 
      A_star = props(4)      ! A0 (s)
      v_star = props(5)       ! shear volume (mm3/mol) 
      mue = props(6)         ! pressure co-efficient (None)
      h_soft = props(7)      ! slope of softening (None)
      D_inf = props(8)       ! saturation volume of softening parameter (None)
      H = props(9)           ! slope of plasticity curve (None)
      T = props(10)          ! temparature (K)
      supPress = props(11)    ! atmospheric pressure   

      ! define the G and K, for elastic property 
      G = E/((R1 + niu)*R2) ! shear Modulus
      K = E/((R1-R2*niu)*R3) ! Bulk Modulus
      tau_star =  (R_gas*T) / v_star   ! charectoristic stress  

      ! initialized importan parameters 
      ! second order tensor SOld 
      SOld = R0
      forall(i =1:3) SOld(i,i) = R1
      !fourth-order identity tensor  
      FOld = R0 
      forall(i=1:3, j=1:3) FOld(i,j,i,j) = R1       
      ! transpose of fourth-order identity tensor 
      FTran = R0  
      forall(i=1:3,j=1:3) FTran(i,j,j,i) = R1
      ! Fourth-order systemtri projecion tensor 
      FOSym = (R1/R2)*(FOld+FTran)
      call TensorProduct(SOld,SOld,FODiagTrace)
      FODevProjSym = FOSym - (R1/R3) * FODiagTrace
      ! determint of deformation gradien
      Je = det(DFGRD1) 
      DTIME = dabs(DTIME)
	  ! initilize the variables for calculation 
      eStrainOVgt = STATEV(1:NTENS) 
      accPlastStrainO = STATEV(NTENS+1) 
      pDevStrainOVgt = STATEV(NTENS+2:NTENS+NTENS+1) 
      viscoO = STATEV(NTENS+NTENS+2)
      call VoigtToTensorStrain(NTENS, eStrainOVgt, eStrainO)
      call VoigtToTensorStrain(NTENS, pDevStrainOVgt, pDevStrainO) 
      call ElasticTrailStrain(NTENS, STRAN, DSTRAN, eTriStrain)
      eTriStrain = eTriStrain - pDevStrainO
c 	  
    ! calculate volumetric trial strain 
      eTrialVolStrain = eTriStrain(1,1)+eTriStrain(2,2) + eTriStrain(3,3)
      trialHydPress = -K*eTrialVolStrain
      eTrialDevStrain = eTriStrain -(R1/R3) *eTrialVolStrain * SOld
      kTrialDevStress = R2 * G * eTrialDevStrain		  	  	  
      kTrialDevStressNorm = norm2(kTrialDevStress)		  	   
	  if (kTrialDevStressNorm == R0) then
        DDSDDE = R0
        do i = 1, NDI 
            do j = 1, NDI
                DDSDDE(i, j) = K
            enddo
            DDSDDE(i,i) = K + R2*G 
        enddo
        do i = NDI+1, NTENS
            DDSDDE(i, i) = G
        enddo
        ! update the stress 
        do i = 1, NTENS
            do j = 1, NTENS
                stress(i) = stress(i) + DDSDDE(i,j) * dstran(j)
            end do 
        end do  
        ! update statev variables             
        STATEV(1:NTENS) = eStrainOVgt
        STATEV(NTENS+1) = accPlastStrainO
        STATEV(NTENS+2:NTENS+NTENS+1) = pDevStrainOVgt
        STATEV(NTENS+NTENS+2) = viscoO  
	  else
	      if (accPlastStrainO < small) then
		      C4 = R0
			    C1 = A_star * dexp((delta_H/(R_gas*T))
     +    +(mue*(supPress + trialHydPress)/tau_star)) 

			    C2 = dsqrt(R1/R2) * kTrialDevStressNorm
			    C3 = dsinh(C2/tau_star)
			    visco = C1*C2/C3		  
		    else
		      visco = viscoO
		  
		    endif
c     ! the mian loop of return mapping        
		    do iter = 1, maxiter 
          ! write(6,*), "--------------------iter------------------", iter
		      C4 = accPlastStrainO 
     +         + (R1/R3)*dsqrt(R3/R2)*(DTIME/(visco+DTIME*G))
     +         * kTrialDevStressNorm
			    C1 = A_star*dexp((delta_H/(R_gas*T)) 
     +         +(mue * (supPress + trialHydPress)/tau_star)-D_inf
     +         + D_inf*dexp(-h_soft*dsqrt(R3/R2)*C4/D_inf))
			    C2 = (visco/(visco+DTIME*G))
     +      * dsqrt(R1/R2) * kTrialDevStressNorm
			    C3 = dsinh(C2/tau_star)
          ! calculate the residual 
			    residual = visco - C1*(C2/C3)	 
          ! update the convergence norm value 
          if (abs(visco - R0) < tolerance) then
            convrgNorm = dabs(residual)
          else
            convrgNorm = dabs(residual/visco)
          end if
          ! check the convergence cirterion  
          if (convrgNorm < tolerance) then
            exit
          endif
          ! calculate the value of deriative of Residual 
          C4Der = -(DTIME/(visco+DTIME*G)**2)*(R1/R3)
     +            *dsqrt(R3/R2)*kTrialDevStressNorm

          C1Der = -C1*dsqrt(R3/R2)*h_soft*C4Der
     +     *dexp(-h_soft*dsqrt(R3/R2)*C4/D_inf)

          C2Der = (DTIME*G/(visco+DTIME*G)**2) 
     +           * dsqrt(R1/R2) * kTrialDevStressNorm

          C3Der = dcosh(C2/tau_star) * (C2Der/tau_star)
c         ! update the risidual derivative 
			    residualDer = R1-(C1/C3)*C2Der-(C2/C3)*C1Der+(C1*C2/C3**2)*C3Der	  

          if (residualDer == R0) then
            residualDer = R1
          endif	  
c            ! update the visco function 
			    visco = visco - residual/residualDer 	
	      enddo

        ! update the output variables 
        ! compute the deviatoric Kirchhoff stress tensor 
        kDevStress = (visco/(visco+DTIME*G))*kTrialDevStress
        ! update the elastic deviatoric strian of the next step 
		    eDevStrain = kDevStress/(R2*G)
        ! update the logrithmic plastic deviatoric strain tensor   
		    pDevStrain = pDevStrainO + (DTIME/(R2*visco+DTIME*G)) * kTrialDevStress
        ! update the accumulated plastic strain 
		    accPlastStrain =accPlastStrainO+(R1/R3)*dsqrt(R3/R2)
     +    *(DTIME/(visco+DTIME*G))*kTrialDevStressNorm
        ! update total logarithmic strain 
		    devStrain = eDevStrain + pDevStrain
        ! update elastic strain 
        eStrain = eTriStrain - pDevStrain + pDevStrainO 
        ! update the driving force  
		    kStress = kDevStress - trialHydPress * SOld
        ! update the hardening stress 
		    kHStress = H * devStrain
        ! update the total kirchhoff stress  
		    kTotalStress = kStress + kHStress
c
	    ! changing tensor notation to voigt notation
	      call TensorToVoigtStrain(NTENS,eStrain,eStrainVgt)
		    call TensorToVoigtStrain(NTENS,pDevStrain,pDevStrainVgt)
c
	      !  updating the state variables
	      STATEV(1:NTENS) = eStrainVgt
        STATEV(NTENS+1) = accPlastStrain
        STATEV(NTENS+2:NTENS+NTENS+1) = pDevStrainVgt
        STATEV(NTENS+NTENS+2) = visco
        CauchyStress = (R1/Je) * kTotalStress
        call TensorToVoigtStress(NTENS,CauchyStress,STRESS)
c			  			  
!########################### consistent tangent DMat matrix#############
        ! update norm and kdevstress 
        ! kTrialDevStress = ((visco+DTIME*G)/visco)*kDevStress
        ! kTrialDevStressNorm = norm2(kTrialDevStress)
        C4Der = -(DTIME/(visco+DTIME*G)**2)*(R1/R3)
     +            *dsqrt(R3/R2)*kTrialDevStressNorm

        C1Der = -C1*dsqrt(R3/R2)*h_soft*C4Der
     +     *dexp(-h_soft*dsqrt(R3/R2)*C4/D_inf)

        C2Der = (DTIME*G/(visco+DTIME*G)**2) 
     +           * dsqrt(R1/R2) * kTrialDevStressNorm

        C3Der = dcosh(C2/tau_star) * (C2Der/tau_star)
c     ! update the risidual derivative 
        residualDer =R1-(C1/C3)*C2Der-(C2/C3)*C1Der+(C1*C2/C3**2)*C3Der	  
        if (residualDer == R0) then
          residualDer = R1
        endif	  
C     ! begin to calculate the F1, F2 and F3, which are critial   
        F1 = dsqrt(R2)*G*(C1/C3) * (visco/(visco+DTIME*G)) 
     +      * (R1/kTrialDevStressNorm)
        F2 = C2*(R1/dsqrt(R2))*((h_soft*DTIME/visco)) 
     +     *dexp(-dsqrt(R3/R2) * (h_soft/D_inf) * C4)
        F3 = (R1/tau_star) * (C2/C3) * dcosh(C2/tau_star)
c		  		  		  
        call TensorProduct(kDevStress,kDevStress,FOAux1)
        call TensorProduct(kDevStress,SOld,FOAux2)
c
        FO1=(DTIME*G/(visco**2))*(R1/residualDer)*F1*(R1-F2-F3)*FOAux1	  
        FO2 = -(DTIME*G/(visco*(visco+DTIME*G)))*(R1/residualDer)* 
     +          C1*(C2/C3)*(mue*K/tau_star)*FOAux2
        FO3 = (visco/(visco+DTIME*G))*R2*G*FODevProjSym
        FO4 = K * FODiagTrace
        FO5 = H * FODevProjSym
        DEVPTensor = FO1+FO2+FO3+FO4+FO5 
        call FOTensorToMatrix(NTENS, DEVPTensor,DMat)
        DDSDDE = (R1/Je)*(DMat)
	  endif

	  return
	  end subroutine

C ----------------------------------------------------------------------
C  det: Calculate determinant of 3x3 matrix
C  Input: A (3x3 matrix)
C ----------------------------------------------------------------------
      function det(A)
          implicit none
          real*8 det,A
          dimension A(3,3)
          det = A(1,1)*A(2,2)*A(3,3)
     1        - A(1,2)*A(2,1)*A(3,3)
     2        - A(1,1)*A(2,3)*A(3,2)
     3        + A(1,3)*A(2,1)*A(3,2)
     4        + A(1,2)*A(2,3)*A(3,1)
     5        - A(1,3)*A(2,2)*A(3,1)
          return
      end 
C ----------------------------------------------------------------------
C ElasticTrailStrain: tensor multiply 3x3 matrix
C input:  vector(6)
C output: tensor(3x3)
C ----------------------------------------------------------------------
      subroutine ElasticTrailStrain(NTENS, Strain, DStrain,eTrialStrain) 
        implicit none 
        real*8 Strain, DStrain, eTrialStrain, R0, R2
        integer NTENS
        dimension Strain(NTENS), DStrain(NTENS)
        dimension eTrialStrain(3,3)
        parameter (R0=0.0D0, R2 = 2.0D0)

        if (NTENS == 4) then

          eTrialStrain(1,1) = Strain(1) + DStrain(1)
          eTrialStrain(2,2) = Strain(2) + DStrain(2)
          eTrialStrain(3,3) = Strain(3) + DStrain(3)
          eTrialStrain(1,2) = Strain(4) + DStrain(4)/ R2
          eTrialStrain(1,3) = R0
          eTrialStrain(2,3) = R0
        else 
          eTrialStrain(1,1) = Strain(1) + DStrain(1)
          eTrialStrain(2,2) = Strain(2) + DStrain(2)
          eTrialStrain(3,3) = Strain(3) + DStrain(3)
          eTrialStrain(1,2) = Strain(4) + DStrain(4)/ R2
          eTrialStrain(1,3) = Strain(5) + DStrain(5)/ R2
          eTrialStrain(2,3) = Strain(6) + DStrain(6)/ R2

        endif

        ! complete the sysmetric tensor 
        eTrialStrain(2,1) = eTrialStrain(1,2)
        eTrialStrain(3,1) = eTrialStrain(1,3)
        eTrialStrain(3,2) = eTrialStrain(2,3) 	
        ! write(6, *), 'etrial', eTrialStrain

        return 
        end subroutine

C ----------------------------------------------------------------------
C VoigtToTensorStress: tensor multiply 3x3 matrix
C input:  vector(6)
C output: tensor(3x3)
C ----------------------------------------------------------------------
      subroutine VoigtToTensorStress(NTENS,vector,tensor)
        implicit none
        real*8 vector,tensor, R0
        integer NTENS
        dimension tensor(3,3), vector(NTENS) 
        parameter (R0=0.0D0) 

        if (NTENS == 4) then 
          tensor(1,1) = vector(1) 
          tensor(2,2) = vector(2)
          tensor(3,3) = vector(3)
          tensor(1,2) = vector(4)
          tensor(1,3) = R0
          tensor(2,3) = R0
          ! write(6, *), 'stress_tensor', tensor
        else 
          tensor(1,1) = vector(1)
          tensor(2,2) = vector(2)
          tensor(3,3) = vector(3)
          tensor(1,2) = vector(4)
          tensor(1,3) = vector(5)
          tensor(2,3) = vector(6)
          ! write(6, *), 'stress_tensor', tensor
        endif 
        
        !complete the systemtric tensor  
        tensor(2,1) = tensor(1,2)
        tensor(3,1) = tensor(1,3)
        tensor(3,2) = tensor(2,3)
        ! write(6, *), 'stress_tensor', tensor

        return
	    end subroutine
C ----------------------------------------------------------------------
C VoigtToTensorStrain: tensor multiply 3x3 matrix
C input:  vector(6)
C output: tensor(3x3)
C ----------------------------------------------------------------------	   
      subroutine VoigtToTensorStrain(NTENS,vector,tensor)
        implicit none
        real*8 vector,tensor, R0, R2 
        integer NTENS
        dimension tensor(3,3), vector(NTENS) 	  
        parameter (R0 = 0.0D0, R2 = 2.0D0)

        if (NTENS == 4)  then 
          tensor(1,1) = vector(1) 
          tensor(2,2) = vector(2)
          tensor(3,3) = vector(3)
          tensor(1,2) = vector(4)/R2 
          tensor(1,3) = R0
          tensor(2,3) = R0 
        !  write(6, *), 'VoigtToTensorStrain', tensor
        else 
          tensor(1,1) = vector(1)
          tensor(2,2) = vector(2)
          tensor(3,3) = vector(3)
          tensor(1,2) = vector(4)/R2
          tensor(1,3) = vector(5)/R2
          tensor(2,3) = vector(6)/R2
        endif 
        tensor(2,1) = tensor(1,2)
        tensor(3,1) = tensor(1,3)
        tensor(3,2) = tensor(2,3) 
        ! write(6, *), 'VoigtToTensorStrain', tensor

        return
      end subroutine
C ----------------------------------------------------------------------
C TensorToVoigtStress: tensor multiply 3x3 matrix
C input:  tensor(3x3)
C output: vector(6) 
C ----------------------------------------------------------------------	  
      subroutine TensorToVoigtStress(NTENS,tensor,vector)
        implicit none
        real*8 vector,tensor
        integer NTENS
        dimension tensor(3,3), vector(NTENS)	  
        if (NTENS == 4 ) then 
          vector(1) = tensor(1,1)
          vector(2) = tensor(2,2)
          vector(3) = tensor(3,3)
          vector(4) = tensor(1,2)
        else 
          vector(1) = tensor(1,1)
          vector(2) = tensor(2,2)
          vector(3) = tensor(3,3)
          vector(4) = tensor(1,2)
          vector(5) = tensor(1,3)
          vector(6) = tensor(2,3)
        endif 
        ! write(6, *), 'TensorToVoigtStress', vector
        return
      end subroutine	 
C ----------------------------------------------------------------------
C TensorToVoigtStrain: tensor multiply 3x3 matrix
C input:  tensor(3x3)
C output: vector(6) 
C ----------------------------------------------------------------------	
      subroutine TensorToVoigtStrain(NTENS,tensor,vector)
        implicit none
        real*8 vector,tensor, R2
        integer NTENS
        dimension tensor(3,3), vector(NTENS)	
        parameter (R2 = 2.0D0)
        
        if (NTENS == 4) then 
          vector(1) = tensor(1,1)
          vector(2) = tensor(2,2)
          vector(3) = tensor(3,3)
          vector(4) = R2*tensor(1,2) 
        else 
          vector(1) = tensor(1,1)
          vector(2) = tensor(2,2)
          vector(3) = tensor(3,3)
          vector(4) = R2 * tensor(1,2) 
          vector(5) = R2 * tensor(1,3)
          vector(6) = R2 * tensor(2,3)
        endif
        ! write(6, *), 'TensorToVoigtStrain', vector
        return
      end subroutine
C ----------------------------------------------------------------------
C TensorProduct: tensor multiply 3x3 matrix
C input:  tensor1(3x3), tensor(3x3)
C output: tensor(3x3x3x3) 
C ----------------------------------------------------------------------
      subroutine TensorProduct(a,b,c)
        implicit none
        real*8 a,b,c
        dimension a(3,3),b(3,3),c(3,3,3,3)
        integer i,j,k,l
        do i=1,3
            do j=1,3
                do k=1,3
                    do l = 1,3
                      c(i,j,k,l) = a(i,j) * b(k,l)
                    enddo	  
                enddo
            enddo
        enddo
        return
      end
C ----------------------------------------------------------------------
C FOTensorToMatrix: tensor multiply 3x3 matrix
C input:  tensor(3x3x3x3) 
C output: matrix(6x6)
C ----------------------------------------------------------------------	  
      subroutine FOTensorToMatrix(NTENS, FO_tensor,matrix)
        implicit none
        integer i,j,k,l,irow,icol,NTENS
        real*8 FO_tensor,matrix,indexMatrix
        dimension FO_tensor(3,3,3,3)
        dimension matrix(NTENS,NTENS),indexMatrix(2,NTENS)

        if (NTENS == 6) then 
          indexMatrix = reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,NTENS/))

        else 
          indexMatrix = reshape((/1,1,2,2,3,3,1,2/),(/2,NTENS/)) 
        endif 
        do irow = 1,NTENS
            do icol = 1,NTENS
                i = indexMatrix(1,irow)
                j = indexMatrix(2,irow)				  
                k = indexMatrix(1,icol)
                l = indexMatrix(2,icol)
                matrix(irow,icol) = FO_tensor(i,j,k,l)
            enddo
        enddo

        return
      end

