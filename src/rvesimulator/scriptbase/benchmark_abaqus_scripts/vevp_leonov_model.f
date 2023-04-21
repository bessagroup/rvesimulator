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

C    
      REAL(8), dimension(3,3) :: SOId 
      ! second oer identity matrix 
      REAL(8), dimension(3,3) :: etristrain 
      ! trial strain 
      REAL(8), dimension(3,3) :: devstrain, devstrain0 
      !deviatoric strain 
      REAL(8), dimension(3,3) :: evolstrain, evolstrain0
      ! elastic volumetric strain 
      REAL(8), dimension(3,3) :: estrain, estrain0
      ! elastic strain 
      REAL(8), dimension(3,3) :: strain, strain0
      ! elastic strain of previous step 
      REAL(8), dimension(3,3) :: edevstrain, edevstrain0
      ! elastic deviatoric strain 
      ! REAL(8), dimension(3,3) :: etridevstrain 
      ! ! elastic trial deviatoric strain 
      REAL(8), dimension(3,3) :: pdevstrain,pdevstrain0 
      ! plastic deviatotic strain ( previous step)
      REAL(8), dimension(3,3) :: ktotalstress 
      ! total kirchoff stress 
      REAL(8), dimension(3,3) :: kdevstress, kdevstressvisc0, kdevstressvisc
      ! deviatoric kirchoff stress 
      REAL(8), dimension(3,3) :: ktridevstress 
      ! trial deviatoric kirchoff stress 
      REAL(8), dimension(3,3) :: ktristress 
      ! trial kirchoff stress 
      REAL(8), dimension(3,3) :: khstress 
      ! hardening stress 
      REAL(8), dimension(3,3) :: deltastrain 
      ! the delta strain to calculate the equivalent strain rate 
      REAL(8), dimension(3,3) :: kstress, kstress0
      ! kirchoff driving stress 
      REAL(8), dimension(3,3) :: cauchystress 
      ! cauchy stress 
      REAL(8), dimension(3,3) :: kdevstresstotalvisc0 
      ! visco elastic deviatoric stress for all relaxation modes 
      REAL(8), dimension(3,3) :: kstresselast0 
      ! total elastic stress from elastic spring 
      REAL(8), dimension(3,3) :: deltadevstrain
      ! delta deriatoric strain 
      REAL(8), dimension(3,3) :: dC1, dC2, dC3, dC4, dresidual, dC4hard
      ! Derivative of residual term w.r.t total strain
      REAL(8), dimension(3,3) :: dviscosity
      ! Derivative of viscosity w.r.t total strain 
      REAL(8), dimension(3,3) :: deqstrainrate
      ! Derivative of equivalent strain rate w.r.t total strain
      REAL(8), dimension(3,3) :: dgmodvisc
      ! Derivative of visco-elastic shear modulus w.r.t total strain
      REAL(8), dimension(3,3) :: dktridevstressnorm
      ! drivative norm for kichoff trial stress 
C     
      ! REAL(8), dimension(3,3,3,3) :: detensor 
      ! ! elastic tangent operator 
      REAL(8), dimension(3,3,3,3) :: dvevptensor
      ! Viscoelastic-viscoplastic consistent tangent operator 
      REAL(8), dimension(3,3,3,3) :: dktridevstress 
      ! Derivative of Kirchhoff deviatoric trial stress w.r.t total strain     
      REAL(8), dimension(3,3,3,3) :: dkstress 
      ! derivative of kichoff driving stress wrt total strain 
      REAL(8), dimension(3,3,3,3) :: dkhstress 
      ! derivative of kichoff hardening stress wrt total strain 
      REAL(8), dimension(3,3,3,3) :: FOId 
      ! Fourth-order identity tensor
      REAL(8), dimension(3,3,3,3) :: FOtrans 
      ! fourth-order transpose tensor 
      REAL(8), dimension(3,3,3,3) :: FOdiagtrace 
      ! fourth-order 'diagonal trace' tensor 
      REAL(8), dimension(3,3,3,3) :: FOsym 
      ! fourth-order sysmetric projection tensor 
      REAL(8), dimension(3,3,3,3) :: FOdevprojsym 
      ! fourth-order  tensor 
      REAL(8), dimension(3,3) :: dgeq
      !! Derivative of visco-elastic relaxation modes relaxation times w.r.t total strain
      REAL(8), dimension(3,3) :: dw1
      !! Derivative of visco-elastic coefficient 1 w.r.t total strain  
      REAL(8), dimension(3,3,3,3) :: dktridevstress_temp
      !! used for as temperary output for tensor product 
      REAL(8), dimension(3,3,3,3) :: dg_temp 
      ! dg times deltadevstrain
      REAL(8), dimension(3,3,3,3) :: dktristress_temp
      ! dktristress_temp 
      REAL(8), dimension(3,3,3,3) :: daccplasticstrain_temp
C     
      REAL(8), dimension(NTENS) :: pdevstrain0vgt, pdevstrainvgt 
      ! plastic deviotoric strain in voigt notation
      REAL(8), dimension(NTENS) :: estrain0vgt, estrainvgt 
      ! elastic strain in voigt notation
      REAL(8), dimension(NTENS) :: kstress0vgt, kstressvgt , khstressvgt
      ! kirchoff stress in vgt format 
      REAL(8), dimension(NTENS,NTENS) :: dmat 
      ! the DDSDDE 
      REAL(8), dimension(NTENS) :: etristrainvgt
      ! elastic trial strain in the vgt format 
      REAL(8), dimension(NTENS) :: kdevstressviscmodesvgt
      ! all deviatoric kirchoff visco-elastic stress 
      REAL(8), dimension(NTENS) :: kdevstressvisc0vgt
      ! deriatoric stress of visco-elastic 
      REAL(8), dimension(NTENS) :: deltastrainvgt
      ! delta strain vgt 
      REAL(8), dimension(NTENS) :: kdevstressviscvgt
      ! kirchoff deviatoric stress vico-elastic vgt
      REAL(8), dimension(NTENS) :: strainvgt, strain0vgt,kstresselast0vgt
C     
      REAL(8) :: etrivolstrain 
      ! trial volumetric strain
      REAL(8) :: ktridevstressnorm 
      ! trial deviotoric kirchoff stress norm
      REAL(8) :: visco, visco0
      ! viscocity 
      REAL(8) :: residual, residualder
      ! residual and residual derivative   
      REAL(8) :: trihydpress 
      ! trial hydrostatic pressure
      REAL(8) :: accplaststrain, accplaststrain0
      ! equivalent plastic strain
      REAL(8) :: convrgnorm 
      ! norm for convergence check 
      REAL(8) :: jacobian 
      ! jacobian(for transferring kirchoff stress to cauchy stress)
      REAL(8) :: gmod, gmodvisc 
      ! shear modulus of the elastic spring 
      REAL(8) :: kmodvisc 
      ! bulk modulus of the elastic spring 
      REAL(8) :: eqstrainrate 
      ! equivalent strain rate 
      REAL(8) :: geq 
      ! relaxation time of visco-elastic relaxation mode 1 
      REAL(8) :: w1, w2
      ! visco-elastic coefficients 
      REAL(8) :: deltavolstrain 
      ! volume delta strain  
      REAL(8) :: C1, C2, C3, C4 
      ! components of viscosity function 
      REAL(8) :: C1der, C2der, C3der, C4der       
      ! derivative components for viscosity function 
      REAL(8), parameter :: Rgas = 8314.4598484848D0 
      ! constant of gas 
C   
      REAL(8) :: youngs 
      ! youngs modulus 
      REAL(8) :: poisson 
      ! poisson ratio 
      REAL(8) :: T  
      ! temperature  
      REAL(8) :: gt 
      ! temperature-dependent relaxation time of visco-elastic mode 1 
      REAL(8) :: lamda 
      ! visco-elastic relaxation time sensitivity 
      ! REAL(8) :: gviscratio 
      ! ! visco-elastic shear moduli ratio 
      ! REAL(8) :: relaxtimeratio 
      ! ! visco-elastic relaxation time ratio 
      REAL(8) :: G1
      ! shear modulus of visco-elastic relaxation mode 1 
      REAL(8) :: deltaH 
      ! enthalpy 
      REAL(8) :: Astar 
      ! A star 
      REAL(8) :: taustar
      ! tau star 
      REAL(8) :: vstar 
      ! v*
      REAL(8) :: mu  
      ! characteristic stress 
      REAL(8) :: hsoft 
      ! slope of softening 
      REAL(8) :: Dinf 
      ! saturation volume of softening parameter
      REAL(8) :: H 
      ! slope of plasticity curve 
      REAL(8) :: suppress 
      ! atmospheric pressure   
C     
      INTEGER iter,i,j
      REAL*8 R0, R1, R2, R3,small,tolerance,maxiter
      parameter(R0=0.0D0, R1 = 1.0D0, R2=2.0D0, R3=3.0D0)
      parameter(small = 1.0D-8,tolerance =1.0D-6,maxiter = 50)		   
      !                                                    Get arguments 
      !================================================================= 
      youngs = props(1)       ! youngs modulus 
      poisson = props(2)      ! poisson ratio 
      T = props(3)            ! temperature 
      suppress = props(4)     ! atmospheric pressure  
      gt = props(5)           ! gt, for calculating g1
      lamda = props(6)        ! lamda  
      G1 = props(7)           ! G1
      deltaH = props(8)       ! enthalpy
      Astar = props(9)        ! A0
      vstar = props(10)       !  v* shear volume 
      mu = props(11)          ! pressure co-efficient
      hsoft = props(12)       ! slope of softening
      Dinf = props(13)        ! saturation volume of softening parameter
      H = props(14)           ! slope of plasticity curve

      !                                         Define needed parameters
      !================================================================= 
      ! second order identity matrix
      SOId = R0
      forall(i=1:3) SOId(i,i) = R1
      ! fourth-order identity tensor 
      FOId = R0
      forall(i=1:3,j=1:3) FOId(i,j,i,j) = R1 
      ! tranpose of fourth-order identity tensor 
      FOtrans = R0
      forall(i=1:3,j=1:3) FOtrans(i,j,j,i) = R1
      ! fourth-order systemtric projection tensor 
      FOsym = (R1/R2)*(FOId+FOtrans)
      ! Set fourth-order 'diagonal trace' tensor
      call TensorProduct(SOId,SOId,FOdiagtrace)      
      ! Set fourth-order deviatoric projection tensor
      FOdevprojsym = FOsym - (R1/R3) * FOdiagtrace
      ! calculate the increment of time 
      dtime = abs(DTIME)    
      ! calculate the jacobian 
      jacobian = det(DFGRD1) 
      !calculate taustar 
      taustar = (Rgas * T) / vstar 
      ! calculate the shear and bulk modulus  
      gmod = youngs /(R2 *(R1 + poisson)) 
      kmodvisc = youngs /(R3 *(R1 - R2*poisson))

      !                    Get state variables from previous update step 
      !=================================================================      
      ! get values of previous step (vgt format) 
      estrain0vgt = STATEV(1:NTENS) 
      strain0vgt = STATEV(NTENS+1:2*NTENS) 
      kstress0vgt = STATEV(2*NTENS+1:3*NTENS)
      pdevstrain0vgt = STATEV(3*NTENS+1:4*NTENS) 
      accplaststrain0 = STATEV(4*NTENS+1) 
      visco0 = STATEV(4*NTENS+2) 
      kdevstressvisc0vgt = STATEV(4*NTENS+3:5*NTENS+2) 
      ! get the tensor format of those variable 
      call VoigtToTensorStrain(NTENS, estrain0vgt, estrain0)
      call VoigtToTensorStrain(NTENS, pdevstrain0vgt, pdevstrain0) 
      call VoigtToTensorStress(NTENS, kstress0vgt, kstress0) 
      call VoigtToTensorStrain(NTENS, strain0vgt, strain0)

      !                                           Equivalent strain rate
      !================================================================= 
      ! calculate the elastic trial strain 
      call ElasticTrailStrain(NTENS, STRAN, DSTRAN, etristrain)     
      etristrain = etristrain - pdevstrain0
      ! get the delatstrain 
      deltastrain = etristrain-estrain0
      ! calculate the equivalent strain rate 
      eqstrainrate = (R1/dtime)*dsqrt(R2/R3)*norm2(deltastrain) 

      !                                              Material properties
      !=================================================================  
      ! calculate the geq(g1)
      if (eqstrainrate == R0) then
        geq = gt*R1**(lamda)
      else 
        geq = gt*eqstrainrate**(-lamda)
      endif
      ! initialize visco-elastic shear modulus 
      w1 = dexp(-dtime/geq)
      w2 = (R1-w1)/(dtime/geq)
      ! update visco-elastic shear modulus 
      ! (equation 6.97)  
      gmodvisc = gmod  + w2*G1

      !                             Real and algorithmic state variables
      !=================================================================       
      ! Initialize the previous dev kirchoff visco-elastic stress
      call VoigtToTensorStress(NTENS,kdevstressvisc0vgt, kdevstressvisc0)
      kstresselast0 = kstress0 - kdevstressvisc0
      kdevstresstotalvisc0 = w1*kdevstressvisc0

      !                                              Elastic trial state
      !================================================================= 
      ! calculate volumetric incremental strain 
      deltavolstrain= deltastrain(1,1)+ deltastrain(2,2)+ deltastrain(3,3) 
      ! calculate deviatoric incremental strain 
      deltadevstrain = deltastrain-(R1/R3)*deltavolstrain*SOId          
      ! calculate kirchoff trial stress 
      ktristress = (kstresselast0 + kdevstresstotalvisc0) +  
     +         kmodvisc*deltavolstrain*SOId + R2*gmodvisc*deltadevstrain
      ! calculate trial hydrostatic pressure 
      trihydpress = -(R1/R3)*(ktristress(1,1)+ktristress(2,2) + ktristress(3,3))  
      ! Compute Kirchhoff deviatoric trial stress
      ktridevstress = ktristress + trihydpress*SOId       
      ktridevstressnorm = norm2(ktridevstress) 
  
      !                              Viscosity initial (iterative) guess
      !================================================================= 
      if (ktridevstressnorm == R0) then
        DDSDDE = R0
        do i = 1, NDI 
            do j = 1, NDI
                DDSDDE(i, j) = kmodvisc
            enddo
            DDSDDE(i,i) = kmodvisc + R2*gmodvisc 
        enddo
        do i = NDI+1, NTENS
            DDSDDE(i, i) = gmodvisc
        enddo
        ! update the stress 
        do i = 1, NTENS
            do j = 1, NTENS
                stress(i) = stress(i) + ddsdde(i,j) * dstran(j)
            enddo 
        enddo  
        ! update state variables (vgt format) 
        STATEV(1:NTENS)= estrain0vgt
        STATEV(NTENS+1:2*NTENS)= strain0vgt
        STATEV(2*NTENS+1:3*NTENS) = kstress0vgt 
        STATEV(3*NTENS+1:4*NTENS) = pdevstrain0vgt 
        STATEV(4*NTENS+1) = accplaststrain0 
        STATEV(4*NTENS+2)  = visco0
        STATEV(4*NTENS+3:5*NTENS+2) = kdevstressvisc0vgt
      else 
        if (abs(accplaststrain0) < small) then 
          C4 = R0 
          C1 = Astar * dexp((deltaH/(Rgas*T)) +
     +       (mu*(suppress + trihydpress)/taustar)
     +       - Dinf + Dinf*dexp(-hsoft*dsqrt(R3/R2)*C4/Dinf))
          C2 = dsqrt(R1/R2) * ktridevstressnorm
          C3 = dsinh(C2/taustar)        
          visco = C1*(C2/C3)    
        else 
          visco = visco0
        endif 
        ! the main loop of return mapping 
        do iter = 1, maxiter 
          C4 = accplaststrain0
     +         + (R1/R3)*dsqrt(R3/R2)*(dtime/(visco+dtime*gmodvisc))
     +         * ktridevstressnorm
          C1 = Astar * dexp((deltaH/(Rgas*T)) +
     +       (mu*(suppress + trihydpress)/taustar)
     +       -Dinf + Dinf*dexp(-hsoft*dsqrt(R3/R2)*C4/Dinf))     
          C2 = (visco/(visco+dtime*gmodvisc))
     +      * dsqrt(R1/R2) * ktridevstressnorm
          C3 = dsinh(C2/taustar)           
          ! calculate the residual 
          residual = visco - C1*(C2/C3)             
          ! update the convergence norm value 
          if (abs(visco - R0) < tolerance) then
            convrgnorm = dabs(residual)
          else
            convrgnorm = dabs(residual/visco)
          endif
          ! check the convergence cirterion  
          if (convrgnorm < tolerance) then
            exit
          endif
          ! calculate the value of deriative of Residual 
          C4der = -(dtime/(visco+dtime*gmodvisc)**2)*(R1/R3)
     +          *dsqrt(R3/R2)*ktridevstressnorm 
          C1der = -C1*dsqrt(R3/R2)*hsoft*C4der
     +          *dexp(-hsoft*dsqrt(R3/R2)*C4/Dinf)
          C2der = ((dtime*gmodvisc)/(visco+dtime*gmodvisc)**2) 
     +           * dsqrt(R1/R2) * ktridevstressnorm     
          C3der = dcosh(C2/taustar) * (C2der/taustar)
          residualder = R1-(C1/C3)*C2der-(C2/C3)*C1der+(C1*C2/C3**2)*C3der
          if (residualder == R0) then
            residualder = R1
          endif	 
          visco = visco - residual/residualder      
        enddo 
        
        ! compute the deviatoric Kirchhoff stress tensor 
        kdevstress = (visco/(visco+dtime*gmodvisc))*ktridevstress
        ! update the logrithmic plastic deviatoric strain tensor   
        pdevstrain = pdevstrain0 + 
     +             (R1/R2)*(dtime/(visco+dtime*gmodvisc)) * ktridevstress
        ! update the accumulated plastic strain 
        accplaststrain = accplaststrain0+(R1/R3)*dsqrt(R3/R2)
     +            *(dtime/(visco+dtime*gmodvisc))*ktridevstressnorm 
        ! update elastic strain 
        estrain = etristrain - pdevstrain + pdevstrain0
        ! Compute logarithmic strain tensor
        strain = estrain + pdevstrain
        ! Compute logarithmic deviatoric strain tensor
        devstrain = strain - 
     +        (R1/R3)*(strain(1,1) + strain(2,2) + strain(3,3))*SOId 
        devstrain0 = strain0 - 
     +        (R1/R3)*(strain0(1,1) + strain0(2,2) + strain0(3,3))*SOId 
        ! kstress 
        kstress  = kdevstress - trihydpress*SOId  
        ! update kirchoff driving strain 
        evolstrain = estrain(1,1) + estrain(2,2) + estrain(3,3) 
        edevstrain = estrain - (R1/R3)*evolstrain*SOId 
        ! calculate  previous converged logarithmic elastic deviatoric strain
        evolstrain0 = estrain0(1,1) + estrain0(2,2) + estrain0(3,3)
        edevstrain0 = estrain0 - (R1/R3)*evolstrain0*SOId
        ! Compute Kirchhoff driving (deviatoric) visco-elastic stress tensor
        kdevstressvisc = w1*kdevstressvisc0 + 
     +       w2*R2*G1*(edevstrain - edevstrain0)
        ! Store Kirchhoff driving (deviatoric) visco-elastic stress tensor
        call TensorToVoigtStress (NTENS, kdevstressvisc, kdevstressviscvgt)             
        ! Compute Kirchhoff hardening stress tensor accplaststrain
        khstress = H*accplaststrain*devstrain 
        ! Compute Kirchhoff total stress tensor
        ktotalstress = kstress + khstress

        ! update the cauchy stress as the output for abaqus 
        cauchystress = (R1/jacobian)*ktotalstress 
        call TensorToVoigtStress(NTENS, cauchystress, STRESS)
        ! convert variables to vgt format 
        call TensorToVoigtStrain(NTENS, estrain, estrainvgt)
        call TensorToVoigtStrain(NTENS, pdevstrain, pdevstrainvgt) 
        call TensorToVoigtStress(NTENS, kstress, kstressvgt)
        call TensorToVoigtStrain(NTENS, strain, strainvgt )
        call TensorToVoigtStress(NTENS, khstress, khstressvgt)
        call TensorToVoigtStress(NTENS, kstresselast0, kstresselast0vgt)
        ! update state variables (vgt format) 
        ! note the last three state variables are used to track the 
        ! update of different stress component of the devloped model 
        STATEV(1:NTENS)= estrainvgt
        STATEV(NTENS+1:2*NTENS)= strainvgt
        STATEV(2*NTENS+1:3*NTENS) = kstressvgt 
        STATEV(3*NTENS+1:4*NTENS) = pdevstrainvgt 
        STATEV(4*NTENS+1) = accplaststrain 
        STATEV(4*NTENS+2)  = visco 
        STATEV(4*NTENS+3:5*NTENS+2) = kdevstressviscvgt
        STATEV(5*NTENS+3:6*NTENS+2) = khstressvgt
        STATEV(6*NTENS+3:7*NTENS+2) = kstresselast0vgt 
        
  !=====================consistent tangent DMat matrix==================
  ! this is very import for abaqus subroutine, which basically is used  
  ! to caluclate the DDSDDE
  !=====================================================================
        ! calculate the value of deriative of Residual 
        ! update the ktridevstress and its norm 
        ktridevstress = ((visco + dtime*gmodvisc)/visco)*kdevstress
        ktridevstressnorm = norm2(ktridevstress) 
        C4der = -(dtime/(visco+dtime*gmodvisc)**2)*(R1/R3)
     +          *dsqrt(R3/R2)*ktridevstressnorm
        C1der = -C1*dsqrt(R3/R2)*hsoft*C4der
     +          *dexp(-hsoft*dsqrt(R3/R2)*C4/Dinf)
        C2der = ((dtime*gmodvisc)/(visco+dtime*gmodvisc)**2) 
     +           * dsqrt(R1/R2) * ktridevstressnorm     
        C3der = dcosh(C2/taustar) * (C2der/taustar)
        residualder = R1-(C1/C3)*C2der-(C2/C3)*C1der+(C1*C2/C3**2)*C3der
        if (residualder == R0) then
          residualder = R1
        endif	  
   !                     Elasto-viscoplastic consistent tangent operator
   !====================================================================        
        ! calculate the derivative of equivalent strain rate 
        if (eqstrainrate == R0) then 
          dvevptensor = R2*gmodvisc*(FOsym-(R1/R3)*FOdiagtrace) 
     +                 + kmodvisc*FOdiagtrace 
          call FOTensorToMatrix(NTENS, dvevptensor,dmat)      
          DDSDDE = (R1/jacobian)*(dmat)
        else 
          deqstrainrate = (R1/(dtime**2))*(R2/R3)*((strain-strain0)/eqstrainrate)         
          ! loop over the visco-elastic relaxation modes to calculate dgodvics 
          dgeq= -lamda*(geq/eqstrainrate)*deqstrainrate
          dw1 = (dtime/(geq**2))* dexp(-dtime/geq)*dgeq
          dgmodvisc = (G1/dtime)*((R1-dexp(-dtime/geq))*dgeq- geq*dw1)          
          ! loop over visco-elastic relaxation modes to calculate derivatives 
          ! of trial kirchoff stress 
          call TensorProduct((devstrain-devstrain0),dgmodvisc,dg_temp) 
          dktridevstress = R2*gmodvisc*FOdevprojsym + R2*dg_temp   
          call TensorProduct(kdevstressvisc, dw1, dktridevstress_temp)
          dktridevstress = dktridevstress + dktridevstress_temp
          ! Compute derivative of Kirchhoff deviatoric trial stress tensor
          call doublecontraction(ktridevstress/ktridevstressnorm,dktridevstress, 
     +    dktridevstressnorm)      
          ! Compute residual terms derivatives w.r.t total strain
          dC2 = dsqrt(R1/R2)*((visco/(visco + gmodvisc*dtime))*dktridevstressnorm 
     +        - ktridevstressnorm*(dtime /(visco + gmodvisc*dtime)**2)*visco*dgmodvisc)
          dC3 = dcosh((R1/taustar)*C2)*(R1/taustar)*dC2  
          dC4 = (R1/R3)*dsqrt(R3/R2)*((dtime/(visco+gmodvisc*dtime))*dktridevstressnorm 
     +       -((dtime/(visco+gmodvisc*dtime))**2)*ktridevstressnorm*dgmodvisc)
          dC1 = C1*(-(mu/taustar)*kmodvisc*SOId - 
     +    dsqrt(R3/R2)*hsoft*dexp(-dsqrt(R3/R2)*(hsoft/Dinf)*C4)*dC4)
          ! Compute derivative of residual w.r.t total strain
          dresidual = -(C2/C3)*dC1 - (C1/C3)*dC2 + ((C1*C2)/C3**2)*dC3        
          ! Compute derivative of viscosity w.r.t total strain
          dviscosity = -(R1/residualder)*dresidual       
          ! Compute derivative of Kirchhoff driving stress tensor w.r.t total strain
          call TensorProduct(ktridevstress, 
     +   (dtime/(visco + gmodvisc*dtime)**2)*(gmodvisc*dviscosity -visco*dgmodvisc), 
     +   dktristress_temp)    
          dkstress = kmodvisc*FOdiagtrace + dktristress_temp +  
     +            (visco/(visco+gmodvisc*dtime))*dktridevstress
          ! Compute residual terms derivative w.r.t total strain
          dC4hard = (R1/R3)*dsqrt(R3/R2)*((dtime/(visco + 
     +   gmodvisc*dtime))*dktridevstressnorm - (dtime/((visco + 
     +   gmodvisc*dtime))**2)*ktridevstressnorm*(dgmodvisc*dtime+dviscosity))
          ! Compute derivative of Kirchhoff hardening stress tensor w.r.t total strain
          call TensorProduct(devstrain, dC4hard, daccplasticstrain_temp)
          dkhstress = H*(daccplasticstrain_temp + accplaststrain*FOdevprojsym)
          ! Compute viscoelastic-viscoplastic consistent tangent operator
          dvevptensor = dkstress + dkhstress 
          call FOTensorToMatrix(NTENS, dvevptensor, dmat)      
          DDSDDE = (R1/jacobian)*(dmat)
        endif 
      endif

	  return
	  end 

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
C doublecontraction: tensor multiply (3,3) and (3,3,3,3) matrix
C input:  vector(6)
C output: tensor(3x3)
C ----------------------------------------------------------------------
      subroutine doublecontraction(a, b, c)
        implicit none 
        !!## Summary
        !! Double contraction between a 2nd order tensor and a 4th 
        !! order tensor
        !! (kl = ij:ijkl).
        ! ARGUMENTS ====================================================
        real(8), dimension(3,3)  :: a
          !! Second-order tensor (ij)
        real(8), dimension(3,3,3,3) :: b
          !! Fourth-order tensor (ijkl)
        real(8), dimension(3,3) :: c
          !! Double contraction
        ! LOCAL VARIABLES ==============================================
        integer :: i, j, k, l
        real(8) :: R0 
        parameter (R0=0.0D0)
          !! Auxiliary variables
        !===============================================================
        !                                      doublecontraction
        ! Compute double contraction
        c = R0
        do k = 1, size(b,3)
          do l = 1, size(b,4)
            do i = 1, size(a,1)
              do j = 1, size(a,2)
                c(k,l) = c(k,l) + a(i,j)*b(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

        return
        !                                      doublecontraction_ij_ijkl
        !===============================================================
      end subroutine
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
        else 
          tensor(1,1) = vector(1)
          tensor(2,2) = vector(2)
          tensor(3,3) = vector(3)
          tensor(1,2) = vector(4)
          tensor(1,3) = vector(5)
          tensor(2,3) = vector(6)
        endif 
        
        !complete the systemtric tensor  
        tensor(2,1) = tensor(1,2)
        tensor(3,1) = tensor(1,3)
        tensor(3,2) = tensor(2,3)

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
      end subroutine
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
      end subroutine

