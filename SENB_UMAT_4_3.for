************************************************************************
!     A finite element implementation of a large deformation gradient-
!     damage phase-field theory for fracture using Abaqus user material subroutines
!     Keven Alkhoury, Shawn A. Chester, and Vikas Srivastava
!     Keven Alkhoury, September 2025, Abaqus 2024, In review
************************************************************************
!     State Variables
!     --------------------------------------------------------------
!     statev(1:9) = Fp ---------- Plastic distortion
!     statev(10)  = gBarP ------- Equiv. plastic shear strain
!     statev(11)  = nuP --------- Equiv. plastic shear strain rate
!     statev(12)  = Y  ---------- Isotropic def. resistance
!     statev(13) = Histmax_tau ------ Referential history function
!     statev(14) = detF_tau ------ Determinant of the deformation gradient
!     statev(15) = F_tau(1,1) ------ Deformation gradient (1,1)
!     statev(16) = F_tau(2,2) ------ Deformation gradient (2,2)
!     statev(17) = F_tau(3,3) ------ Deformation gradient (3,3)
!     statev(18) = F_tau(1,2) ------ Deformation gradient (1,2)
!     statev(19) = F_tau(1,3) ------ Deformation gradient (1,3)
!     statev(20) = F_tau(2,1) ------ Deformation gradient (2,1)
!     statev(21) = F_tau(2,3) ------ Deformation gradient (2,3)
!     statev(22) = F_tau(3,1) ------ Deformation gradient (3,1)
!     statev(23) = F_tau(3,2) ------ Deformation gradient (3,2)
!     --------------------------------------------------------------
!     Material Properties Vector
!     --------------------------------------------------------------
!     Mechanical (UMAT):
!     --------------------------------------------------------------
!     Eyoung  = props(1)  ! Elastic modulus
!     poisson = props(2)  ! Poisson ratio
!     Y0      = props(3)  ! Initial resistance
!     H       = props(4)  ! Hardening modulus
!     Ysat    = props(5)  ! Saturation resistance
!     mRate   = props(6)  ! Rate sensitivity
!     nu0     = props(7)  ! Reference shear rate
!     erff = props(8)  -- Fracture energy/unit reference volume
!     lc = props(9) ------- Intrinsic characteristic length scale
!     zeta = props(10) ------ Viscous regularization parameter
!     psi_cr = proprs(11) ---- Critical fracture energy
!     --------------------------------------------------------------
!     Thermal (UMATHT):
!     --------------------------------------------------------------
!      erff = props(1) --- Fracture energy/unit reference volume
!      lc = props(2) ------ Intrinsic characteristic length scale
!      zeta = props(3) -- Viscous regularization parameter
!      specht = props(4) ---- Specific heat
!**********************************************************************

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      include 'aba_param.inc'

      dimension stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character*80 cmname,file1
      character*256 jobName,outDir,fileName

      integer i,j,k,l,iterError,lenJobName,lenOutDir

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),theta_t,theta_tau,T_tau(3,3)
      real*8 MatTan(3,3,3,3),Y_t,Y_tau,Fp_t(3,3),Fp_tau(3,3),gBarP_t
      real*8 gBarP_tau,nuP_t,nuP_tau,gBarLmt,plasticwork
      real*8 StressTempJac(3,3),PlWrkDefJac(3,3),dplasticwork
      real*8 Histmax_t,Histmax_tau,detF_t,detF_tau,rpl


      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0)


      ! Obtain old and new deformation gradients
      !
      F_t = dfgrd0
      F_tau = dfgrd1
      !
      ! Obtain old and new temperature
      !
      theta_t = temp
      theta_tau = temp + dtemp
      !
      ! Obtain state variables from the previous increment
      !
      if(time(2).eq.zero) then
          Fp_t = zero
          Fp_t(1,1) = one
          Fp_t(2,2) = one
          Fp_t(3,3) = one
          gBarP_t = zero
          nuP_t = zero
          Y_t = props(3)
          Histmax_t = zero
          detF_t = one
          F_t(1,1) = one
          F_t(2,2) = one
          F_t(3,3) = one
          F_t(1,2) = zero
          F_t(1,3) = zero
          F_t(2,1) = zero
          F_t(2,3) = zero
          F_t(3,1) = zero
          F_t(3,2) = zero
      else
          Fp_t(1,1) = statev(1)
          Fp_t(2,2) = statev(2)
          Fp_t(3,3) = statev(3)
          Fp_t(2,3) = statev(4)
          Fp_t(3,2) = statev(5)
          Fp_t(1,3) = statev(6)
          Fp_t(3,1) = statev(7)
          Fp_t(1,2) = statev(8)
          Fp_t(2,1) = statev(9)
          gBarP_t   = statev(10)
          nuP_t     = statev(11)
          Y_t       = statev(12)
          Histmax_t = statev(13)
          detF_t = statev(14)
          F_t(1,1) = statev(15)
          F_t(2,2) = statev(16)
          F_t(3,3) = statev(17)
          F_t(1,2) = statev(18)
          F_t(1,3) = statev(19)
          F_t(2,1) = statev(20)
          F_t(2,3) = statev(21)
          F_t(3,1) = statev(22)
          F_t(3,2) = statev(23)
      endif
      !
      ! Perform the time integration implicitly, compute the stress
      !  and the state at the end of the increment
      !
      call integ(props,nprops,dtime,
     +     F_tau,theta_tau,
     +     Fp_t,gBarP_t,nuP_t,Y_t,
     +     Fp_tau,gBarP_tau,nuP_tau,Y_tau,
     +     Histmax_t,Histmax_tau,detF_tau,
     +     T_tau,MatTan,rpl,drpldt)
      !
      ! Update state variables at the end of the increment
      !
      statev(1)  = Fp_tau(1,1)
      statev(2)  = Fp_tau(2,2)
      statev(3)  = Fp_tau(3,3)
      statev(4)  = Fp_tau(2,3)
      statev(5)  = Fp_tau(3,2)
      statev(6)  = Fp_tau(1,3)
      statev(7)  = Fp_tau(3,1)
      statev(8)  = Fp_tau(1,2)
      statev(9)  = Fp_tau(2,1)
      statev(10) = gBarP_tau
      statev(11) = nuP_tau
      statev(12) = Y_tau
      statev(13) = Histmax_tau
      statev(14) = detF_tau
      statev(15) = F_tau(1,1)
      statev(16) = F_tau(2,2)
      statev(17) = F_tau(3,3)
      statev(18) = F_tau(1,2)
      statev(19) = F_tau(1,3)
      statev(20) = F_tau(2,1)
      statev(21) = F_tau(2,3)
      statev(22) = F_tau(3,1)
      statev(23) = F_tau(3,2)
      !
      ! Time stepping algorithim based on the constitutive response
      !
      ! Return Abaqus/Standard the Cauchy stress
      !
      if(ntens.eq.6) then
         !
         ! 3D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
         stress(5) = T_tau(1,3)
         stress(6) = T_tau(2,3)
      elseif(ntens.eq.4) then
         !
         ! 2D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
      endif
      !
      ! Return Abaqus/Standard the stress-deformation jacobian
      !
      if(ntens.eq.6) then
         call jac3D(MatTan,ddsdde)
      elseif(ntens.eq.4) then
         call jac2D(MatTan,ddsdde)
      endif
      !
      ! Return the source term and its derivative (Appendix A.3)
      !
      rpl = rpl
      !
      drpldt = drpldt
      !
      return
      end subroutine umat
         
***********************************************************************
      ! Check Appendix A and B
      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,statev,temp,
     $     dtemp,dtemdx,time,dtime,predef,dpred,cmname,ntgrd,nstatv,
     $     props,nprops,coords,pnewdt,noel,npt,layer,kspt,kstep,kinc)
c
      include 'aba_param.inc'
c
      character*80 cmname
c
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     $     dfdg(ntgrd,ntgrd),statev(nstatv),dtemdx(ntgrd),time(2),
     $     predef(1),dpred(1),props(nprops),coords(3)
      real*8 F(3,3),B(3,3),erff,lc,zeta,aux,detF,B_tau(3,3),detF_tau
      !
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,six=6.d0,eight=8.d0)
c
c
      !
      ! Obtain material properties
      !
      erff = props(1)
      lc = props(2)
      zeta = props(3)
      specht = props(4)
      !
      detF_tau = statev(14)
      !
      ! Used to update conduction term (Appendix A.2)
      !
      aux = (erff*(lc**two))
      !
      ! Remember, these are F_tau as UMAT is called before UMATHT
      !
      F(1,1) = statev(15)
      F(2,2) = statev(16)
      F(3,3) = statev(17)
      F(1,2) = statev(18)
      F(1,3) = statev(19)
      F(2,1) = statev(20)
      F(2,3) = statev(21)
      F(3,1) = statev(22)
      F(3,2) = statev(23)
      !
      B_tau = matmul(F,transpose(F))
      !
      ! Update transient term (Appendix A.1)
      !
      dudt = specht/detF_tau
      du = dudt*dtemp
      u = u+du
      !
      ! Update conduction term (Appendix A.2)
      !
      flux = zero
      do i=1, ntgrd
          do j =1, ntgrd
         flux(i) = flux(i) - (B_tau(i,j)/detF_tau)*dtemdx(j)
          end do
      end do
      flux = flux * aux
      !
      dfdg = zero
      do i=1, ntgrd
          do j=1, ntgrd
          dfdg(i,j) = dfdg(i,j) - (B_tau(i,j)/detF_tau)
          end do
      end do
      dfdg = dfdg * aux
      !
      return
      end
***********************************************************************

***********************************************************************

      subroutine integ(props,nprops,dtime,
     +     F_tau,theta_tau,
     +     Fp_t,gBarP_t,nuP_t,Y_t,
     +     Fp_tau,gBarP_tau,nuP_tau,Y_tau,
     +     Histmax_t,Histmax_tau,detF,
     +     T_tau,MatTan,rpl,drpldt)
      

      implicit none
      
      integer i,j,k,l,m,n,nprops,nargs
      parameter(nargs=10)
      
      real*8 props(nprops),dtime,F_tau(3,3),Fp_t(3,3),Fp_tau(3,3),nuP_t
      real*8 gBarP_t,gBarP_tau,T_tau(3,3),MatTan(3,3,3,3),Y_tau,Y_t,H
      real*8 nuP_tau,T_ta9(3,3),plasticwork,Eyoung,Kbulk,Iden(3,3),mRate
      real*8 poisson,Gshear,Lambda,omega,Y0,Ysat,Fp_t_inv(3,3),theta0
      real*8 Fe_tr(3,3),Re_tr(3,3),Ue_tr(3,3),Ee_tr(3,3),trEe_tr,alpha
      real*8 Ee0_tr(3,3),Me_tr(3,3),Me0_tr(3,3),pbar,tauBar_tr,Np(3,3)
      real*8 theta_tau,args(nargs),Dp_tau(3,3),Dp_vec(3,3),Dp_eig(3),nu0
      real*8 expdtDP(3,3),Me_tau(3,3),Fe_tau(3,3),Re_tau(3,3),tmp,detF
      real*8 Ue_tau(3,3),Ee_tau(3,3),Fp_tau_inv(3,3),Ctilde(3,3,3,3)
      real*8 dYdnu,dYdgBar,factor1,factor2,ratio,tauBar_tau,psi_cr
      real*8 StressTempJac(3,3),PlWrkDefJac(3,3),kk,erff,lc,zeta
      real*8 trEe_tau,Ee0_tau(3,3),trEe0_tau_Ee0_tau,dplasticwork
      real*8 Ee0_tau_Ee0_tau(3,3),dum1,Histmax_tau,Histmax_t
      real*8 rpl,drpldt

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third,six
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0,
     +     six=6.d0)

      !
      ! Identity matrix
      !
      call onem(Iden)
      !
      ! Read material properties
      !
      Eyoung   = props(1)
      poisson  = props(2)
      Y0       = props(3)
      H        = props(4)
      Ysat     = props(5)
      mRate    = props(6)
      nu0      = props(7)
      erff = props(8)
      lc = props(9)
      zeta = props(10)
      psi_cr = props(11)
      !
      ! Calculate elastic modulii
      !
      Gshear = Eyoung/(two*(one+poisson))
      Kbulk  = Eyoung/(three*(one-two*poisson))
      Lambda = Kbulk - two_third*Gshear
      !
      kk     = 1d-4
      !
      ! Compute the trial elastic distortion
      !
      call m3inv(Fp_t,Fp_t_inv)
      Fe_tr = matmul(F_tau,Fp_t_inv)
      !
      ! Compute the volume ratio
      !
      call mdet(F_tau,detF)
      !
      ! Perform kinematical calculations for the trial state
      !
      call skinem(Fe_tr,Re_tr,Ue_tr,Ee_tr)
      !
      ! Compute the trace of the trial strain, and the deviator
      !  of the trial strain
      !
      trEe_tr = Ee_tr(1,1) + Ee_tr(2,2) + Ee_tr(3,3)
      Ee0_tr = Ee_tr - third*trEe_tr*Iden
      !
      ! Compute the trial Mandel stress
      !
      Me_tr = ((one-theta_tau)**two+kk)*(
     +      two*Gshear*Ee0_tr + Kbulk*(trEe_tr)*Iden)
      !
      ! Compute the mean normal pressure (this is the value
      !  at the end of the increment as well as the trial)
      !
      pbar = -third*(Me_tr(1,1) + Me_tr(2,2) + Me_tr(3,3))
      !
      ! Compute the trial Mandel stress deviator
      !
      Me0_tr = Me_tr + pbar*Iden
      !
      ! Compute the trial equiv. shear stress in mechanism 1
      !
      tauBar_tr = dsqrt(one/two)*dsqrt(sum(Me0_tr*Me0_tr))
      !
      ! Compute the direction of plastic flow in mechanism 1
      !  (this is the trial, as well as the tau value)
      !
      if(tauBar_tr.gt.zero) then
         Np = Me0_tr/(dsqrt(two)*tauBar_tr)
      else
         Np = Iden
      endif
      !
      ! Solve the implicit relation for the equiv. plastic
      !  shearing rate in mechanism 1 at the end of the increment
      !
      if(tauBar_tr.le.zero) then
         nuP_tau = zero
      else
         args(1)   = tauBar_tr
         args(2)   = Gshear
         args(3)   = dtime
         args(4)   = Y_t
         args(5)   = mRate
         args(6)   = nu0
         args(7)   = H
         args(8)   = Ysat
         args(9)   = theta_tau
         args(10)  = kk
         !
         ! call the ``fail-safe'' solver
         !
         call solveRate(nuP_t,args,nargs,nuP_tau)
         !
      endif
      !
      ! Compute the plastic stretching at the end of the increment
      !
      Dp_tau = (one/dsqrt(two))*nuP_tau*Np
      !
      ! Compute the equiv. plastic shear strain at the
      !  end of the increment
      !
      gBarP_tau = gBarP_t + dtime*nuP_tau
      !
      ! Compute the plastic distortion at the end of
      !  the increment using the exponential map
      !
      if(nuP_tau.le.zero) then
         Fp_tau = Fp_t
      else
         call spectral(dtime*Dp_tau,Dp_eig,Dp_vec)
         expdtDp = zero
         expdtDp(1,1) = dexp(Dp_eig(1))
         expdtDp(2,2) = dexp(Dp_eig(2))
         expdtDp(3,3) = dexp(Dp_eig(3))
         expdtDp = matmul(matmul(Dp_vec,expdtDp),transpose(Dp_vec))
         Fp_tau = matmul(expdtDp,Fp_t)
      endif
      !
      ! Check to make sure that det(Fp_tau)>0
      !
      call mdet(Fp_tau,tmp)
      if(tmp.le.zero) then
         write(*,*) 'det(Fp_tau).le.zero in INTEG'
         call xit
      endif
      !
      ! Compute the elastic distortion at the end
      !  of the increment
      !
      call m3inv(Fp_tau,Fp_tau_inv)
      Fe_tau = matmul(F_tau,Fp_tau_inv)
      !
      ! Perform kinematical calculations for the trial state
      !
      call skinem(Fe_tau,Re_tau,Ue_tau,Ee_tau)
      !
      trEe_tau = Ee_tau(1,1) + Ee_tau(2,2) + Ee_tau(3,3)
      Ee0_tau = Ee_tau - third*trEe_tau*Iden
      !
      ! Compute the useful quantity E0:E0
      !
      Ee0_tau_Ee0_tau= matmul(Ee0_tau,Ee0_tau)
      trEe0_tau_Ee0_tau = Ee0_tau_Ee0_tau(1,1) +
     +   Ee0_tau_Ee0_tau(2,2) + Ee0_tau_Ee0_tau(3,3)
      !
      ! Update History function
      !
      dum1 = Gshear * trEe0_tau_Ee0_tau
     +    +   half * Kbulk * trEe_tau**two
      !
      if ((dum1.gt.Histmax_t).and.(detF.gt.one)) then
          !
          ! If current history function is maximum, update the history function
          !
          Histmax_tau = dum1
      else 
          Histmax_tau = Histmax_t
      end if
      !
      ! Compute the Mandel stress at the end of the increment
      !
      Me_tau = Me_tr - two*Gshear*dtime*Dp_tau
      !
      ! Compute the equiv. shear stress at the end of
      !  the increment, and its ratio to the trial
      !
      tauBar_tau = tauBar_tr - Gshear*dtime*nuP_tau
      if(tauBar_tr.gt.zero) then
         ratio = tauBar_tau/tauBar_tr
      else
         ratio = one
      endif
      !
      ! Compute the Cauchy stress at the end of the increment
      !
      T_tau = matmul(Re_tau,matmul(Me_tau,transpose(Re_tau)))/detF
      !
      ! Compute the resistance at the end of the increment
      !
      call GetState(Y_t,dtime,H,Ysat,nuP_tau,Y_tau,dYdnu,dYdgBar)
      !
      ! Compute the modified elasticity tensor, Ctilde
      !
      Ctilde = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Ctilde(i,j,k,l) = Ctilde(i,j,k,l) 
     +                 + Gshear*ratio*
     +                 (Iden(i,k)*Iden(j,l) + Iden(i,l)*Iden(j,k))
     +                 + (Kbulk-(two/three)*Gshear)*Iden(i,j)*Iden(k,l)
               enddo
            enddo
         enddo
      enddo
      !
      ! Compute the material tangent
      !
      if((tauBar_tr.gt.zero).and.(nuP_tau.gt.zero)) then
         !
         factor1 = two*Gshear*(one - ratio)
         factor2 = (two*Gshear*Gshear*dtime)/
     +        (
     +        Gshear*dtime 
     +        + dYdnu*((nuP_tau/nu0)**mRate) 
     +        + (Y_tau*mRate/nuP_tau)*((nuP_tau/nu0)**mRate)
     +        )
         !
         MatTan = zero
         do i=1,3
            do j=1,3
               do k=1,3
                  do l=1,3
                     MatTan(i,j,k,l) = Ctilde(i,j,k,l)
     +                    + factor1*Np(i,j)*Np(k,l)
     +                    - factor2*Np(i,j)*Np(k,l)
                  enddo
               enddo
            enddo
         enddo
      else
         MatTan = Ctilde
      endif
      MatTan = ((one-theta_tau)**two+kk) * MatTan
      !
      ! Compute the source term and its derivative (Appendix A.3)
      !
      if (Histmax_tau.gt.psi_cr) then
      rpl = (two*(one - theta_tau)*Histmax_tau)/detF
     +    - (erff*theta_tau)/detF
      drpldt = (-two*Histmax_tau) - (erff)
      !
      else
          rpl = zero
          drpldt = zero
      endif
      !
      return
      end subroutine integ

***********************************************************************

      subroutine getState(Y_t,dtime,H,Ysat,nuP,Y_tau,dYdnu,dYdgBar)

      implicit none

      integer i,j,k,l

      real*8 Y_t,Y_tau,nuP,H,dtime,one,two,dYdnu,dYdgBar,Ysat
      parameter(one=1.d0,two=2.d0)




      ! The deformation resistance
      !
      Y_tau = (Y_t + dtime*H*Ysat*nuP)/(one + dtime*H*nuP)
      !
      ! The derivative with respect to nuP
      !
      dYdnu = (dtime*H*(Ysat - Y_t))/((one + dtime*H*nuP)**two)
      !
      ! The derivative with respect to DgBar=dt*nuP
      !
      dYdgBar = (H*(Ysat - Y_t))/((one + dtime*H*nuP)**two)
      !

      return
      end subroutine getState


***********************************************************************

      subroutine solveRate(rootOld,args,nargs,root)

      ! This subroutine will numerically solve for the equiv.
      !  plastic shear rate. See numerical recipies RTSAFE.

      implicit none

      integer maxit,i,j,nargs

      real*8 xacc,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,args(nargs),zero,rootOld,temp,rootMax,rootMin

      parameter(maxit=50)
      parameter(xacc=1.d-6,zero=0.d0,one=1.d0)

      rootMax = 0.d0
      rootMin = 1.d8

      x1 = rootMin
      x2 = rootMax
      call gfunc(x1,FL,DF,args,nargs)
      call gfunc(x2,FH,DF,args,nargs)

      if(fl*fh.ge.zero) then
         write(*,*) 'FYI, root not bracketed on nuP'
         write(*,*) 'fl=',fl
         write(*,*) 'x1=',x1
         write(*,*) 'fh=',fh
         write(*,*) 'x2=',x2
         write(*,*) 'args='
         do i=1,nargs
            write(*,*) args(i)
         enddo
         call xit
         return
      endif

C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C
      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
      if(rootOld.lt.rootMin) rootOld = rootMin
      if(rootOld.gt.rootMax) rootOld = rootMax
      ROOT = rootOld !0.5D0 *( X1 + X2)
      DXOLD = DABS(X2 - X1)
      DX = DXOLD
      
      call gfunc(root,F,DF,args,nargs)

C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GE. 0.D0
     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN

            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            IF( XL .EQ. ROOT ) THEN
C
C			CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX
            IF( TEMP .EQ. ROOT) THEN
C
C			 CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         END IF
C
C		CONVERVEGENCE CRITERION
C
         IF( DABS(DX) .LT. XACC) RETURN

C
C			THE ONE NEW FUNCTION EVALUATION PER ITERATION
C
         call gfunc(root,F,DF,args,nargs)

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE

      WRITE(*,'(/1X,A)') 'solveRate EXCEEDING MAXIMUM ITERATIONS'
      WRITE(80,'(/1X,A)') 'solveRate EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine solveRate

****************************************************************************

      subroutine gfunc(nuP,f,df,args,nargs)

      implicit none

      integer nargs

      real*8 args(nargs),f,df,tauBar_tr,Gshear,Y_t,Y_tau,dtime,mRate,H
      real*8 Ysat,dYdnu,dYdgBar,nu0,nuP,aux,theta_tau,kk

      real*8 zero,one,two,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)


      if(nuP.le.zero) then
         f  = one
         df = zero
         return
      endif

      
      ! Obtain relevant quantities
      !
      tauBar_tr = args(1)
      Gshear    = args(2)
      dtime     = args(3)
      Y_t       = args(4)
      mRate     = args(5)
      nu0       = args(6)
      H         = args(7)
      Ysat      = args(8)
      theta_tau = args(9)
      kk = args(10)


      ! Compute the resistance for this nuP
      !
      call getState(Y_t,dtime,H,Ysat,nuP,Y_tau,dYdnu,dYdgBar)


      ! Compute the residual
      !
!      f = tauBar_tr - dtime*Gshear*nuP - Y_tau*((nuP/nu0)**mRate)
      f = tauBar_tr - dtime*Gshear*nuP - Y_tau*((nuP/nu0)**mRate) *
     +     ((one-theta_tau)**two+kk)


      ! Compute the tangent
      !
      if((nuP/nu0).gt.zero) then
!         df = -Gshear*dtime - dYdnu*((nuP/nu0)**mRate) 
         df = -Gshear*dtime - dYdnu*((nuP/nu0)**mRate) *
     +    ((one-theta_tau)**two+kk)
      else
         df = -dtime*Gshear
      endif


      return
      end subroutine gfunc

****************************************************************************

      subroutine jac2D(SpTanMod,ddsdde)

      implicit none

      real*8 SpTanMod(3,3,3,3),ddsdde(4,4)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)

      end subroutine jac2D

***********************************************************************

      subroutine jac3D(SpTanMod,ddsdde)

      implicit none

      real*8 SpTanMod(3,3,3,3),ddsdde(6,6)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)
      ddsdde(1,5) = SpTanMod(1,1,1,3)
      ddsdde(1,6) = SpTanmod(1,1,2,3)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)
      ddsdde(2,5) = SpTanMod(2,2,1,3)
      ddsdde(2,6) = SpTanmod(2,2,2,3)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)
      ddsdde(3,5) = SpTanMod(3,3,1,3)
      ddsdde(3,6) = SpTanmod(3,3,2,3)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)
      ddsdde(4,5) = SpTanMod(1,2,1,3)
      ddsdde(4,6) = SpTanmod(1,2,2,3)

      ddsdde(5,1) = SpTanMod(1,3,1,1)
      ddsdde(5,2) = SpTanMod(1,3,2,2)
      ddsdde(5,3) = SpTanMod(1,3,3,3)
      ddsdde(5,4) = SpTanMod(1,3,1,2)
      ddsdde(5,5) = SpTanMod(1,3,1,3)
      ddsdde(5,6) = SpTanmod(1,3,2,3)

      ddsdde(6,1) = SpTanMod(2,3,1,1)
      ddsdde(6,2) = SpTanMod(2,3,2,2)
      ddsdde(6,3) = SpTanMod(2,3,3,3)
      ddsdde(6,4) = SpTanMod(2,3,1,2)
      ddsdde(6,5) = SpTanMod(2,3,1,3)
      ddsdde(6,6) = SpTanmod(2,3,2,3)

      end subroutine jac3D

***********************************************************************
c
c
c  The following are all utility routines used in fortran codes
c
c
c
C**********************************************************************
C	THE NEXT SUBROUTINE CALCULATES VARIOUS KINEMATICAL QUANTITIES 
C	ASSOCIATED WITH THE DEFORMATION GRADIENT
C**********************************************************************
	SUBROUTINE SKINEM(F,R,U,E)

C	THIS SUBROUTINE PERFORMS THE RIGHT POLAR DECOMPOSITION
C	[F] = [R][U] OF THE DEFORMATION GRADIENT [F] INTO
C	A ROTATION [R] AND THE RIGHT  STRETCH TENSOR [U].
C	THE EIGENVALUES AND EIGENVECTORS OF [U] AND
C	THE LOGARITHMIC STRAIN [E] = LN [U]
C	ARE ALSO RETURNED.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION F(3,3),FTRANS(3,3), C(3,3), OMEGA(3),
     +	          UEIGVAL(3),EIGVEC(3,3), EIGVECT(3,3), 
     +            U(3,3),E(3,3),UINV(3,3),R(3,3),TEMPM(3,3)
     
	COMMON/ERRORINFO/UMERROR     

C	F(3,3)	-- THE DEFORMATION GRADIENT MATRIX WHOSE
C		   POLAR DECOMPOSITION IS DESIRED.
C	DETF	-- THE DETRMINANT OF [F]; DETF > 0.
C	FTRANS(3,3)	-- THE TRANSPOSE OF [F].
C	R(3,3)	-- THE ROTATION MATRIX; [R]^T [R] = [I];
C		   OUTPUT.
C	U(3,3)	-- THE RIGHT STRETCH TENSOR; SYMMETRIC
C		   AND POSITIVE DEFINITE; OUTPUT.
C	UINV(3,3)	-- THE INVERSE OF [U].
C	C(3,3)	-- THE RIGHT CAUCHY-GREEN TENSOR = [U][U];
C		   SYMMETRIC AND POSITIVE DEFINITE.
C	OMEGA(3)-- THE SQUARES OF THE PRINCIPAL STRETCHES.
C 	UEIGVAL(3)	-- THE PRINCIPAL STRETCHES; OUTPUT.
C	EIGVEC(3,3)	-- MATRIX OF EIGENVECTORS OF [U];OUTPUT.
C	EIGVECT(3,3)    -- TRANSPOSE OF THE ABOVE.
C	E(3,3)	-- THE LOGARITHMIC STRAIN TENSOR, [E]=LN[U];
C		   OUTPUT.
C**********************************************************************

C	STORE THE IDENTITY MATRIX IN  [R], [U], AND [UINV]

	CALL ONEM(R)
	CALL ONEM(U)
	CALL ONEM(UINV)

C	STORE THE ZERO MATRIX IN [E]

	CALL ZEROM(E)

C      	CHECK IF THE DETERMINANT OF [F] IS GREATER THAN ZERO.
C	IF NOT, THEN PRINT DIAGNOSTIC AND STOP.

        CALL MDET(F,DETF)
        IF (DETF .LE. 0.D0) THEN
          WRITE(*,100)
          UMERROR=5.
          RETURN
        ENDIF

C      	CALCULATE THE RIGHT CAUCHY GREEN STRAIN TENSOR [C]

        CALL  MTRANS(F,FTRANS)
        CALL  MPROD(FTRANS,F,C)
 
C	CALCULATE THE EIGENVALUES AND EIGENVECTORS OF  [C]

	CALL SPECTRAL(C,OMEGA,EIGVEC)

C	CALCULATE THE PRINCIPAL VALUES OF [U] AND [E]

	UEIGVAL(1) = DSQRT(OMEGA(1))
	UEIGVAL(2) = DSQRT(OMEGA(2))
	UEIGVAL(3) = DSQRT(OMEGA(3))

	U(1,1) = UEIGVAL(1)
	U(2,2) = UEIGVAL(2)
	U(3,3) = UEIGVAL(3)

	E(1,1) = DLOG( UEIGVAL(1) )
	E(2,2) = DLOG( UEIGVAL(2) )
	E(3,3) = DLOG( UEIGVAL(3) )

C	CALCULATE THE COMPLETE TENSORS [U] AND [E]

	CALL MTRANS(EIGVEC,EIGVECT)
	CALL MPROD(EIGVEC,U,TEMPM)
	CALL MPROD(TEMPM,EIGVECT,U)
	CALL MPROD(EIGVEC,E,TEMPM)
	CALL MPROD(TEMPM,EIGVECT,E)

C	CALCULATE [UINV]

	CALL M3INV(U,UINV)

C	CALCULATE [R]

	CALL MPROD(F,UINV,R)
100     FORMAT(5X,'--ERROR IN KINEMATICS-- THE DETERMINANT OF [F]',
     +         ' IS NOT GREATER THAN 0')

	RETURN
	END

C**********************************************************************
C	THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
C	DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
C**********************************************************************
	SUBROUTINE SPECTRAL(A,D,V)
C
C	THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C	A SYMMETRIC 3 BY 3 MATRIX [A]. 
C
C	THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C	EIGENVALUES IN ASCENDING ORDER, AND
C	A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
C	EIGENVECTORS.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER(NP=3)
	DIMENSION D(NP),V(NP,NP)
	DIMENSION A(3,3),E(NP,NP)

	DO 2 I = 1,3
          DO 1 J= 1,3
            E(I,J) = A(I,J)
1	  CONTINUE
2	CONTINUE

	CALL JACOBI(E,3,NP,D,V,NROT)
	CALL EIGSRT(D,V,3,NP)

	RETURN
	END

C**********************************************************************
	SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C	COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C	MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
C	NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
C	ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C	AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C	VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
C	[V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C	[A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
C	EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
C	WHICH WERE REQUIRED.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NMAX =100)
	DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C	INITIALIZE [V] TO THE IDENTITY MATRIX

	DO 12 IP = 1,N	
	  DO 11 IQ = 1,N
	    V(IP,IQ) = 0.D0
11        CONTINUE
          V(IP,IP) = 1.D0
12	CONTINUE

C	INITIALIZE [B] AND [D] TO THE DIAGONAL OF [A], AND Z TO ZERO.
C	THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C	IN EQUATION (11.1.14)

	DO 13 IP = 1,N
	  B(IP) = A(IP,IP)
	  D(IP) = B(IP)
	  Z(IP) = 0.D0
13	CONTINUE
C
	NROT = 0
	DO 24 I = 1,50

C	SUM OFF-DIAGONAL ELEMENTS

          SM = 0.D0
          DO 15 IP = 1, N-1
            DO 14 IQ = IP + 1, N
	      SM = SM + DABS ( A(IP,IQ ))
14          CONTINUE
15        CONTINUE

C	IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
C	WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
C	UNDERFLOW.

          IF ( SM .EQ. 0.D0) RETURN
C
C	  IF ( SM .LT. 1.0D-15) RETURN

C	IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C	|A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, 
C	SEE EQUATION (11.1.25). THEREAFTER TRESH = 0.

          IF ( I .LT. 4) THEN
            TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = 0.D0
          ENDIF
C
          DO 22 IP = 1, N-1
            DO 21 IQ = IP+1,N
              G = 100.D0*DABS(A(IP,IQ))

C	AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C	OFF-DIAGONAL ELEMENT IS SMALL.

	      IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP)))
     +            .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G .EQ. DABS(H)) THEN

C	T = 1./(2.*THETA), EQUATION(11.1.10)

	          T =A(IP,IQ)/H
	        ELSE
	          THETA = 0.5D0*H/A(IP,IQ)
	          T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
	          IF (THETA .LT. 0.D0) T = -T
	        ENDIF
	        C = 1.D0/DSQRT(1.D0 + T**2)
	        S = T*C
	        TAU = S/(1.D0 + C)
	        H = T*A(IP,IQ)
	        Z(IP) = Z(IP) - H
	        Z(IQ) = Z(IQ) + H
	        D(IP) = D(IP) - H
	        D(IQ) = D(IQ) + H
	        A(IP,IQ) = 0.D0

C	CASE OF ROTATIONS 1 <= J < P
				
	        DO 16 J = 1, IP-1
	          G = A(J,IP)
	          H = A(J,IQ)
	          A(J,IP) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
16	        CONTINUE

C	CASE OF ROTATIONS P < J < Q

	        DO 17 J = IP+1, IQ-1
	          G = A(IP,J)
	          H = A(J,IQ)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
17	        CONTINUE

C	CASE OF ROTATIONS Q < J <= N

	        DO 18 J = IQ+1, N
                  G = A(IP,J)
	          H = A(IQ,J)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(IQ,J) = H + S*(G - H*TAU)
18	        CONTINUE
	        DO 19 J = 1,N
	          G = V(J,IP)
	          H = V(J,IQ)
	          V(J,IP) = G - S*(H + G*TAU)
	          V(J,IQ) = H + S*(G - H*TAU)
19	        CONTINUE
	        NROT = NROT + 1
              ENDIF
21	    CONTINUE
22	  CONTINUE

C	UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

	  DO 23 IP = 1, N
	    B(IP) = B(IP) + Z(IP)
	    D(IP) = B(IP)
	    Z(IP) = 0.D0
23	  CONTINUE
24	CONTINUE

C	IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
C	THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
C	AND STOP.

	WRITE (*,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'

	RETURN
	END

C**********************************************************************
	SUBROUTINE EIGSRT(D,V,N,NP)

C	GIVEN THE EIGENVALUES [D] AND EIGENVECTORS [V] AS OUTPUT FROM
C	JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
C	AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(NP),V(NP,NP)

	DO 13 I = 1,N-1
	  K = I
	  P = D(I)
	  DO 11 J = I+1,N
	    IF (D(J) .GE. P) THEN
	      K = J
	      P = D(J)
	    END IF
11	  CONTINUE
	  IF (K .NE. I) THEN
	    D(K) = D(I)
	    D(I) = P
	    DO 12 J = 1,N
	      P = V(J,I)
	      V(J,I) = V(J,K)
	      V(J,K) = P
12	    CONTINUE
  	  ENDIF
13	CONTINUE

	RETURN
	END

C**********************************************************************
C	THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
C**********************************************************************
	SUBROUTINE ZEROV(V,SIZE)

C	THIS SUBROUTINE STORES THE ZERO VECTOR IN A VECTOR V
C	OF SIZE SIZE.
C**********************************************************************

	INTEGER SIZE
	REAL*8 V(0:SIZE-1)

	DO 1 I=0,SIZE
	  V(I) = 0.D0
1	CONTINUE
	
	RETURN
	END

C**********************************************************************
      	SUBROUTINE ZEROM(A)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.D0.
C**********************************************************************

        REAL*8 A(3,3)

	DO 1 I=1,3
	  DO 1 J=1,3
	    A(I,J) = 0.D0
1	CONTINUE
C	
	RETURN
	END

C**********************************************************************
	SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

	DO 1 I=1,3
	  DO 1 J=1,3
	    IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MTRANS(A,ATRANS)
 
C	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

	REAL*8 A(3,3),ATRANS(3,3)

	DO 1 I=1,3
 	  DO 1 J=1,3
	    ATRANS(J,I) = A(I,J)
1	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MPROD(A,B,C)
 
C 	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C(3,3)

	DO 2 I = 1, 3
	  DO 2 J = 1, 3
	    C(I,J) = 0.D0
	    DO 1 K = 1, 3
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE
C
	RETURN
	END

C**********************************************************************
	SUBROUTINE MPROD4(A,B,C)
 
C	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(4,4),B(4,4),C(4,4)

	DO 2 I = 1, 4
   	  DO 2 J = 1, 4
	    C(I,J) = 0.D0
	    DO 1 K = 1, 4
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE DOTPM(A,B,C)

C	THIS SUBROUTINE CALCULATES THE SCALAR PRODUCT OF TWO
C	3 BY 3 MATRICES [A] AND [B] AND STORES THE RESULT IN THE
C	SCALAR C.
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C

	C = 0.D0
	DO 1 I = 1,3
	  DO 1 J = 1,3
            C = C + A(I,J)*B(I,J)
1	CONTINUE
C
	RETURN
	END

C**********************************************************************
	SUBROUTINE MDET(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].
C**********************************************************************

	REAL*8  A(3,3), DET

	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	        + A(1,2)*A(2,3)*A(3,1)
     +	        + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************
	SUBROUTINE M3INV(A,AINV)

C 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	[A] AND PLACES THE RESULT IN [AINV]. 
C 	IF DET(A) IS ZERO, THE CALCULATION
C 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	DET		-- THE COMPUTED DETERMINANT OF [A].
C	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C			   IS CALLED THE COFACTOR OF A(I,J).
C	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C			   OBTAINED BY REPLACING EACH ELEMENT OF
C			   [A] BY ITS COFACTOR, AND THEN TAKING
C			   TRANSPOSE OF THE RESULTING MATRIX.
C	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C			   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

	CALL MDET(A,DET)
	IF ( DET .EQ. 0.D0 ) THEN
	  write(*,10)
	  STOP
	ENDIF
	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE
10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

	RETURN
	END

C**********************************************************************
	SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

	REAL*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************
        SUBROUTINE INVAR(A,IA,IIA,IIIA)

C	THIS SUBROUTINE CALCULATES THE PRINCIPAL INVARIANTS 
C	IA, IIA, IIIA OF A TENSOR [A].
C**********************************************************************

        REAL*8 A(3,3), AD(3,3),AD2(3,3), DETA, IA,IIA,IIIA

        DO 1 I=1,3
          DO 1 J=1,3
            AD(I,J) = A(I,J)
1       CONTINUE
        IA = AD(1,1) + AD(2,2) + AD(3,3)

C	CALCULATE THE SQUARE OF [AD]

        CALL MPROD(AD,AD,AD2)
        IIA =0.5D0 * ( IA*IA - ( AD2(1,1) + AD2(2,2) + AD2(3,3) ) )

        CALL  MDET(AD,DETA)
        IIIA = DETA

        RETURN
        END

C**********************************************************************
	SUBROUTINE TRACEM(A,TRA)

C	THIS SUBROUTINE CALCULATES THE TRACE OF A 3 BY 3 MATRIX [A]
C	AND STORES THE RESULT IN THE SCALAR TRA
C**********************************************************************

	REAL*8 A(3,3),TRA

	TRA = A(1,1) + A(2,2) + A(3,3)

	RETURN 
	END

C**********************************************************************
	SUBROUTINE DEVM(A,ADEV)

C	THIS SUBROUTINE CALCULATES THE DEVIATORIC PART OF A
C	3 BY 3 MATRIX [A]
C**********************************************************************

	REAL*8 A(3,3),TRA,ADEV(3,3),IDEN(3,3)

	CALL TRACEM(A,TRA)
	CALL ONEM(IDEN)
	CALL ZEROM(ADEV)

	DO 1 I = 1,3
	  DO 1 J = 1,3
	    ADEV(I,J) = A(I,J) - (1.D0/3.D0)*TRA*IDEN(I,J)
1	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE EQUIVS(S,SB)

C	THIS SUBROUTINE CALCULATES THE EQUIVALENT SHEAR STRESS SB
C	CORRESPONDING TO A 3 BY 3 STRESS MATRIX [S]
C**********************************************************************

	REAL*8 S(3,3),SDEV(3,3),SDOTS,SB

	SB = 0.D0
	SDOTS = 0.D0

	CALL DEVM(S,SDEV)
	CALL DOTPM(SDEV,SDEV,SDOTS)
	SB = DSQRT(1.5D0*SDOTS)

	RETURN
	END
C **********************************************************************
        SUBROUTINE PRESSURE(A,PRA)
C
C       THIS SUBROUTINE CALCULATES THE MEAN NORMAL PRESSURE
C       OF A 3 BY 3 MATRIX [A]
C       AND STORES THE RESULT IN THE SCALAR PRA
C ----------------------------------------------------------------------
C       VARIABLES
C
        REAL*8 A(3,3),PRA

        PRA = -(1.D0 / 3.D0)*( A(1,1) + A(2,2) + A(3,3) )

        RETURN 
        END
C**********************************************************************
	SUBROUTINE PRTMAT(A,M,N)
C**********************************************************************

	INTEGER M,N
	REAL*8 A(M,N)	  

	DO 10 K=1,M
	  WRITE(80,'(2X,6E12.4,2X)') (A(K,L), L=1,N)
10      CONTINUE

        RETURN
        END

C**********************************************************************
	SUBROUTINE PRTVEC(A,M)
C**********************************************************************

	INTEGER M
	REAL*8 A(M)	  

	WRITE(80,'(2X,6E12.4,2X)') (A(K), K=1,M)

        RETURN
	END
C*************************************************************************	  
c***********************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
C
C	Given an NxN matrix [A], with physical dimension NP, this 
C	routine replaces it by the LU decomposition of a row-wise 
C	permutation of itself. [A] and N are input. [A] is output, 
C	arranged in LU form. INDX is an output vector which records
C	the row permutation effected by the partial pivoting; 
C	D is output as +1 or -1 depending on wheter the nuber of
C	row interchanges was even or odd, respectively. This routine
C	is used in combination with LUBKSB to solve linear equations 
C	or invert a matrix.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	
	PARAMETER (NMAX=100,TINY=1.0E-20)
	DIMENSION A(NP,NP),INDX(N),VV(NMAX)
	D=1.
	DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12	CONTINUE
	DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16	CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19	CONTINUE
	IF(A(N,N).EQ.0.)A(N,N)=TINY
	RETURN
	END
C**********************************************************
       SUBROUTINE LUBKSB(A,N,NP,INDX,B)
C
C	Solves the set of N linear equations [A]{X} = {B}. 
C	Here [A] is input, not as the matrix [A], but as its LU 
C	decomposition, determined by the routine LUDCMP. INDX
C	is input as the permutation vector returned by LUDCMP. {B}
C	is input as the right-hand side vector {B}, and returns
C	with the solution vector {X}. [A], N, NP, INDX are not 
C	modified by this routine, and can be left in place
C	for succesive calls with different right-hand sides {B}.
C	This routine takes into account that {B} will begin with
C	many zero elements, so it is efficient for use in matrix 
C	inversion.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION A(NP,NP),INDX(N),B(N)
       II=0
       DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12     CONTINUE
       DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14     CONTINUE
       RETURN
       END
c**********************************************************
