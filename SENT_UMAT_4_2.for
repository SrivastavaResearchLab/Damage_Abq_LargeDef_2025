************************************************************************
!     A finite element implementation of a large deformation gradient-
!     damage phase-field theory for fracture using Abaqus user material subroutines
!     Keven Alkhoury, Shawn A. Chester, and Vikas Srivastava
!     Keven Alkhoury, September 2025, Abaqus 2024, In review
************************************************************************
!     State Variables
!     --------------------------------------------------------------
!     statev(1) = Histmax_tau ------ Referential history function
!     statev(2) = detF_tau ------ Determinant of the deformation gradient
!     statev(3) = F_tau(1,1) ------ Deformation gradient (1,1)
!     statev(4) = F_tau(2,2) ------ Deformation gradient (2,2)
!     statev(5) = F_tau(3,3) ------ Deformation gradient (3,3)
!     statev(6) = F_tau(1,2) ------ Deformation gradient (1,2)
!     statev(7) = F_tau(1,3) ------ Deformation gradient (1,3)
!     statev(8) = F_tau(2,1) ------ Deformation gradient (2,1)
!     statev(9) = F_tau(2,3) ------ Deformation gradient (2,3)
!     statev(10) = F_tau(3,1) ------ Deformation gradient (3,1)
!     statev(11) = F_tau(3,2) ------ Deformation gradient (3,2)
!     --------------------------------------------------------------
!     Material Properties Vector
!     --------------------------------------------------------------
!     Mechanical (UMAT):
!     --------------------------------------------------------------
!     Youngs= props(1)  ------ Young's moudulus
!     Poisson = props(2) ------ Poisson ratio
!     erff=props(3)  -- Fracture energy/unit reference volume
!     lc=props(4) ------- Intrinsic characteristic length scale
!     zeta=props(5) ------ Viscous regularization parameter
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

      integer i,j,ii,jj,iterError,stat

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),theta_t,theta_tau,Eyoung
      real*8 poisson,Gshear,Kbulk,lambda,H_tau(3,3),E_tau(3,3),trE
      real*8 E0_tau(3,3),T_tau(3,3),plasticwork,Cmat(6,6)
      real*8 erff,lc,zeta,trE0_tau_E0_tau,E0_tau_E0_tau(3,3)
      real*8 Finv(3,3),Finv_t(3,3),detF_tau,detF_t,dum1
      real*8 Histmax_tau,Histmax_t,kk,rpl,drpldt
      !
      ! Parameters
      !
      real*8 zero,one,two,half,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0)
      !
      ! Identity matrix
      !
      Iden = zero
      do i=1,3
         Iden(i,i) = one
      enddo
      !
      ! Obtain old and new deformation gradients
      !
      F_t = dfgrd0
      F_tau = dfgrd1
      !
      ! Obtain old and new temperature/damage
      !
      theta_t = temp
      theta_tau = temp + dtemp
      !
      ! Obtain state variables from the previous increment
      !
      if((kinc.le.1).and.(kstep.eq.1)) then
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
          Histmax_t = statev(1)
          detF_t = statev(2)
          F_t(1,1) = statev(3)
          F_t(2,2) = statev(4)
          F_t(3,3) = statev(5)
          F_t(1,2) = statev(6)
          F_t(1,3) = statev(7)
          F_t(2,1) = statev(8)
          F_t(2,3) = statev(9)
          F_t(3,1) = statev(10)
          F_t(3,2) = statev(11)
      endif
      !
      ! Obtain material properties
      !
      Eyoung = props(1)
      poisson = props(2)
      erff = props(3)
      lc = props(4)
      zeta = props(5)
      !
      Gshear = Eyoung/(two*(one+poisson))
      Kbulk = Eyoung/(three*(one-two*poisson))
      lambda = Kbulk - (two/three)*Gshear
      !
      kk = 1d-4
      !
      ! Compute the total strain (small strain) and its deviator
      !
      H_tau = F_tau - Iden
      E_tau = half*(H_tau + transpose(H_tau))
      trE = E_tau(1,1) + E_tau(2,2) + E_tau(3,3)
      E0_tau = E_tau - third*trE*Iden
      !
      ! Compute the inverse of F, its determinant, and its transpose
      !
      call matInv3D(F_tau,Finv,detF_tau,stat)
      call matInv3D(F_t,Finv_t,detF_t,stat)
      !
      ! Compute the useful quantity E0:E0
      !
      E0_tau_E0_tau= matmul(E0_tau,E0_tau)
      trE0_tau_E0_tau = E0_tau_E0_tau(1,1) +
     +   E0_tau_E0_tau(2,2) + E0_tau_E0_tau(3,3)
      !
      ! Update History function
      !
      dum1 = Gshear * trE0_tau_E0_tau
     +    +   half * Kbulk * trE**two
      !
      if (dum1.gt.Histmax_t) then
          !
          ! If current history function is maximum, update the history function
          !
          Histmax_tau = dum1
      else 
          Histmax_tau = Histmax_t
      end if
      !
      ! Cauchy stress
      !
      T_tau = ((one-theta_tau)**two+kk) *
     +    (two*Gshear*E0_tau + Kbulk*trE*Iden)
      !
      ! Compute the source term and its derivative (Appendix A.3)
      !
      rpl = (two*(one - theta_tau)*Histmax_tau)/detF_tau
     +    - (erff*theta_tau)/detF_tau
      !
      drpldt = (-two*Histmax_tau) - (erff)
      !
      ! Update state variables at the end of the increment
      !
      statev(1) = Histmax_tau
      statev(2) = detF_tau
      statev(3) = F_tau(1,1)
      statev(4) = F_tau(2,2)
      statev(5) = F_tau(3,3)
      statev(6) = F_tau(1,2)
      statev(7) = F_tau(1,3)
      statev(8) = F_tau(2,1)
      statev(9) = F_tau(2,3)
      statev(10) = F_tau(3,1)
      statev(11) = F_tau(3,2)
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
      ! Calculate the isotropic elastic stiffness matrix and use
      !  it for the stress-deformation tangent
      !
      Cmat = 0.d0
      do ii=1,3
         Cmat(ii,ii) = (Eyoung*(one-poisson))/
     +        ((one+poisson)*(one-two*poisson))
      enddo
      do ii=4,6
         Cmat(ii,ii) = Eyoung/(two*(one+poisson))
      enddo
      Cmat(1,2) = (Eyoung*poisson)/((one+poisson)*(one-two*poisson))
      Cmat(1,3) = Cmat(1,2)
      Cmat(2,1) = Cmat(1,2)
      Cmat(2,3) = Cmat(1,2)
      Cmat(3,1) = Cmat(1,2)
      Cmat(3,2) = Cmat(1,2)
      !
      Cmat = ((one-theta_tau)**two+kk)*Cmat
      !
      ! Return Abaqus/Standard the stress-deformation tangent matrix
      !
      ddsdde = Cmat
      !
      return
      end subroutine umat

***********************************************************************
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
      detF_tau = statev(2)
      !
      ! Used to update conduction term (Appendix A.2)
      !
      aux = (erff*(lc**two))
      !
      ! Remember, these are F_tau as UMAT is called before UMATHT
      !
      F(1,1) = statev(3)
      F(2,2) = statev(4)
      F(3,3) = statev(5)
      F(1,2) = statev(6)
      F(1,3) = statev(7)
      F(2,1) = statev(8)
      F(2,3) = statev(9)
      F(3,1) = statev(10)
      F(3,2) = statev(11)
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
      !
!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************
