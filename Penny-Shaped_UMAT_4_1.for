************************************************************************
!     Alkhoury, K., Chester, S. A., & Srivastava, V. (2026). A finite
!	  element implementation of a large deformation gradient-damage
!	  theory for fracture with Abaqus user material subroutines.
! 	  Engineering Fracture Mechanics, 331, 111677. 
!     https://doi.org/10.1016/j.engfracmech.2025.111677
!     Keven Alkhoury, November 2025, Abaqus 2024
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
!     Gshear= props(1)  ------ Shear moudulus
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

      dimension stress(ntens+3),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character*80 cmname,file1
      character*256 jobName,outDir,fileName

      integer i,j,k,l,iterError,lenJobName,lenOutDir

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),theta_t,theta_tau,T_tau(3,3)
      real*8 Gshear,effStr,detF,dTRdF(3,3,3,3),SpTanMod(3,3,3,3),trB
      real*8 Finv(3,3),Bdis(3,3),trBdis,Bdis0(3,3),B_tau(3,3),Poisson
      real*8 effStr_t,effStr_tau,psi_t,psi_tau,H_t,H_tau,g,beta
      real*8 B_t(3,3),trB_t,TR_tau(3,3),DTKDF(3,3,3,3),dum1
      real*8 FinvT(3,3),Histmax_tau,Histmax_t,kk,Gc,xl,trB_tau
      real*8 Gshear0,detF_t,detF_tau,Finv_t(3,3),SpHistmax_t
      real*8 SpHistmax_tau,StressTempJac(3,3),erff,lc,JBAR,KBULK,DKBULK
      !
      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine
      !
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,six=6.d0,eight=8.d0)
      !
      ! Identity matrix
      !
      call onem(Iden)
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
      Gshear0 = props(1)
      Poisson = props(2)
      erff = props(3)
      lc = props(4)
      zeta = props(5)
      !
      beta = two*Poisson/(one-two*Poisson)
      !
      kk     = 1d-4
      !
      ! Compute the inverse of F, its determinant, and its transpose
      !
      call matInv3D(F_tau,Finv,detF_tau,stat)
      FinvT = transpose(Finv)
      !
      call matInv3D(F_t,Finv_t,detF_t,stat)
      !
      ! Left Green Cauchy tensor
      !
      B_tau = matmul(F_tau,transpose(F_tau))
      trB_tau = B_tau(1,1) + B_tau(2,2) + B_tau(3,3)
      !
      B_t = matmul(F_t,transpose(F_t))
      trB_t = B_t(1,1) + B_t(2,2) + B_t(3,3)
      !
      ! Update History function
      !
      dum1 = Gshear0/two*(trB_t-three)
     +    + Gshear0/beta*(detF_t**(-beta)-one)
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
      T_tau = ((one-theta_tau)**two+kk)*Gshear0*
     +         (B_tau-detF_tau**(-beta)*Iden)/detF_tau
      !
      ! First Piola stress
      !
      TR_tau = ((one-theta_tau)**two+kk)*Gshear0*
     +         (F_tau-detF_tau**(-beta)*FinvT)
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
      ! Compute material tangent modulus, dTRdF
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3          
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l) + 
     +                ((one-theta_tau)**two+kk)*(
     +                Gshear0*(Iden(i,k)*Iden(j,l)
     +              + beta*detF_tau**(-beta)*Finv(l,k)*Finv(j,i)
     +              + detF_tau**(-beta)*Finv(l,i)*Finv(j,k)))
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate dTKdF based on dTRdF, where dTKdF is the tangent
      !  modulus of the Kirchhoff stress with respect to the deformation
      !  gradient
      !
      dTKdF = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     dTKdF(i,j,k,l) = dTKdF(i,j,k,l) +
     +                    dTRdF(i,m,k,l)*F_tau(j,m) +
     +                    TR_tau(i,m)*Iden(j,k)*Iden(m,l)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate material jacobian
      ! Calculate the equilibrium tangent modulus based on dTKdF
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) +
     +                    (half/detF_tau)*
     +                    (
     +                    F_tau(l,m)*dTKdF(i,j,k,m) + 
     +                    F_tau(k,m)*dTKdF(i,j,l,m)
     +                    )
                  enddo
               enddo
            enddo
         enddo
      enddo
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
         call jac3D(SpTanMod,ddsdde)
      elseif(ntens.eq.4) then
         call jac2D(SpTanMod,ddsdde)
      endif
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
      subroutine jac2D(SpTanMod,ddsdde)

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
