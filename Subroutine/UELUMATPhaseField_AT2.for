! UEL & UMAT subroutines for the phase field model for fracture  
! Linear elasticity version, suitable for both plane strain and 3D     

! If using this code for research or industrial purposes, please cite:
! Wei Tan, Emilio Martinez-Paneda. 
! Phase field predictions of microscopic fracture and R-curve behaviour of fibre-reinforced composites
! (2021):108539 
! doi: https://doi.org/10.1016/j.compscitech.2020.108539

! Emilio Martinez-Paneda (e.martinez-paneda@imperial.ac.uk)
! Imperial College London

! Wei Tan (wei.tan@qmul.ac.uk)
! Queen Mary University of London



      module ktransf
      implicit none
      real*8 UserVar(4,2,400000) !CPE4 or CPE8R 
      !real*8 UserVar(8,2,400000) !C3D8 or C3D20R
      !real*8 UserVar(4,2,400000) !C3D10
      integer nelem,kincK,kkiter,kflagE
      save
      end module
      
!***********************************************************************
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      include 'aba_param.inc' !implicit real(a-h o-z)
      dimension time(2)
      if (lop.eq.0) then !start of analysis
       call mutexinit(1)
      endif
      return
      end
      
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

      use ktransf
      include 'aba_param.inc' !implicit real(a-h o-z)

      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

      parameter(ndim=2,ninpt=4,nsvint=2,ndof=1) !CPE4 or CPE8R
      !parameter(ndim=3,ninpt=8,nsvint=2,ndof=1) !C3D8 or C3D20R
      !parameter(ndim=3,ninpt=4,nsvint=2,ndof=1) !C3D10
      
      dimension wght(ninpt),dN(nnode,1),dNdz(ndim,nnode),
     1 dNdx(ndim,nnode),statevLocal(nsvint)
        
!     initialising
      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0
      if (ninpt.eq.4.and.nnode.eq.10.and.ndim.eq.3) then
       wght=0.25d0/6.d0
      else
       wght=1.d0
      endif     
      
!     reading parameters
      xlc=props(1)
      Gc=props(2)
      zeta=props(3)
      kflag=props(4)

!     AT1 or AT2 flag (if flag=1, then AT1, else AT2)
      kflagAT=2
      
!     viscous dissipation and iteration counter
      if (kflag.eq.1) then
       if (jelem.eq.1) then
        if (kinc.ne.kincK) then    
         kincK=kinc
         kkiter=1
        else
         kkiter=kkiter+1 
        endif
       endif
      endif 

      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       dvol=wght(kintk)*djac

!     recover and assign state variables       
       call kstatevar(kintk,nsvint,svars,statevLocal,1)
       Hn=statevLocal(1)
       H=statevLocal(2)

!     calculate phase field and incremental strains from nodal values
       phi=0.d0
       dphi=0.d0
       do inod=1,nnode
        phi=phi+dN(inod,1)*u(inod)
        dphi=dphi+dN(inod,1)*du(inod,1)
       end do
       
!      AT1 / AT2 variables
       if (kflagAT.eq.1) then ! AT1
        Hmin=3.d0*Gc/(16.d0*xlc)
        ATpar=1.d0
        c0 = 8.d0/3.d0
       else ! AT2
        Hmin=0.d0
        ATpar=0.d0
        c0= 2.d0
       endif

!     enforcing Karush-Kuhn-Tucker conditions
       H=max(H,Hmin)
       if (H.lt.Hn) H=Hn
       
!     collect information from UMAT and store state variables
       statevLocal(1)=H
       statevLocal(2)=UserVar(kintk,2,jelem-nelem)
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
       
       if (kflag.eq.1) then
        if (kkiter.le.3) zeta=0.d0
       elseif (kflag.eq.2) then
        if (dphi.ge.0.d0) zeta=0.d0           
       endif
       
       dalph=ATpar+2.d0*(1.d0-ATpar)*phi
       ddalph=2.d0*(1.d0-ATpar)
       dGdeg=-2.d0+2.d0*phi
       
!     form and assemble stiffness matrix and internal force vector
       amatrx(1:ndofel,1:ndofel)=amatrx(1:ndofel,1:ndofel)+
     1 dvol*(matmul(transpose(dNdx),dNdx)*2.d0*Gc*xlc/c0+
     2 matmul(dN,transpose(dN))*(Gc/xlc*ddalph/c0+zeta/dtime+2.d0*H))

       rhs(1:ndofel,1)=rhs(1:ndofel,1)-dvol*
     1 (matmul(transpose(dNdx),matmul(dNdx,u(1:ndofel)))*2.d0*Gc*xlc/c0
     2 +dN(:,1)*(dGdeg*H+Gc/xlc*dalph/c0+zeta*dphi/dtime))

!     information transfer to UMAT
       UserVar(kintk,1,jelem-nelem)=phi
      
      end do       ! end loop on material integration points

      RETURN
      END

!***********************************************************************      
      subroutine kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.577350269d0)
      parameter (gca=0.5854101966d0, gcb=0.1381966012d0)      
      dimension dN(nnode,1),dNdz(ndim,*),coord24(2,4),coord34(3,4),
     * coord38(3,8)
      
      data  coord24 /-1.d0, -1.d0,
     2                1.d0, -1.d0,
     3               -1.d0,  1.d0,
     4                1.d0,  1.d0/  
      
      data  coord34 /gcb, gcb, gcb,
     1               gca, gcb, gcb,
     2               gcb, gca, gcb,
     3               gcb, gcb, gca/       
      
      data  coord38 /-1.d0, -1.d0, -1.d0,
     2                1.d0, -1.d0, -1.d0,
     3               -1.d0,  1.d0, -1.d0,
     4                1.d0,  1.d0, -1.d0, 
     5               -1.d0, -1.d0,  1.d0,
     6                1.d0, -1.d0,  1.d0,
     7               -1.d0,  1.d0,  1.d0,
     8                1.d0,  1.d0,  1.d0/
      
      if (ninpt.eq.4.and.nnode.eq.4.and.ndim.eq.2) then ! CPE4
          
!     determine (g,h)
       g=coord24(1,kintk)*gaussCoord
       h=coord24(2,kintk)*gaussCoord

!     shape functions 
       dN(1,1)=(1.d0-g)*(1.d0-h)/4.d0
       dN(2,1)=(1.d0+g)*(1.d0-h)/4.d0
       dN(3,1)=(1.d0+g)*(1.d0+h)/4.d0
       dN(4,1)=(1.d0-g)*(1.d0+h)/4.d0

!     derivative d(Ni)/d(g)
       dNdz(1,1)=-(1.d0-h)/4.d0
       dNdz(1,2)=(1.d0-h)/4.d0
       dNdz(1,3)=(1.d0+h)/4.d0
       dNdz(1,4)=-(1.d0+h)/4.d0

!     derivative d(Ni)/d(h)
       dNdz(2,1)=-(1.d0-g)/4.d0
       dNdz(2,2)=-(1.d0+g)/4.d0
       dNdz(2,3)=(1.d0+g)/4.d0
       dNdz(2,4)=(1.d0-g)/4.d0
       
      elseif (ninpt.eq.4.and.nnode.eq.8.and.ndim.eq.2) then ! CPE8R 
          
!     determine (g,h,r)
       g=coord24(1,kintk)*gaussCoord
       h=coord24(2,kintk)*gaussCoord

!     shape functions 
       dN(1,1)=-0.25d0*(1.d0-g)*(1.d0-h)*(1.d0+g+h)
       dN(2,1)=0.25d0*(1.d0+g)*(1.d0-h)*(g-h-1.d0)
       dN(3,1)=0.25d0*(1.d0+g)*(1.d0+h)*(g+h-1.d0)
       dN(4,1)=0.25d0*(1.d0-g)*(1.d0+h)*(h-g-1.d0)
       dN(5,1)=0.5d0*(1.d0-g*g)*(1.d0-h)
       dN(6,1)=0.5d0*(1.d0+g)*(1.d0-h*h)
       dN(7,1)=0.5d0*(1.d0-g*g)*(1.d0+h)
       dN(8,1)=0.5d0*(1.d0-g)*(1.d0-h*h)        

!     derivative d(Ni)/d(g)
       dNdz(1,1)=0.25d0*(1.d0-h)*(2.d0*g+h)
       dNdz(1,2)=0.25d0*(1.d0-h)*(2.d0*g-h)
       dNdz(1,3)=0.25d0*(1.d0+h)*(2.d0*g+h)
       dNdz(1,4)=0.25d0*(1.d0+h)*(2.d0*g-h)
       dNdz(1,5)=-g*(1.d0-h)
       dNdz(1,6)=0.5d0*(1.d0-h*h)
       dNdz(1,7)=-g*(1.d0+h)
       dNdz(1,8)=-0.5d0*(1.d0-h*h)      

!     derivative d(Ni)/d(h)
       dNdz(2,1)=0.25d0*(1.d0-g)*(g+2.d0*h)
       dNdz(2,2)=0.25d0*(1.d0+g)*(2.d0*h-g)
       dNdz(2,3)=0.25d0*(1.d0+g)*(2.d0*h+g)
       dNdz(2,4)=0.25d0*(1.d0-g)*(2.d0*h-g)
       dNdz(2,5)=-0.5d0*(1.d0-g*g) 
       dNdz(2,6)=-(1.d0+g)*h 
       dNdz(2,7)=0.5d0*(1.d0-g*g)
       dNdz(2,8)=-(1.d0-g)*h           
          
      elseif (ninpt.eq.4.and.nnode.eq.10.and.ndim.eq.3) then ! C3D10
       f=coord34(1,kintk) ! Multiply here by for consistency
       g=coord34(2,kintk)
       h=coord34(3,kintk)
       
!     shape functions 
       dN(1,1)=(1.d0-f-g-h)*(1.d0-2.d0*f-2.d0*g-2.d0*h)
       dN(2,1)=f*(2.d0*f-1.d0)
       dN(3,1)=g*(2.d0*g-1.d0)
       dN(4,1)=h*(2.d0*h-1.d0)
       dN(5,1)=4.d0*f*(1.d0-f-g-h)
       dN(6,1)=4.d0*f*g
       dN(7,1)=4.d0*g*(1.d0-f-g-h)
       dN(8,1)=4.d0*h*(1.d0-f-g-h)
       dN(9,1)=4.d0*h*f
       dN(10,1)=4.d0*g*h

!     derivative d(Ni)/d(f)       
       dNdz(1,1)=-3.d0+4.d0*f+4.d0*g+4.d0*h
       dNdz(1,2)=4.d0*f-1.d0
       dNdz(1,3)=0.d0
       dNdz(1,4)=0.d0
       dNdz(1,5)=4.d0*(1.d0-f-g-h)-4.d0*f
       dNdz(1,6)=4.d0*g
       dNdz(1,7)=-4.d0*g
       dNdz(1,8)=-4.d0*h
       dNdz(1,9)=4.d0*h
       dNdz(1,10)=0.d0

!     derivative d(Ni)/d(g)
       dNdz(2,1)=-3.d0+4.d0*f+4.d0*g+4.d0*h
       dNdz(2,2)=0.d0
       dNdz(2,3)=4.d0*g-1.d0
       dNdz(2,4)=0.d0
       dNdz(2,5)=-4.d0*f
       dNdz(2,6)=4.d0*f
       dNdz(2,7)=4.d0*(1.d0-f-g-h)-4.d0*g
       dNdz(2,8)=-4.d0*h
       dNdz(2,9)=0.d0
       dNdz(2,10)=4.d0*h
      
!     derivative d(Ni)/d(h)
       dNdz(3,1)=-3.d0+4.d0*f+4.d0*g+4.d0*h
       dNdz(3,2)=0.d0
       dNdz(3,3)=0.d0
       dNdz(3,4)=4.d0*h-1.d0
       dNdz(3,5)=-4.d0*f
       dNdz(3,6)=0.d0
       dNdz(3,7)=-4.d0*g
       dNdz(3,8)=4.d0*(1.d0-f-g-h)-4.d0*h
       dNdz(3,9)=4.d0*f
       dNdz(3,10)=4.d0*g   
          
      elseif (ninpt.eq.8.and.nnode.eq.8.and.ndim.eq.3) then ! C3D8

!     determine (g,h,r)     
       f=coord38(1,kintk)*gaussCoord
       g=coord38(2,kintk)*gaussCoord
       h=coord38(3,kintk)*gaussCoord

!     shape functions
       dN(1,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0-h)
       dN(2,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0-h)
       dN(3,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0-h)
       dN(4,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0-h)
       dN(5,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0+h)
       dN(6,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0+h)
       dN(7,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0+h)
       dN(8,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0+h)     

!     derivative d(Ni)/d(f)
       dNdz(1,1)=-0.125d0*(1.d0-g)*(1.d0-h)
       dNdz(1,2)= 0.125d0*(1.d0-g)*(1.d0-h)
       dNdz(1,3)= 0.125d0*(1.d0+g)*(1.d0-h)
       dNdz(1,4)=-0.125d0*(1.d0+g)*(1.d0-h)
       dNdz(1,5)=-0.125d0*(1.d0-g)*(1.d0+h)
       dNdz(1,6)= 0.125d0*(1.d0-g)*(1.d0+h)
       dNdz(1,7)= 0.125d0*(1.d0+g)*(1.d0+h)
       dNdz(1,8)=-0.125d0*(1.d0+g)*(1.d0+h)

!     derivative d(Ni)/d(g)
       dNdz(2,1)=-0.125d0*(1.d0-f)*(1.d0-h)
       dNdz(2,2)=-0.125d0*(1.d0+f)*(1.d0-h)
       dNdz(2,3)= 0.125d0*(1.d0+f)*(1.d0-h)
       dNdz(2,4)= 0.125d0*(1.d0-f)*(1.d0-h)
       dNdz(2,5)=-0.125d0*(1.d0-f)*(1.d0+h)
       dNdz(2,6)=-0.125d0*(1.d0+f)*(1.d0+h)
       dNdz(2,7)= 0.125d0*(1.d0+f)*(1.d0+h)
       dNdz(2,8)= 0.125d0*(1.d0-f)*(1.d0+h)
      
!     derivative d(Ni)/d(h)
       dNdz(3,1)=-0.125d0*(1.d0-f)*(1.d0-g)
       dNdz(3,2)=-0.125d0*(1.d0+f)*(1.d0-g)
       dNdz(3,3)=-0.125d0*(1.d0+f)*(1.d0+g)
       dNdz(3,4)=-0.125d0*(1.d0-f)*(1.d0+g)
       dNdz(3,5)= 0.125d0*(1.d0-f)*(1.d0-g)
       dNdz(3,6)= 0.125d0*(1.d0+f)*(1.d0-g)
       dNdz(3,7)= 0.125d0*(1.d0+f)*(1.d0+g)
       dNdz(3,8)= 0.125d0*(1.d0-f)*(1.d0+g)
       
      elseif (ninpt.eq.8.and.nnode.eq.20.and.ndim.eq.3) then ! C3D20R       
       
!     determine (g,h,r)  
       f=coord38(1,kintk)*gaussCoord
       g=coord38(2,kintk)*gaussCoord
       h=coord38(3,kintk)*gaussCoord

!     shape functions       
       dN(9,1)= 0.25d0*(1.d0-f**2.d0)*(1.d0-g)*(1.d0-h)
       dN(10,1)=0.25d0*(1.d0+f)*(1.d0-g**2.d0)*(1.d0-h)
       dN(11,1)=0.25d0*(1.d0-f**2.d0)*(1.d0+g)*(1.d0-h)
       dN(12,1)=0.25d0*(1.d0-f)*(1.d0-g**2.d0)*(1.d0-h)
       dN(13,1)=0.25d0*(1.d0-f**2.d0)*(1.d0-g)*(1.d0+h)
       dN(14,1)=0.25d0*(1.d0+f)*(1.d0-g**2.d0)*(1.d0+h)
       dN(15,1)=0.25d0*(1.d0-f**2.d0)*(1.d0+g)*(1.d0+h)
       dN(16,1)=0.25d0*(1.d0-f)*(1.d0-g**2.d0)*(1.d0+h)
       dN(17,1)=0.25d0*(1.d0-f)*(1.d0-g)*(1.d0-h**2.d0)
       dN(18,1)=0.25d0*(1.d0+f)*(1.d0-g)*(1.d0-h**2.d0)
       dN(19,1)=0.25d0*(1.d0+f)*(1.d0+g)*(1.d0-h**2.d0)
       dN(20,1)=0.25d0*(1.d0-f)*(1.d0+g)*(1.d0-h**2.d0)
       dN(1,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0-h)
     & -(dN(9,1)+dN(12,1)+dN(17,1))/2.d0
       dN(2,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0-h)
     & -(dN(9,1)+dN(10,1)+dN(18,1))/2.d0
       dN(3,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0-h)
     & -(dN(10,1)+dN(11,1)+dN(19,1))/2.d0
       dN(4,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0-h)
     & -(dN(11,1)+dN(12,1)+dN(20,1))/2.d0
       dN(5,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0+h)
     & -(dN(13,1)+dN(16,1)+dN(17,1))/2.d0
       dN(6,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0+h)
     & -(dN(13,1)+dN(14,1)+dN(18,1))/2.d0
       dN(7,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0+h)
     & -(dN(14,1)+dN(15,1)+dN(19,1))/2.d0
       dN(8,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0+h)  
     & -(dN(15,1)+dN(16,1)+dN(20,1))/2.d0

!     derivative d(Ni)/d(f)
       dNdz(1,9)=-0.5d0*f*(1.d0-g)*(1.d0-h)
       dNdz(1,10)=0.25d0*(1.d0-g**2.d0)*(1.d0-h)
       dNdz(1,11)=-0.5d0*f*(1.d0+g)*(1.d0-h)
       dNdz(1,12)=-0.25d0*(1.d0-g**2.d0)*(1.d0-h)
       dNdz(1,13)=-0.5d0*f*(1.d0-g)*(1.d0+h)
       dNdz(1,14)=0.25d0*(1.d0-g**2.d0)*(1.d0+h)
       dNdz(1,15)=-0.5d0*f*(1.d0+g)*(1.d0+h)
       dNdz(1,16)=-0.25d0*(1.d0-g**2.d0)*(1.d0+h)
       dNdz(1,17)=-0.25d0*(1.d0-g)*(1.d0-h**2.d0)
       dNdz(1,18)=0.25d0*(1.d0-g)*(1.d0-h**2.d0)
       dNdz(1,19)=0.25d0*(1.d0+g)*(1.d0-h**2.d0)
       dNdz(1,20)=-0.25d0*(1.d0+g)*(1.d0-h**2.d0)
       dNdz(1,1)=-0.125d0*(1.d0-g)*(1.d0-h)
     & -(dNdz(1,9)+dNdz(1,12)+dNdz(1,17))/2.d0
       dNdz(1,2)=0.125d0*(1.d0-g)*(1.d0-h)
     & -(dNdz(1,9)+dNdz(1,10)+dNdz(1,18))/2.d0
       dNdz(1,3)=0.125d0*(1.d0+g)*(1.d0-h)
     & -(dNdz(1,10)+dNdz(1,11)+dNdz(1,19))/2.d0
       dNdz(1,4)=-0.125d0*(1.d0+g)*(1.d0-h)
     & -(dNdz(1,11)+dNdz(1,12)+dNdz(1,20))/2.d0
       dNdz(1,5)=-0.125d0*(1.d0-g)*(1.d0+h)
     & -(dNdz(1,13)+dNdz(1,16)+dNdz(1,17))/2.d0
       dNdz(1,6)=0.125d0*(1.d0-g)*(1.d0+h)
     & -(dNdz(1,13)+dNdz(1,14)+dNdz(1,18))/2.d0
       dNdz(1,7)=0.125d0*(1.d0+g)*(1.d0+h)
     & -(dNdz(1,14)+dNdz(1,15)+dNdz(1,19))/2.d0
       dNdz(1,8)=-0.125d0*(1.d0+g)*(1.d0+h)
     & -(dNdz(1,15)+dNdz(1,16)+dNdz(1,20))/2.d0

!     derivative d(Ni)/d(g)
       dNdz(2,9)=-0.25d0*(1.d0-f**2.d0)*(1.d0-h)
       dNdz(2,10)=-0.5d0*g*(1.d0+f)*(1.d0-h)
       dNdz(2,11)=0.25d0*(1.d0-f**2.d0)*(1.d0-h)
       dNdz(2,12)=-0.5d0*g*(1.d0-f)*(1.d0-h)
       dNdz(2,13)=-0.25d0*(1.d0-f**2.d0)*(1.d0+h)
       dNdz(2,14)=-0.5d0*g*(1.d0+f)*(1.d0+h)
       dNdz(2,15)= 0.25d0*(1.d0-f**2.d0)*(1.d0+h)
       dNdz(2,16)=-0.5d0*g*(1.d0-f)*(1.d0+h)
       dNdz(2,17)=-0.25d0*(1.d0-f)*(1.d0-h**2.d0)
       dNdz(2,18)=-0.25d0*(1.d0+f)*(1.d0-h**2.d0)
       dNdz(2,19)=0.25d0*(1.d0+f)*(1.d0-h**2.d0)
       dNdz(2,20)=0.25d0*(1.d0-f)*(1.d0-h**2.d0)
       dNdz(2,1)=-0.125d0*(1.d0-f)*(1.d0-h)
     & -(dNdz(2,9)+dNdz(2,12)+dNdz(2,17))/2.d0
       dNdz(2,2)=-0.125d0*(1.d0+f)*(1.d0-h)
     & -(dNdz(2,9)+dNdz(2,10)+dNdz(2,18))/2.d0
       dNdz(2,3)=0.125d0*(1.d0+f)*(1.d0-h)
     & -(dNdz(2,10)+dNdz(2,11)+dNdz(2,19))/2.d0
       dNdz(2,4)=0.125d0*(1.d0-f)*(1.d0-h)
     & -(dNdz(2,11)+dNdz(2,12)+dNdz(2,20))/2.d0
       dNdz(2,5)=-0.125d0*(1.d0-f)*(1.d0+h)
     & -(dNdz(2,13)+dNdz(2,16)+dNdz(2,17))/2.d0
       dNdz(2,6)=-0.125d0*(1.d0+f)*(1.d0+h)
     & -(dNdz(2,13)+dNdz(2,14)+dNdz(2,18))/2.d0
       dNdz(2,7)=0.125d0*(1.d0+f)*(1.d0+h)
     & -(dNdz(2,14)+dNdz(2,15)+dNdz(2,19))/2.d0
       dNdz(2,8)=0.125d0*(1.d0-f)*(1.d0+h)
     & -(dNdz(2,15)+dNdz(2,16)+dNdz(2,20))/2.d0
      
!     derivative d(Ni)/d(h)
       dNdz(3,9)=-0.25d0*(1.d0-f**2.d0)*(1.d0-g)
       dNdz(3,10)=-0.25d0*(1.d0+f)*(1.d0-g**2.d0)
       dNdz(3,11)=-0.25d0*(1.d0-f**2.d0)*(1.d0+g)
       dNdz(3,12)=-0.25d0*(1.d0-f)*(1.d0-g**2.d0)
       dNdz(3,13)=0.25d0*(1.d0-f**2.d0)*(1.d0-g)
       dNdz(3,14)=0.25d0*(1.d0+f)*(1.d0-g**2.d0)
       dNdz(3,15)=0.25d0*(1.d0-f**2.d0)*(1.d0+g)
       dNdz(3,16)=0.25d0*(1.d0-f)*(1.d0-g**2.d0)
       dNdz(3,17)=-0.5d0*h*(1.d0-f)*(1.d0-g)
       dNdz(3,18)=-0.5d0*h*(1.d0+f)*(1.d0-g)
       dNdz(3,19)=-0.5d0*h*(1.d0+f)*(1.d0+g)
       dNdz(3,20)=-0.5d0*h*(1.d0-f)*(1.d0+g)
       dNdz(3,1)=-0.125d0*(1.d0-f)*(1.d0-g)
     & -(dNdz(3,9)+dNdz(3,12)+dNdz(3,17))/2.d0
       dNdz(3,2)=-0.125d0*(1.d0+f)*(1.d0-g)
     & -(dNdz(3,9)+dNdz(3,10)+dNdz(3,18))/2.d0
       dNdz(3,3)=-0.125d0*(1.d0+f)*(1.d0+g)
     & -(dNdz(3,10)+dNdz(3,11)+dNdz(3,19))/2.d0
       dNdz(3,4)=-0.125d0*(1.d0-f)*(1.d0+g)
     & -(dNdz(3,11)+dNdz(3,12)+dNdz(3,20))/2.d0
       dNdz(3,5)=0.125d0*(1.d0-f)*(1.d0-g)
     & -(dNdz(3,13)+dNdz(3,16)+dNdz(3,17))/2.d0
       dNdz(3,6)=0.125d0*(1.d0+f)*(1.d0-g)
     & -(dNdz(3,13)+dNdz(3,14)+dNdz(3,18))/2.d0
       dNdz(3,7)=0.125d0*(1.d0+f)*(1.d0+g)
     & -(dNdz(3,14)+dNdz(3,15)+dNdz(3,19))/2.d0
       dNdz(3,8)=0.125d0*(1.d0-f)*(1.d0+g)
     & -(dNdz(3,15)+dNdz(3,16)+dNdz(3,20))/2.d0       
       
      else
       write (6,*) '***ERROR: The shape fuctions cannot be found'   
      endif    
 
      return
      end

!***********************************************************************
      subroutine kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix
!     dNdx - shape functions derivatives w.r.t. global coordinates      
      include 'aba_param.inc'

      dimension xjac(ndim,ndim),xjaci(ndim,ndim),coords(mcrd,nnode),
     1 dNdz(ndim,nnode),dNdx(ndim,nnode)      

      xjac=0.d0

      do inod=1,nnode
       do idim=1,ndim
        do jdim=1,ndim
         xjac(jdim,idim)=xjac(jdim,idim)+
     1        dNdz(jdim,inod)*coords(idim,inod)      
        end do
       end do 
      end do

      if (ndim.eq.3) then
          
       djac=xjac(1,1)*xjac(2,2)*xjac(3,3)+xjac(2,1)*xjac(3,2)*xjac(1,3)
     & +xjac(3,1)*xjac(2,3)*xjac(1,2)-xjac(3,1)*xjac(2,2)*xjac(1,3)
     & -xjac(2,1)*xjac(1,2)*xjac(3,3)-xjac(1,1)*xjac(2,3)*xjac(3,2)
       if (djac.gt.0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=(xjac(2,2)*xjac(3,3)-xjac(2,3)*xjac(3,2))/djac
        xjaci(1,2)=(xjac(1,3)*xjac(3,2)-xjac(1,2)*xjac(3,3))/djac
        xjaci(1,3)=(xjac(1,2)*xjac(2,3)-xjac(1,3)*xjac(2,2))/djac
        xjaci(2,1)=(xjac(2,3)*xjac(3,1)-xjac(2,1)*xjac(3,3))/djac
        xjaci(2,2)=(xjac(1,1)*xjac(3,3)-xjac(1,3)*xjac(3,1))/djac
        xjaci(2,3)=(xjac(1,3)*xjac(2,1)-xjac(1,1)*xjac(2,3))/djac
        xjaci(3,1)=(xjac(2,1)*xjac(3,2)-xjac(2,2)*xjac(3,1))/djac
        xjaci(3,2)=(xjac(1,2)*xjac(3,1)-xjac(1,1)*xjac(3,2))/djac
        xjaci(3,3)=(xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))/djac
       else ! negative or zero jacobian
        write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
       endif 
          
      else if (ndim.eq.2) then 
          
       djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
       if (djac.gt.0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=xjac(2,2)/djac
        xjaci(2,2)=xjac(1,1)/djac
        xjaci(1,2)=-xjac(1,2)/djac
        xjaci(2,1)=-xjac(2,1)/djac
       else ! negative or zero jacobian
        write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
       endif
       
      endif
      
      dNdx=matmul(xjaci,dNdz)
			
      return
      end

!***********************************************************************
      subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)
c
c     Transfer data to/from element-level state variable array from/to
c     material-point level state variable array.
c
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc=(npt-1)*nsvint     ! integration point increment

      if (icopy.eq.1) then ! Prepare arrays for entry into umat
       do i=1,nsvint
        statev_ip(i)=statev(i+isvinc)
       enddo
      else ! Update element state variables upon return from umat
       do i=1,nsvint
        statev(i+isvinc)=statev_ip(i)
       enddo
      end if

      return
      end

!***********************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use ktransf
      include 'aba_param.inc' !implicit real(a-h o-z)

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)
      
      dimension Edev(ntens),eprin(3)
      
!     find number of elements 
      if (npt.eq.1) then
       if ((time(1)-dtime).lt.0.d0) then
        if (kflagE.ne.13) then
         kincK=0 
         kkiter=1 
         nelem=noel
         UserVar=0.d0
         kflagE=13
        else
         CALL MutexLock(1)   
         if (noel.gt.nelem) nelem=noel
         CALL MutexUnlock(1)
        endif 
       endif
      endif

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      xk=1.d-07

!     Amor (kflagD=1) or Miehe [only 2D] (kflagD=2) decomposition flag      
      kflagD=0
      
      phi=statev(ntens+1)

!     Build stiffness matrix
      eg=E/(1.d0+xnu)/2.d0
      elam=E*xnu/((1.d0+xnu)*(1.d0-2.d0*xnu))
      bk=elam+eg*2.d0/real(ntens/2)
      
!     Update stresses
      do k1=1,3
       do k2=1,3
        ddsdde(k2,k1)=elam
       end do
       ddsdde(k1,k1)=2.d0*eg+elam
      end do
      do k1=4,ntens
       ddsdde(k1,k1)=eg
      end do     
      
      stran=stran+dstran
      stress=matmul(ddsdde,stran)
      
!     compute strain energy density
      if (kflagD.eq.1) then
       Edev=stran
       trE=stran(1)+stran(2)+stran(3)
       Edev(1:3)=stran(1:3)-trE/3.d0
       EdevS=Edev(1)**2+Edev(2)**2+Edev(3)**2
       do k1=1,ntens-3
        EdevS=EdevS+Edev(3+k1)**2/2.d0
       end do
       trEp=0.5d0*(trE+abs(trE))
       trEn=0.5d0*(trE-abs(trE))
       psip=0.5d0*bk*trEp**2+eg*EdevS
       psin=0.5d0*bk*trEn**2
      elseif (kflagD.eq.2) then
       eprin(1)=(stran(1)+stran(2))/2+sqrt(((stran(1)-stran(2))/2)**2
     1 +(stran(4)/2)**2)
       eprin(2)=(stran(1)+stran(2))/2-sqrt(((stran(1)-stran(2))/2)**2
     1 +(stran(4)/2)**2)
       eprin(3)=stran(3)
    
       trp1=(eprin(1)+eprin(2)+eprin(3)+abs(eprin(1)+eprin(2)+
     1 eprin(3)))/2.d0
       trn1=(eprin(1)+eprin(2)+eprin(3)-abs(eprin(1)+eprin(2)+
     1 eprin(3)))/2.d0
       trp2=0.d0
       trn2=0.d0
       do i=1,3
        trp2=trp2 + (eprin(i)+abs(eprin(i)))**2.d0/4.d0
        trn2=trn2 + (eprin(i)-abs(eprin(i)))**2.d0/4.d0
       end do 
       psip=props(2)*eg/(1d0-2d0*props(2))*trp1**2d0+eg*trp2
       psin=props(2)*eg/(1d0-2d0*props(2))*trn1**2d0+eg*trn2
      else
       psip=0.d0
       do i=1,ntens
        psip=psip+stress(i)*stran(i)*0.5d0
       end do
      endif
      
!     Degradation function
      g=(1.d0-phi)**2+xk
      stress=stress*g
      ddsdde=ddsdde*g      

!     Collect information from UEL 
      phi=UserVar(npt,1,noel)

!     information transfer to UEL
      UserVar(npt,2,noel)=psip

!     output      
      statev(1:ntens)=stress
      statev(ntens+1)=phi
      statev(ntens+2)=psip
     
      return
      end