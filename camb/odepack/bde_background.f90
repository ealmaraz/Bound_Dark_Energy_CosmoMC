program bde_background

!The following is a simple example problem, with the coding
!needed for its solution by DLSODA.  The problem is from chemical
!kinetics, and consists of the following three rate equations:
!     dy1/dt = -.04*y1 + 1.e4*y2*y3
!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!     dy3/dt = 3.e7*y2**2
! on the interval from t = 0.0 to t = 4.e10, with initial conditions
! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.

! The following coding solves this problem with DLSODA,
! printing results at t = .4, 4., ..., 4.e10.  It uses
! ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
! y2 has much smaller values.
! At the end of the run, statistical quantities of interest are
! printed (see optional outputs in the full description below).

  external bde_equations, jac

  double precision atol, rtol, rwork, t, tout, y
  dimension y(3), atol(3), rwork(70), iwork(23)

  integer neq, jt
  
  neq = 3
  
  y(1) = 0.2729322530
  y(2) = 0.1929922469
  y(3) = 0.3687988897D+26
  
  t = 0
  tout = 0.1290659433E-11

  itol = 2
  rtol = 1.D-04
  
  atol(1) = 1.D-8
  atol(2) = 1.D-8
  atol(3) = 1.D-8
  
  itask = 1
  istate = 1
  iopt = 0
  
  lrw = 70
  liw = 23
  
  jt = 1

  CALL DLSODA(bde_equations,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,jt)
  WRITE(*,*) t,y(1),y(2),y(3)
  write(*,*) 'istate', istate

  tout = tout*1.E+09
  CALL DLSODA(bde_equations,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,jt)
  WRITE(*,*) t,y(1),y(2),y(3)
  write(*,*) 'istate', istate
  
  tout = 0.1290530367E+02
  CALL DLSODA(bde_equations,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,jt)
  WRITE(*,*) t,y(1),y(2),y(3)
  write(*,*) 'istate', istate

  
!     DO 40 IOUT = 1,12
!       CALL DLSODA(bde_equations,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)
!       WRITE(6,20) t,y(1),y(2),y(3)
! 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6)
!       IF (ISTATE .LT. 0) GO TO 80
! 40    TOUT = TOUT*10.
!     WRITE(6,60) iwork(11),iwork(12),iwork(13),iwork(19),iwork(15)
! 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4/&
!        ' Method last used =',I2,'   Last switch was at t =',D12.4)
!     STOP
! 80  WRITE(6,90)istate
! 90  FORMAT(///' Error halt.. ISTATE =',I3)
!     STOP

end program bde_background


     subroutine bde_equations(neq,t,y,ydot)
       double precision t, y, ydot
       dimension y(3), ydot(3)
       double precision ws

       double precision grhog,grhornomass,grhob,grhoc
       double precision alphaBDE,ac

       grhog = 8.254089497837408E-012
       grhornomass = 5.709922326520952E-012
       grhob = 7.517000357645051E-009
       grhoc = 3.917061156685587E-008
       alphaBDE = 0.666666666666666
       ac = 2.484836566392821E-006
                
       ws = (1./3.)*(grhog+grhornomass)/(ac*exp(t)*(grhob+grhoc)+(grhog+grhornomass))
       
       
       ydot(1)=-3.*y(1)+sqrt(3./2.)*y(3)*y(2)**2+(3./2.)*y(1)*(2*y(1)**2+(1.+ws)*(1-y(1)**2-y(2)**2))
       ydot(2)=-sqrt(3./2.)*y(1)*y(2)*y(3)+(3./2.)*y(2)*(2*y(1)**2+(1.+ws)*(1-y(1)**2-y(2)**2))
       ydot(3)=-sqrt(6.)*y(1)*y(3)**2/alphaBDE   

       return
     end subroutine


      subroutine jac(neq,t,y,ml,mu,pd,nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(neq), pd(nrowpd,neq)

       double precision grhog,grhornomass,grhob,grhoc
       double precision alphaBDE,ac

       grhog = 8.254089497837408E-012
       grhornomass = 5.709922326520952E-012
       grhob = 7.517000357645051E-009
       grhoc = 3.917061156685587E-008
       alphaBDE = 0.666666666666666
       ac = 2.484836566392821E-006
                
       ws = (1./3.)*(grhog+grhornomass)/(ac*exp(t)*(grhob+grhoc)+(grhog+grhornomass))
     
      pd(1,1) = -3.+3./2.*(6.*y(1)*y(1)+(1.+ws)*(1.-3.*y(1)*y(1)-y(2)*y(2)))
      pd(1,2) = sqrt(6.)*y(3)*y(2)-3.*y(1)*y(2)*(1+ws)
      pd(1,3) = sqrt(6.)/2.*y(2)*y(2)
      
      pd(2,1) = -sqrt(6.)/2.*y(3)*y(2)+3./2.*y(2)*(4*y(1)-2.*y(1)*(1+ws))
      pd(2,2) = -sqrt(6.)/2.*y(3)*y(1)+3./2.*(2.*y(1)*y(1)+(1.+ws)*(1.-y(1)*y(1)-3.*y(2)*y(2)))
      pd(2,3) = -sqrt(6.)/2.*y(1)*y(2)
      
      pd(3,1) = -3.*sqrt(6.)/2.*y(3)*y(3)
      pd(3,2) = 0.
      pd(3,3) = -3.*sqrt(6.)*y(1)*y(3)
     
      return
      end subroutine

     
