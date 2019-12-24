program test

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

     EXTERNAL FEX
     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y
     DIMENSION Y(3), ATOL(3), RWORK(70), IWORK(23)
     NEQ = 3
     Y(1) = 1.
     Y(2) = 0.
     Y(3) = 0.
     T = 0.
     TOUT = .4
     ITOL = 2
     RTOL = 1.D-4
     ATOL(1) = 1.D-6
     ATOL(2) = 1.D-10
     ATOL(3) = 1.D-6
     ITASK = 1
     ISTATE = 1
     IOPT = 0
     LRW = 70
     LIW = 23
     JT = 2
     DO 40 IOUT = 1,12
       CALL DLSODA(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)
       WRITE(6,20)T,Y(1),Y(2),Y(3)
 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6)
       IF (ISTATE .LT. 0) GO TO 80
 40    TOUT = TOUT*10.
     WRITE(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15)
 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4/&
        ' Method last used =',I2,'   Last switch was at t =',D12.4)
     STOP
 80  WRITE(6,90)ISTATE
 90  FORMAT(///' Error halt.. ISTATE =',I3)
     STOP

end program test


SUBROUTINE FEX (NEQ, T, Y, YDOT)
     DOUBLE PRECISION T, Y, YDOT
     DIMENSION Y(3), YDOT(3)
     YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)
     YDOT(3) = 3.D7*Y(2)*Y(2)
     YDOT(2) = -YDOT(1) - YDOT(3)
     RETURN
END
