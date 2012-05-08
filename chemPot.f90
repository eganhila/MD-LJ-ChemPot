PROGRAM leonardJones
  ! A program to Compute the effects of the Lennard-Jones
  ! potential
  ! To compile: gfortran exercise1.f90 /home/delft03/XPS/lib/libxps.a /usr/lib/libX11.so

IMPLICIT NONE

!=========================================
!|||||||||| PARAMETERS ||||||||||||||||||
!=========================================
INTEGER, PARAMETER :: iteration_len= 2000,avg_len=20,cubes =4, N = 4*cubes**3 !Number of cubes in lattice and number of particles
INTEGER, PARAMETER :: eq_len=500, eq_len_avg=33
REAL*8, PARAMETER :: Pi = datan(1.0d0)*4.0d0

!==========================================
!||||||||||| VARIABLES ||||||||||||||||||||
!==========================================
REAL*8 :: density
REAL*8 :: bin_width,R
REAL*8 :: box_size,min_sep
INTEGER:: i, j,k,m,l,checks, t !loop variables
REAL*8 :: Mag_Dij, T_step = 0.002d0, time = 0, energy, Potij, Kin
REAL*8 :: maxD
REAL*8, DIMENSION(3,N):: Pos,Vel, Ftotal
REAL*8, DIMENSION(3,N):: pos_init, pos_unbound
REAL*8, DIMENSION(N):: rF=0.d0,rf_add
REAL*8, DIMENSION(3):: Fij, Dij, momentum_cur
REAL*8:: Vol


REAL*8 :: P_dot,pressure

!===========================================
!|||||||||||||  DATA VARIABLES |||||||||||||
!===========================================

REAL*8 :: temp,temp_cur, sqdist_cur, temp_avg, temp_avg_var
REAL*8, DIMENSION(iteration_len/avg_len) :: temp_dat
REAL*8, DIMENSION(eq_len/eq_len_avg)::temp_eq_dat
REAL*8, DIMENSION(eq_len*6+700*6+iteration_len):: temp_cur_dat

REAL*8, DIMENSION(iteration_len/avg_len) :: pressure_dat
!====================================
!|||||| Pre-program initialization|||
!====================================

  CALL Inputs
  
  !Set up a few variables
  min_sep =1.122462048d0*(1/density)**(.333)
  box_size=cubes*1.5874010519681996d0*(1/density)**(.333)
  maxD = (box_size)
  P_dot = 0

  CALL Initialize  
  CALL Update_Forces
 
!====================================
!||||| Temperature scaling|||||||||||
!====================================
k=1
do checks= 1,5          !Outer loop for temperature rescaling
  temp_eq_dat = 0
  DO l= 1,200           !First loop for rescaling before equilibration
    time=l*T_step
    Kin = Kinetic(Vel)
    temp_cur = Temperature(Kin)
    temp_cur_dat(k) = temp_cur  
    k=k+1
    CALL MoveParticles
  END DO

  DO l=1,eq_len         !Second loop for rescaling and finding temp average
    Kin= Kinetic(Vel)
    temp_cur = Temperature(Kin)
    temp_cur_dat(k) = temp_cur  
    k=k+1
    CALL average(temp_cur, temp_eq_dat,l,eq_len_avg)
    CALL MoveParticles
  END DO
    CALL rescaletemp(SUM(temp_eq_dat)/SIZE(temp_eq_dat))
    PRINT*, temp_cur
end do 

!Set running total variables to 0 because equilibration phase=garbage data


!==================================
!|||||| Main Loop||||||||||||||||||
!==================================

DO t = 1, iteration_len
  !print*, t
  P_dot=0

  ! Evaluate Equations of motion
  CALL MoveParticles
  !CALL PlotParticles

  !Update Quantities of Interest
  Kin = Kinetic(vel)
  temp_cur = Temperature(Kin)
  energy = Kin + Potij

  CALL getPressure()
  !Write/Average appropriate quantities:
  pressure_dat(t) = Pressure
END DO


!===========================================
!|||||||| Post-Program calculations|||||||||
!===========================================


!=================================
!||||||| WRITE OUTPUTS||||||||||||
!=================================

PRINT*,"Current Pressure: "
PRINT*, pressure
PRINT*,"Averaged Pressure: " 
PRINT*, blockAverage(pressure_dat)
PRINT*, "Done"


!=================================
!||||||||SUBROUTINES||||||||||||||
!=================================

CONTAINS

SUBROUTINE getPressure
REAL*8::Pdot,P_reduced,P
Pdot = sum( Dij*Fij,2)
P_reduced = 1 - 1/(3*N*Temp_cur)*Pdot
Pressure = P_reduced*temp_cur/density

END SUBROUTINE

SUBROUTINE inputs     !User submitted input parameters

  PRINT*, "Input simulation density:"
  READ*, density
  print*,
  
  print*, 'INPUT SIMULATION TEMPERATURE'
  read*, temp
  print*
END SUBROUTINE !Inputs


SUBROUTINE Initialize

  CALL random
  CALL init_positions(Pos, min_sep, N, cubes)
  
  CALL init_velocities(Vel,N,temp)

END SUBROUTINE Initialize


SUBROUTINE MoveParticles

  !Integrator (verlet)
  vel=vel+FTotal*T_step/2.d0
  Pos = Pos + Vel*T_step;
  Pos_unbound = Pos_unbound + Vel*T_step
  call update_forces
  vel=vel+FTotal*T_step/2.0d0
 
  !Boundary Conditions
  Pos = Modulo(Pos, box_size)

END SUBROUTINE MoveParticles


!SUBROUTINE PlotParticles
!  INTEGER :: I
!  DO I=1, N
!    CALL SetPoint(Pos(1,I), Pos(2,I))
!  END DO
!END SUBROUTINE PlotParticles

SUBROUTINE Update_Forces
  Ftotal = 0.d0
  Potij = 0.d0
  Do i = 1,N-1
    Do j = i+1,N
      
      Dij = Pos(:,i) - Pos(:,j)
      Dij = Dij - NINT(dij/box_size)*box_size
      
      Mag_Dij = DOT_PRODUCT(Dij, Dij) 
      Fij = (48.d0/(Mag_Dij**7.0d0)-24.d0/(Mag_Dij**4.0d0))*Dij
        
      Ftotal(:,i) = Ftotal(:,i) + Fij;
      Ftotal(:,j) = Ftotal(:,j) - Fij;
      
      Potij = Potij + Potential(Mag_Dij);
      
      P_dot = P_dot + DOT_PRODUCT(Dij,Fij)

    END DO
  END DO
END SUBROUTINE Update_Forces


subroutine init_positions( p_array, sep0 , N ,cubes)
  !Initialize all positions
  implicit none 
  
  integer i,j,k,N,ll,cubes,ii,cnt
  real*8,dimension(3,N) :: p_array
  real*8,dimension(3,4) :: box
  real*8 sep0,sep

  sep=sep0/dsqrt(2.d0)
  
  
  !small box
  
  box=0.d0
  box(1,2)=sep;box(3,2)=sep
  box(2,3)=sep;box(3,3)=sep
  box(1,4)=sep;box(2,4)=sep

  ll=1
  do i=0,cubes-1
     do j=0,cubes-1
        do k=0,cubes-1
           cnt=1
           do ii=ll,ll+3
            p_array(1,ii)=box(1,cnt)+k*2.d0*sep
            p_array(2,ii)=box(2,cnt)+j*2.d0*sep
            p_array(3,ii)=box(3,cnt)+i*2.d0*sep
            cnt=cnt+1
           end do 
           ll=ll+4
        end do 
     end do 
  end do 
end subroutine init_positions

subroutine init_velocities( v_array, n,temp)
  !Use gaussian maxwell boltzman distribution to assign initial velocities
 implicit none
 
 integer n,i,j
 real*8 x,y,z0,z1,ch,temp
 real*8,dimension(3,n) :: v_array
 !real*8:: Pi
   
 do i=1,n
  do j=1,3
 call random_number(x)
 call random_number(y)
 
 z0=dSqrt(-2.0d0*dlog(x))*dcos(2.0d0*Pi*y)
 z1=dsqrt(-2.0d0*dlog(x))*dsin(2.0d0*Pi*y)
 
 call random_number(ch)
 if (ch > .5) then 
  v_array(j,i)=dSqrt(temp)*z0
 else
  v_array(j,i)=dsqrt(temp)*z1
 end if 
 
 end do 
 end do 
 
 do i=1,3
  v_array(i,:)=v_array(i,:)-sum(v_array(i,:))/dfloat(n)
 end do   

end subroutine !init_velocities
 
subroutine random !Random number generation
 implicit none 
 
 integer seed,i
 real*8 x
 
 open(unit=26,file='seed')
  read(26,*) seed
 close(26)
 
 do i=1,seed
  call random_number(x)
 end do 
 
 open(unit=26,status='replace',file='seed')
  write(26,*) ceiling(1000*x)
 close(26)
 
end subroutine random

subroutine rescaletemp(avg_temp)
  REAL*8::avg_temp
  
  vel = (Temp/avg_temp)**(0.5d0)*vel;
end subroutine

FUNCTION Temperature(Kin_Energy)
  IMPLICIT NONE
  REAL*8 :: Temperature
  REAL*8, INTENT(IN) :: Kin_Energy
  
  Temperature = 2.d0*Kin_Energy/(3.d0*dfloat(N-1))
END FUNCTION !Temperature

FUNCTION Kinetic(Velocities)
  REAL*8 :: Kinetic
  REAL*8, INTENT(IN), DIMENSION(3,N) :: Velocities
  
  Kinetic = sum(Velocities**2.0d0)/2.d0
END FUNCTION !Kinetic

FUNCTION Potential(Distance_sq) ! Potential of two particles
  REAL*8 :: Potential
  REAL*8:: Distance_sq
  
  Potential = 4.d0/(Distance_sq**6.0d0) - 4.d0/(Distance_sq**3.0d0)
END FUNCTION !Potential

FUNCTION blockAverage(dat)
  IMPLICIT NONE
  REAL*8, DIMENSION(:):: dat
  REAL*8 :: blockAverage
  INTEGER :: corTime, numBlocks
  REAL*8 :: avg
  INTEGER :: i,j,k
  
  corTime = 20
  numBlocks = SIZE(dat)/corTime
  k = 0

  DO i = 1,numBlocks
    DO j=1,corTime
      k=k+1
      avg = avg+dat(k)/corTime
    END DO
    blockAverage = blockAverage + avg/numBlocks
    avg = 0
  END DO

END FUNCTION

FUNCTION var(dat)
  IMPLICIT NONE
  REAL*8, DIMENSION(:):: dat
  REAL*8 :: avg, var
  INTEGER :: iter

  var=0
  avg= SUM(dat)/SIZE(dat)
  
  !print*, dat(1), avg, dat(1)-avg
  DO iter = 1, SIZE(dat)
    
    var = var+(dat(iter)-avg)**2.0d0
    !print*, var
  END DO

  var=var/DFLOAT(SIZE(dat))

END FUNCTION !var

end program leonardJones
