PROGRAM exercise1
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

!===========================================
!|||||||||||||  DATA VARIABLES |||||||||||||
!===========================================

REAL*8 :: temp,temp_cur, sqdist_cur, temp_avg, temp_avg_var
REAL*8 :: running_press_p
REAL*8:: running_press_vir, cur_press_vir, cur_press_p, cur_press_id
REAL*8:: press_vir_avg, press_vir_avg_var
REAL*8, DIMENSION(iteration_len/avg_len):: temp_dat, press_vir_dat, press_p_dat, press_id_dat
REAL*8, DIMENSION(800)::pair_hist,pair_hist_ideal
REAL*8, DIMENSION(iteration_len)::sqdist_dat
REAL*8, DIMENSION(eq_len/eq_len_avg)::temp_eq_dat
REAL*8, DIMENSION(eq_len*6+700*6+iteration_len):: temp_cur_dat

!====================================
!|||||| Pre-program initialization|||
!====================================

  CALL Inputs
  
  !Set up a few variables
  min_sep =1.122462048d0*(1/density)**(.333)
  box_size=cubes*1.5874010519681996d0*(1/density)**(.333)
  maxD = (box_size)
 


  CALL Initialize  
  
  CALL Update_Forces

  !pair_hist_ideal = pair_hist*1000
 
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
pair_hist = 0
pos_init=pos
running_press_vir=0

!==================================
!|||||| Main Loop||||||||||||||||||
!==================================

DO t = 1, iteration_len

  ! Evaluate Equations of motion
  CALL MoveParticles
  !CALL PlotParticles

  !Update Quantities of Interest
  Kin = Kinetic(vel)
  temp_cur = Temperature(Kin)
  energy = Kin + Potij
  cur_press_id = press_id(temp_cur)
  sqdist_cur = sum((pos_unbound-pos_init)**2.0d0)/DFLOAT(N)
  cur_press_vir= (running_press_vir+N*temp_cur)/(box_size**3.0d0)


  !Write/Average appropriate quantities:
  
  !write(22,*) Kin/N, Potij/N, (Kin+Potij)/N
  !CALL average(temp_cur,temp_dat,t,avg_len)
  CALL average(cur_press_vir, press_vir_dat, t, avg_len)
  CALL average(cur_press_p, press_p_dat,t,avg_len)
  !CALL average(cur_press_id, press_id_dat,t,avg_len)
  !sqdist_dat(t) =sqdist_cur
  !temp_cur_dat(k) = temp_cur  
  !k=k+1
  !running_press_vir=0

END DO

!CALL EndPlot()

!===========================================
!|||||||| Post-Program calculations|||||||||
!===========================================
temp_avg = SUM(temp_dat)/SIZE(temp_dat)
temp_avg_var = var(temp_dat)

press_vir_avg = SUM(press_vir_dat)/SIZE(press_vir_dat)
press_vir_avg_var = var(press_vir_dat)

cur_press_p = running_press_p/(iteration_len)
!cur_press_id = N*temp_cur/(box_size**3.0d0)

!CALL Fix_pair_hist

!=================================
!||||||| WRITE OUTPUTS||||||||||||
!=================================


WRITE(23,fmt="(4(f12.6,1x))") press_vir_avg, press_vir_avg_var, temp_avg, temp_avg_var
!WRITE(24,fmt="(1(f12.6,1x))") press_id_dat
WRITE(25,fmt="(1(f12.6,1x))") press_p_dat
!WRITE(26,fmt="(1(f12.6,1x))") pair_hist
!WRITE(27,fmt="(1(f12.6,1x))") sqdist_dat
!WRITE(28,fmt="(1(f12.6,1x))") pair_hist_ideal
!WRITE(31,fmt="(1(f12.6,1x))") temp_cur_dat
!WRITE(33,fmt="(1(f12.6,1x))") temp_dat
!WRITE(32,*) temp_avg, temp_avg_var

PRINT*, "Done"


!=================================
!||||||||SUBROUTINES||||||||||||||
!==========================

CONTAINS


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

  !Open relevant ouput files
  !  open(unit=22,status='replace',file='energy.dat')  
  !  open(unit=23,status='old',file='press_vir.dat',position="Append") 
  !  open(unit=24,status='replace',file='press_id.dat') 
  !  open(unit=25,status='replace',file='press_p.dat')
  !  open(unit=26,status='replace',file='pair_dist2.dat') 
  !  open(unit=27,status='replace',file='sqdist.dat')
  !  open(unit=28,status='replace',file='pair_dist2_ideal.dat')
  !  OPEN(unit=31,status='replace',file='temp.dat')
  !  OPEN(unit=33,status='replace',file='temp_avg.dat')
  !  OPEN(unit=33,status='replace',file='temp_avg2.dat')

  !Initialize graphics
  !  CALL InitPlotF('lightblue', 800, 800, 'out.ps', 1)
  !  CALL PutStopButton()
  !  CALL Framing(-0.1D0*Box_Size, -0.1D0*Box_Size, 1.1D0*Box_Size, 1.1D0*Box_Size)

END SUBROUTINE Initialize


SUBROUTINE MoveParticles

  !Integrator (verlet)
  vel=vel+FTotal*T_step/2.d0
  Pos = Pos + Vel*T_step;
  Pos_unbound = Pos_unbound + Vel*T_step
  call update_forces
  vel=vel+FTotal*T_step/2.0d0
 
  !cur_press_p = press_p(vel,pos)
  running_press_p = running_press_p + press_p(vel,pos) !2.0d0*dot_product(vel(1,:),int(Pos(1,:)/(box_size)))

  !Boundary Conditions
  Pos = Modulo(Pos, box_size)

END SUBROUTINE MoveParticles


SUBROUTINE PlotParticles
  INTEGER :: I
  DO I=1, N
    CALL SetPoint(Pos(1,I), Pos(2,I))
  END DO
END SUBROUTINE PlotParticles

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
   
      !CALL update_pair_hist(pair_hist, Mag_Dij,maxD)      
      !CALL press_vir(Dij, Fij,running_press_vir)

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

subroutine average(observable, datastore, it_count, length_average) !takes observables Ai, divides it by lengt_average and stores it in array "datastore". 
 REAL*8 :: observable 
 REAL*8, DIMENSION(100) :: datastore
 INTEGER :: loc,it_count, length_average
  loc = (it_count-1)/length_average+1
  datastore(loc) = datastore(loc) + observable / length_average
end subroutine average

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

FUNCTION sq_dist(p_new,p_old)
  IMPLICIT NONE 
  REAL*8 :: sq_dist 
  REAL*8,DIMENSION(3):: p_new,p_old
  
  sq_dist=sum((P_new-P_old)**2.0d0)
END FUNCTION

SUBROUTINE Get_Momentum(Total_Momentum, Velocities)
  REAL*8, INTENT(OUT), DIMENSION(3) :: Total_Momentum
  REAL*8, INTENT(IN), DIMENSION(3,N) :: Velocities
  
  Total_Momentum(1) = sum(Velocities(1,:))
  Total_Momentum(2) = sum(Velocities(2,:))
  Total_Momentum(3) = sum(Velocities(3,:))
END SUBROUTINE !Get_Momentum

FUNCTION Potential(Distance_sq) ! Potential of two particles
  REAL*8 :: Potential
  REAL*8:: Distance_sq
  
  Potential = 4.d0/(Distance_sq**6.0d0) - 4.d0/(Distance_sq**3.0d0)
END FUNCTION !Potential

SUBROUTINE press_vir(Rij, Forceij, total)
  IMPLICIT NONE
  REAL*8,DIMENSION(:):: Rij, Forceij
  REAL*8::total

  total = total +DOT_PRODUCT( Rij,Forceij)/(3.0d0)
END SUBROUTINE

FUNCTION press_p(velocities,positions)
  IMPLICIT NONE
  REAL*8:: press_p
  REAL*8:: ptransfer
  REAL*8, DIMENSION(:,:) :: Velocities, positions

  ptransfer =2.0d0*sum(Velocities(1,:)*INT(positions(1,:)/box_size))
 
  press_p =ptransfer/(box_size**2.0d0*T_step)

END FUNCTION !press_p

FUNCTION press_id(temperature)
  IMPLICIT NONE
  REAL*8::press_id,temperature
 
  press_id = N*temperature/(box_size**3.0d0)
END FUNCTION !press_id

SUBROUTINE update_pair_hist(hist,Rij2, maxR)
  IMPLICIT NONE
  REAL*8, DIMENSION(:):: hist
  REAL*8:: maxR, Rij2
  INTEGER:: hist_len, loc

  hist_len = SIZE(hist)
  loc = CEILING((SQRT(Rij2)*hist_len/(maxR)))+1
 
 !print*, loc
  hist(loc) = hist(loc) + .000001d0

END SUBROUTINE !update_pair_hist  

SUBROUTINE fix_pair_hist
  IMPLICIT NONE
  bin_width = maxD/SIZE(Pair_hist)
  R=0

  DO t=1,SIZE(Pair_hist) 
    R = R+bin_width
    Vol = 4.0d0*Pi*bin_width*R**2.0d0
    Pair_hist(t)=Pair_hist(t)/Vol
    pair_hist_ideal(t)= pair_hist_ideal(t)/Vol 
  END DO 

END SUBROUTINE

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

end program exercise1
