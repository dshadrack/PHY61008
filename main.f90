program Sun_Earth_Jupiter
!Start Program
!Version 2
!Variable Declaration
Implicit none
!Define Variables
DoublePrecision :: r(0:7,0:2),v(0:7,0:2),F(0:7,0:2),a_0(0:7,0:2),a_1(0:7,0:2)
doubleprecision :: G, M(0:7),s(0:7,0:7), d(0:7,0:7,0:2), AU, E0,E1,Vabs_Squared(0:7),COM(0:2),Mtot,Cov(0:2)
doubleprecision :: P(0:7), Theta(1:6), pi, Msol,StarMass, TimeFactor, CurrentTime,Step, year, InitSep
Integer(8) :: i,j,k,t,done, tmax,length,x, n, Logging,StarDist, StarV
!Assign constant values
n = 8
year = 3600*24*365
AU = 1.496e11
G = 1*6.67e-11
pi = 4*atan(1.)
!Randomise Starting angles in xy
Call init_random_seed()
Call random_number(Theta)
Theta = Theta * 2 * pi
Msol = 1.99e30
Timefactor = 1.
!Indexing of Planets In the following order: Sun, Earth, Jupiter, Mars, Saturn, Uranus, Neptune, Test Star...
CurrentTime = 0
!Masses of bodies
M(0) = 1*Msol
M(1) = 1*5.97e24
M(2) = 1*1.898e27
M(3) = 1*6.4171e23
M(4) = 5.683e26
M(5) = 8.681e25
M(6) = 1*1.024e26
M(7) = 0

!Starting Velocities
v(0,0) = 0
v(0,1) = 0
v(0,2) = 0
v(1,0) = 29780 * -1. * cos(Theta(1))
v(1,1) = 29780 * -1. * sin(Theta(1))
v(1,2) = 0
v(2,0) = 13060 * -1. *cos(Theta(2))
v(2,1) = 13060 * -1. *sin(Theta(2))
v(2,2) = 0
v(3,0) = 24070 * -1. *cos(Theta(3))
v(3,1) = 24070 * -1. *sin(Theta(3))
v(3,2) = 0
v(4,0) = 9680 * -1. *cos(Theta(4))
v(4,1) = 9680 * -1. *sin(Theta(4))
v(4,2) = 0
v(5,0) = 6800 * -1. *cos(Theta(5))
v(5,1) = 6800 * -1. *sin(Theta(5))
v(5,2) = 0
v(6,0) = 5430 * -1. *cos(Theta(6))
v(6,1) = 5430 * -1. *sin(Theta(6))
v(6,2) = 0
v(7,0) = 0
v(7,1) = 0
v(7,2) = 0





r(0,0) = 0
r(0,1) = 0
r(0,2) = 0
r(1,0) = AU * sin(Theta(1)) * -1.
r(1,1) = AU * cos(Theta(1))
r(1,2) = 0
r(2,0) = 5.2*AU * sin(Theta(2)) * -1.
r(2,1) = 5.2*AU * cos(Theta(2))
r(2,2) = 0
r(3,0) = 1.52*AU * sin(Theta(3)) * -1.
r(3,1) = 1.52*AU * cos(Theta(3))
r(3,2) = 0
r(4,0) = 9.583*AU * sin(Theta(4)) * -1.
r(4,1) = 9.583*AU * cos(Theta(4))
r(4,2) = 0
r(5,0) = 19.201*AU * sin(Theta(5)) * -1.
r(5,1) = 19.201*AU * cos(Theta(5))
r(5,2) = 0
r(6,0) = 30.07 * AU * sin(Theta(6)) * -1.
r(6,1) = 30.07 * AU * cos(Theta(6))
r(6,2) = 0
r(7,0) = 0
r(7,1) = -700 * AU
r(7,2) = 0


Write(6,*) ("Enter Encounter Star Mass")
Read *,StarMass
M(7) = StarMass * Msol
Write(6,*) ("Enter Encounter Star Perpendicular Distance")
Read *,StarDist
r(7,0) = Stardist * AU
Write(6,*) ("Enter Encounter Star Velocity (m/s)")
Read *,StarV
v(7,1) = StarV * 1.



!Empty Arrays
F = 0
s = 0
d = 0
P = 0
Vabs_Squared = 0
COM = 0
Mtot = 0
Cov = 0





!Calculate centre of mass / centre of velocity information
Do i = 0,n-2
        Do k = 0,2
                COM(k) = COM(k) + r(i,k)*M(i)
                Cov(k) = Cov(k) + v(i,k)*M(i)
        End Do
Mtot = Mtot + M(i)
End Do

!Correct for COM and COV
Do i = 0,n-1
    r(i,:) = r(i,:) - COM/Mtot
    v(i,:) = v(i,:) - Cov/Mtot
    Do k = 0,2
        Vabs_Squared(i) = Vabs_Squared(i) + v(i,k)**2
    End Do
End Do


!Calculate Initial Conditions
Do i = 0,n-1
        !For each body
        Do j = 0,n-1
            !For each other body
            If (i /= j) then
            !Ignore self
                !Calculate separation from vectors
                s(i,j) = ((r(i,0)-r(j,0))**2 + (r(i,1) -r(j,1))**2 + (r(i,2) - r(j,2))**2)**0.5
                !Running total GPE
                P(i) = P(i) - G*M(i)*M(j)/s(i,j)
                    Do k = 0,2
                    !Find vectors i -> j
                    d(i,j,k) = r(j,k) - r(i,k)
                    !Calculate forces in dimension k, add to force vector array
                    F(i,k) = F(i,k) + (G*M(i)*M(j)/(s(i,j))**2)*d(i,j,k)/s(i,j)
                    End Do
            End If
        End Do
End Do
P = P/2
!Calculate Initial Acceleration
Do i = 0,n-1
    a_0(i,:) = F(i,:)/M(i)
End Do

!State initial conditions of the simulation
Write(6,*) "Initial Conditions"
Write(6,*) "Position Vectors xyz (m)"
Write(6,*) ""
!Print Initial coordinates of all bodies
Do i = 0,n-1
    Write(6,*) r(i,:)
End Do
Write(6,*) ""
Write(6,*) "Velocities xyz (m/s)"
Write(6,*) ""
!Print initial velocities of all bodies
Do i = 0,n-1
    Write(6,*) v(i,:)
End Do
Write(6,*) ""

!Report initial energy of system
write(6,*) "Initial Energy of System (J)"
Write(6,*) ""
E0 = 0
Do i = 0,n-1
    E0 = E0 + 0.5*M(i)*Vabs_Squared(i) + P(i)
End Do
Write(6,*) E0
Write(6,*) ""

!Ask for user defined parameters

Write(6,*) ""
Write(6,*) "Choose Timestep factor"
read *, TimeFactor
Write(6,*) ""
Write(6,*) "Log Data to Files? (Type 1 for yes, 0 for no)"
read *, Logging
Write(6,*) ""
step = (1/ (maxval(a_0**2))**0.25) * TimeFactor
Write(6,*) Step

!Open Log Files
If (Logging == 1) then
    open(2,file = 'Sun Motion.csv')
    open(3,file = 'Earth Motion.csv')
    open(4,file = 'Jupiter Motion.csv')
    open(7,file = 'Mars Motion.csv')
    open(8,file = 'Saturn Motion.csv')
    open(9,file = 'Uranus Motion.csv')
    open(10,file = 'Neptune Motion.csv')
    open(11, file = 'star motion.csv')
End IF
!Begin Simulation

!Determine initial separation between Sun and test star
InitSep = s(0,7)
!Run until test star has passed as far away from the Sun as it was initially
Do while (s(0,7) <= InitSep)

    F = 0
    s = 0
    d = 0
    !Calculate dynamic timestep based on max acceleration
    step = (1/ (maxval(a_0**2))**0.25) * TimeFactor
    !Find new r pos after time-step
    r = r + v*step + 0.5*a_0*step**2
    !For each body
    !!$OMP PARALLEL
    !!$OMP DO PRIVATE(i,j,k)
    Do i = 0,n-1
        !For each other body
        Do j = 0,n-1
        If (i /= j) then
            !Find absolute distance between bodies i and j
            s(i,j) = ((r(i,0)-r(j,0))**2 + (r(i,1) -r(j,1))**2 + (r(i,2) - r(j,2))**2)**0.5
            !For each dimension xyz
            Do k = 0,2
                !Find vectors i -> j
                d(i,j,k) = r(j,k) - r(i,k)
                !Calculate forces in dimension k, add to force vector array
                F(i,k) = F(i,k) + (G*M(i)*M(j)/(s(i,j))**2)*d(i,j,k)/s(i,j)
            End Do
        End if
        End Do
    End Do
    !!$OMP END DO
    !!$OMP END PARALLEL
    !Find New Accelerations
    Do i = 0,n-1
        a_1(i,:) = F(i,:)/M(i)
    End do



    !Find new velocities after time-step
    v = v + 0.5*(a_0 + a_1)*step
    a_0 = a_1

    !Write Every 250th step to log files
    If ((Mod(t,500) == 0) .and. (Logging == 1)) then
        Write(2,*) CurrentTime,',', d(0,0,0),',',d(0,0,1),',',d(0,0,2),',', s(0,0)
        Write(3,*) CurrentTime,',', d(0,1,0),',',d(0,1,1),',',d(0,1,2),',', s(0,1)
        Write(4,*) CurrentTime,',', d(0,2,0),',',d(0,2,1),',',d(0,2,2),',', s(0,2)
        Write(7,*) CurrentTime,',', d(0,3,0),',',d(0,3,1),',',d(0,3,2),',', s(0,3)
        Write(8,*) CurrentTime,',', d(0,4,0),',',d(0,4,1),',',d(0,4,2),',', s(0,4)
        Write(9,*) CurrentTime,',', d(0,5,0),',',d(0,5,1),',',d(0,5,2),',', s(0,5)
        Write(10,*) CurrentTime,',', d(0,6,0),',',d(0,6,1),',',d(0,6,2),',', s(0,6)
        Write(11,*) CurrentTime ,',', s(0,7)
    End If
!Keep track of how much time has elapsed
CurrentTime = CurrentTime + Step
t = t + 1
End Do
!Close Log files
If (Logging == 1) then
    close(2)
    close(3)
    close(4)
    close(7)
    close(8)
    close(9)
    close(10)
    close(11)
End if

!State Final Conditions
Write(6,*) ""
Write(6,*) "Final Conditions"
Write(6,*) ""
Write(6,*) "Time Elapsed", CurrentTime / year, "Years"
Write(6,*) ""
Write(6,*) "Absolute distance to Sun (AU)"
Write(6,*) s(0,0)/AU
Write(6,*) s(0,1)/AU
Write(6,*) s(0,2)/AU
Write(6,*) s(0,3)/AU
Write(6,*) s(0,4)/AU
Write(6,*) s(0,5)/AU
Write(6,*) s(0,6)/AU
Write(6,*) s(0,7)/AU


!Report Final Energy of System / Compare to Initial

!Loop through each body to find GPE of each
P = 0
Vabs_Squared = 0
!For each body
Do i = 0,n-1
        !For each other body
        Do j = 0,n-1
        !Ignore self
        If (i /= j) then
            !Find absolute separation from position vectors
            s(i,j) = ((r(i,0)-r(j,0))**2 + (r(i,1) -r(j,1))**2 + (r(i,2) - r(j,2))**2)**0.5
            !Sum calculated GPE
            P(i) = P(i) - G*M(i)*M(j)/s(i,j)
        End if
        End Do
        !Sum kinetic energy in each dimension (k)
        Do k = 0,2
            Vabs_Squared(i) = Vabs_Squared(i) + v(i,k)**2
        End Do
End Do
!Half to account for doubling up of counting in loop
P = P/2
E1 = 0
!Sum total energy
Do i = 0,n-1
    E1 = E1 +   0.5*M(i)*Vabs_Squared(i) + P(i)
End Do

!Report final energy and accuracy test
Write(6,*) ""
write(6,*) "Final Energy of System (J)"
write(6,*) E1
Write(6,*) "Percentage of Initial Energy Retained (Indicative of Sim Accuracy)"
write(6,*) E1/E0 * 100

!Stops program from closing automatically
read *, done




!Subroutine for seeding the random number generator
!From https://gcc.gnu.org/onlinedocs/gcc-4.6.1/gfortran/RANDOM_005fSEED.html
Contains
subroutine init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
end
end program Sun_Earth_Jupiter

