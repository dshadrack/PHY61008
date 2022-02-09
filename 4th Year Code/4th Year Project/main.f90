
program NBody
!Start Program
!Version 2
!Variable Declaration
Implicit none
!Define Variables
doubleprecision :: r(0:6,0:2,0:2),v(0:6,0:2,-6:2),a(0:6,0:2,-6:2), ai(0:6,0:2,0:1), vi(0:6,0:2,0:1),pcu, munit,runit,tunit,vunit
doubleprecision :: G, M(0:6),s(0:6,0:6), d(0:6,0:6,0:2), AU, E0,E1,Vabs_Squared(0:6),COM(0:2),Mtot,Cov(0:2), Aerr, Verr
doubleprecision :: P(0:6), Theta(1:6), pi, Msol, TimeFactor, CurrentTime,Step, year, Small, relerr,length,PCerror,start,finish
Integer(8) :: i,j,k,t,done, n, Logging, counter, frac,systime
!Assign constant values
Do
n = 7
t = 0
year = 3600.*24.*365.25
AU = 1.496e11
G = 1*6.67408e-11
pi = 4.*atan(1.)
counter = 0
Small = 1.0e-10
relerr = 5.0e-12
frac = 0
!
!pcu=3.086d18
!munit = 1.99e33
!runit = AU
!tunit=SQRT((runit*pcu**3)/(G*munit*Msol))/year
!vunit=runit*pcu/(tunit*year*1.d5)
!write(6,*) 'tunit/yrs', tunit
!write(6,*) 'vunit/km/s', vunit

!Randomise Starting angles in xy
Call init_random_seed()
Call random_number(Theta)
Theta = Theta * 2 * pi
Msol = 1.99e30
Timefactor = 1.
!Indexing of Planets In the following order: Sun, Earth, Jupiter, Mars, Saturn, Uranus, Neptune
CurrentTime = 0.
!Masses of bodies
M(0) = 1.*Msol
M(1) = 1.*5.97e24
M(2) = 1.*1.898e27
M(3) = 1.*6.4171e23
M(4) = 5.683e26
M(5) = 8.681e25
M(6) = 1.*1.024e26

v = 0.
!Starting Velocities
v(0,0,0) = 0.
v(0,1,0) = 0.
v(0,2,0) = 0.
v(1,0,0) = 29780. * -1. * cos(Theta(1))
v(1,1,0) = 29780. * -1. * sin(Theta(1))
v(1,2,0) = 0.
v(2,0,0) = 13060. * -1. *cos(Theta(2))
v(2,1,0) = 13060. * -1. *sin(Theta(2))
v(2,2,0) = 0.
v(3,0,0) = 24070. * -1. *cos(Theta(3))
v(3,1,0) = 24070. * -1. *sin(Theta(3))
v(3,2,0) = 0.
v(4,0,0) = 9680. * -1. *cos(Theta(4))
v(4,1,0) = 9680. * -1. *sin(Theta(4))
v(4,2,0) = 0.
v(5,0,0) = 6800. * -1. *cos(Theta(5))
v(5,1,0) = 6800. * -1. *sin(Theta(5))
v(5,2,0) = 0.
v(6,0,0) = 5430. * -1. *cos(Theta(6))
v(6,1,0) = 5430. * -1. *sin(Theta(6))
v(6,2,0) = 0.

!Starting positions
r = 0.
r(0,0,0) = 0.
r(0,1,0) = 0.
r(0,2,0) = 0.
r(1,0,0) = AU * sin(Theta(1)) * -1.
r(1,1,0) = AU * cos(Theta(1))
r(1,2,0) = 0.
r(2,0,0) = 5.2*AU * sin(Theta(2)) * -1.
r(2,1,0) = 5.2*AU * cos(Theta(2))
r(2,2,0) = 0.
r(3,0,0) = 1.52*AU * sin(Theta(3)) * -1.
r(3,1,0) = 1.52*AU * cos(Theta(3))
r(3,2,0) = 0.
r(4,0,0) = 9.583*AU * sin(Theta(4)) * -1.
r(4,1,0) = 9.583*AU * cos(Theta(4))
r(4,2,0) = 0.
r(5,0,0) = 19.201*AU * sin(Theta(5)) * -1.
r(5,1,0) = 19.201*AU * cos(Theta(5))
r(5,2,0) = 0.
r(6,0,0) = 30.07 * AU * sin(Theta(6)) * -1.
r(6,1,0) = 30.07 * AU * cos(Theta(6))
r(6,2,0) = 0.

!Ask for length of time sim should run for
Write(6,*) ("Enter Sim Duration (Years)")
Read *, length



!Empty Arrays
a = 0.
s = 0.
d = 0.
P = 0.
Vabs_Squared = 0.
COM = 0.
Mtot = 0.
Cov = 0.





!Calculate centre of mass / centre of velocity information
Do i = 0,n-1
        Do k = 0,2
                COM(k) = COM(k) + r(i,k,0)*M(i)
                Cov(k) = Cov(k) + v(i,k,0)*M(i)
        End Do
Mtot = Mtot + M(i)
End Do

!Correct for COM and COV
Do i = 0,n-1
    r(i,:,0) = r(i,:,0) - COM/Mtot
    v(i,:,0) = v(i,:,0) - Cov/Mtot
    Do k = 0,2
        Vabs_Squared(i) = Vabs_Squared(i) + v(i,k,0)**2
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
                    s(i,j) = ((r(i,0,0)-r(j,0,0))**2 + (r(i,1,0) -r(j,1,0))**2 + (r(i,2,0) - r(j,2,0))**2)**0.5
                    !Running total GPE
                    P(i) = P(i) - G*M(i)*M(j)/s(i,j)
                    !Find vectors i -> j
                    d(i,j,:) = r(j,:,0) - r(i,:,0)
                    !Calculate forces in dimension k, add to force vector array
                    a(i,:,0) = a(i,:,0) + (G*M(j)/(s(i,j))**2)*d(i,j,:)/s(i,j)
                End If
            End Do
    End Do
P = P/2.


Write(6,*) ""

!Report initial energy of system
write(6,*) "Initial Energy of System (J):"
E0 = 0
Do i = 0,n-1
    E0 = E0 + 0.5*M(i)*Vabs_Squared(i) + P(i)
End Do
Write(6,*) E0
!Small timestep while bootstrapping
Step = 1
!Bootstrap into ABM using 2nd order taylor expansion
Do k = 0,3
    !Find new r pos after time-step
    r(:,:,1) = r(:,:,0) + v(:,:,0)*step + 0.5*a(:,:,0)*step**2

    !Reset arrays for current loop
    s = 0.
    d = 0.
    a(:,:,1) = 0.

    !Accleration calculation loop
    Do i = 0,n-1
        !For each other body
        Do j = 0,n-1
            If (i /= j) then
                !Find absolute distance between bodies i and j
                s(i,j) = ((r(i,0,1)-r(j,0,1))**2 + (r(i,1,1) -r(j,1,1))**2 + (r(i,2,1) - r(j,2,1))**2)**0.5
                !For each dimension xyz
                !Find vectors i -> j
                d(i,j,:) = r(j,:,1) - r(i,:,1)
                !Calculate field strength in dimension, add to acceleration vector array
                a(i,:,1) = a(i,:,1) + (G*M(j)/(s(i,j))**2)*d(i,j,:)/s(i,j)
            End if
        End Do
    End Do

    !Update velocity using newfound acceleration
    v(:,:,1) = v(:,:,0) + 0.5*(a(:,:,0) + a(:,:,1))*Step

    !Shift all datapoints backwards in time one step
    a(:,:,-6:0) = a(:,:,-5:1)
    v(:,:,-6:0) = v(:,:,-5:1)
    !Keep track of time elapsed
    t = t+1
    CurrentTime = CurrentTime + Step
    counter = counter + 1

End Do

!Ask for user defined parameters

Write(6,*) ""
Write(6,*) "Log Data to Files? (Type 1 for yes, 0 for no)"
read *, Logging
Write(6,*) ""

!Open Log Files
If (Logging == 1) then
    open(2,file = 'Sun Motion.csv')
    open(3,file = 'Earth Motion.csv')
    open(4,file = 'Jupiter Motion.csv')
    open(7,file = 'Mars Motion.csv')
    open(8,file = 'Saturn Motion.csv')
    open(9,file = 'Uranus Motion.csv')
    open(10,file = 'Neptune Motion.csv')
End IF


!Begin Simulation
systime = Time()
Call cpu_time(start)
Do while ((CurrentTime / year) < length)
    s = 0.
    d = 0.

    !Predict r pos using ABM predictor
    r(:,:,1) = r(:,:,0) + (Step/24.) * (-9.* v(:,:,-3) +37.* v(:,:,-2) -59. * v(:,:,-1) +55. * v(:,:,0))
    v(:,:,1) = v(:,:,0) + (Step/24.) * (-9.* a(:,:,-3) +37.* a(:,:,-2) -59. * a(:,:,-1) +55. * a(:,:,0))
    !Reset predicted acceleration values for current calculation
    a(:,:,1) = 0.

    !Acceleration Loop
    !For each body
    Do i = 0,n-1
        !For each other body
        Do j = 0,n-1
            If (i /= j) then
                !Find absolute distance between bodies i and j
                s(i,j) = ((r(i,0,1)-r(j,0,1))**2 + (r(i,1,1) -r(j,1,1))**2 + (r(i,2,1) - r(j,2,1))**2)**0.5
                !For each dimension xyz
                !Find vectors i -> j
                d(i,j,:) = r(j,:,1) - r(i,:,1)
                !Calculate field strength in dimension k, add to acceleration vector array
                a(i,:,1) = a(i,:,1) + (G*M(j)/(s(i,j))**2)*d(i,j,:)/s(i,j)
            End if
        End Do
    End Do

    !Correct using predicted values using ABM corrector
    r(:,:,2) = r(:,:,0) + (Step/24.) * ( v(:,:,-2) - 5.* v(:,:,-1) + 19.*v(:,:,0) + 9.* v(:,:,1) )
    v(:,:,2) = v(:,:,0) + (Step/24.) * ( a(:,:,-2) - 5.* a(:,:,-1) + 19.*a(:,:,0) + 9.* a(:,:,1) )

    a(:,:,2) = 0
    !Corrected Acceleration Loop
    Do i = 0,n-1
        !For each other body
        Do j = 0,n-1
            If (i /= j) then
                !Find absolute distance between bodies i and j
                s(i,j) = ((r(i,0,2)-r(j,0,2))**2 + (r(i,1,2) -r(j,1,2))**2 + (r(i,2,2) - r(j,2,2))**2)**0.5
                !For each dimension xyz
                !Find vectors i -> j
                d(i,j,:) = r(j,:,2) - r(i,:,2)
                !Calculate field strength in dimension k, add to acceleration vector array
                a(i,:,2) = a(i,:,2) + (G*M(j)/(s(i,j))**2)*d(i,j,:)/s(i,j)
            End if
        End Do
    End Do


    !Shift all datapoints one step backwards
    v(:,:,-6:0) = v(:,:,-5:1)
    a(:,:,-6:0) = a(:,:,-5:1)
    r(:,:,0) = r(:,:,2)
    v(:,:,0) = v(:,:,2)
    a(:,:,0) = a(:,:,2)


    !Find the largest discrepancy between predictor and corrector
    PCerror = (19./270.) * (maxval(  abs( r(:,:,2) - r(:,:,1) )/(abs(r(:,:,2))+ Small )  ))
    Aerr = (19./270.) * (maxval(  abs( a(:,:,2) - a(:,:,1) )/(abs(a(:,:,2))+ Small )  ))
    Verr = (19./270.) * (maxval(  abs( v(:,:,2) - v(:,:,1) )/(abs(v(:,:,2))+ Small )  ))

    PCerror = max(Aerr,PCerror)
    PCerror = max(Verr,PCerror)


    !Check if predictor and corrector agree to overly precise level
    If ((PCerror <(relerr * 0.01))) then
        !Ensure there are enough previous points in memory
        If ((counter >= 7)) then

            !Omit every second historic point
            a(:,:,-1) = a(:,:,-2)
            a(:,:,-2) = a(:,:,-4)
            a(:,:,-3) = a(:,:,-6)

            v(:,:,-1) = v(:,:,-2)
            v(:,:,-2) = v(:,:,-4)
            v(:,:,-3) = v(:,:,-6)

            !Double the timestep and reset the counter
            Step = Step * 2.
            counter = 0
            !write(6,*) "Growing dt to ", Step, t, (maxval(abs(r(:,:,2) - r(:,:,1))/ (abs(r(:,:,2)) + Small)) )
        End If
    else

       If (PCerror > relerr) then
            If (counter >= 4) then
                ai(:,:,0) = (1./128.) * (-5.* a(:,:,-4) + 28.* a(:,:,-3) - 70.* a(:,:,-2) + 140.* a(:,:,-1) + 35.* a(:,:,0))
                ai(:,:,1) = (1./128.) * (3.* a(:,:,-4) - 20.* a(:,:,-3) + 90.* a(:,:,-2) + 60.* a(:,:,-1) - 5.* a(:,:,0))

                vi(:,:,0) = (1./128.) * (-5.* v(:,:,-4) + 28.* v(:,:,-3) - 70.* v(:,:,-2) + 140.* v(:,:,-1) + 35.* v(:,:,0))
                vi(:,:,1) = (1./128.) * (3.* v(:,:,-4) - 20.* v(:,:,-3) + 90.* v(:,:,-2) + 60.* v(:,:,-1) - 5.* v(:,:,0))

                a(:,:,-2) = a(:,:,-1)
                a(:,:,-1) = ai(:,:,0)
                a(:,:,-3) = ai(:,:,1)

                v(:,:,-2) = v(:,:,-1)
                v(:,:,-1) = vi(:,:,0)
                v(:,:,-3) = vi(:,:,1)


                Step = Step * 0.5
                counter = 0
                !write(6,*) "Shrinking dt to", Step,t, (maxval(abs(r(:,:,2) - r(:,:,1))/ (abs(r(:,:,2)) + Small)))
            end if
        End If
    End If
!
    counter = counter + 1

    !Write Every 2000th step to log files
    If ((Mod(t,2000) == 0) .and. (Logging == 1)) then

        P = 0.
        Vabs_Squared = 0.
        !For each body
        Do i = 0,n-1
                !For each other body
                Do j = 0,n-1
                !Ignore self
                If (i /= j) then
                    !Find absolute separation from position vectors
                    s(i,j) = ((r(i,0,0)-r(j,0,0))**2 + (r(i,1,0) -r(j,1,0))**2 + (r(i,2,0) - r(j,2,0))**2)**0.5
                    !Sum calculated GPE
                    P(i) = P(i) - G*M(i)*M(j)/s(i,j)
                End if
                End Do
                !Sum kinetic energy in each dimension (k)
                Do k = 0,2
                    Vabs_Squared(i) = Vabs_Squared(i) + v(i,k,0)**2
                End Do
        End Do
        !Half to account for doubling up of counting in loop
        P = P/2.
        E1 = 0.
        !Sum total energy
        Do i = 0,n-1
            E1 = E1 +   0.5*M(i)*Vabs_Squared(i) + P(i)
        End Do

        Write(2,*) CurrentTime,',', d(0,0,0),',',d(0,0,1),',',d(0,0,2),',', s(0,0),',', E1/E0
        Write(3,*) CurrentTime,',', d(0,1,0),',',d(0,1,1),',',d(0,1,2),',', s(0,1)
        Write(4,*) CurrentTime,',', d(0,2,0),',',d(0,2,1),',',d(0,2,2),',', s(0,2)
        Write(7,*) CurrentTime,',', d(0,3,0),',',d(0,3,1),',',d(0,3,2),',', s(0,3)
        Write(8,*) CurrentTime,',', d(0,4,0),',',d(0,4,1),',',d(0,4,2),',', s(0,4)
        Write(9,*) CurrentTime,',', d(0,5,0),',',d(0,5,1),',',d(0,5,2),',', s(0,5)
        Write(10,*) CurrentTime,',', d(0,6,0),',',d(0,6,1),',',d(0,6,2),',', s(0,6)


    End If

!Keep track of how much time has elapsed
CurrentTime = CurrentTime + Step
t = t + 1

!Mechanism for printing out a progress update in 20% increments
!Can probably be done in a better way
If (CurrentTime / (length * year) > 0.2) Then
    If (CurrentTime / (length * year) > 0.4) Then
        If (CurrentTime / (length * year) > 0.6) then
            If (CurrentTime / (length * year) > 0.8) then
                If (frac == 6) then
                    Write(6,*) "80% Done"
                    frac = 8
                End If
            else
            If (frac == 4) then
                    Write(6,*) "60% Done"
                    frac = 6
                End If
            End If
        else
            If (frac == 2) then
                    Write(6,*) "40% Done"
                    frac = 4
            End If
        End If
    else
        If (frac == 0) then
            Write(6,*) "20% Done"
            frac = 2
        End If
    End If
End If
End Do
Write(6,*) "100% Done!"
!Close Log files
If (Logging == 1) then
    close(2)
    close(3)
    close(4)
    close(7)
    close(8)
    close(9)
    close(10)
End if

!State Final Conditions
Write(6,*) ""
Write(6,*) "Final Conditions"
Write(6,*) ""
Write(6,*) "Time Elapsed", CurrentTime / year, "Years"
Write(6,*) ""
Write(6,*) "Absolute distance to Sun (AU)"
Write(6,*) "Sun",s(0,0)/AU
Write(6,*) "Earth",s(0,1)/AU
Write(6,*) "Mars",s(0,3)/AU
Write(6,*) "Jupiter",s(0,2)/AU
Write(6,*) "Saturn",s(0,4)/AU
Write(6,*) "Uranus",s(0,5)/AU
Write(6,*) "Neptune",s(0,6)/AU


!Report Final Energy of System / Compare to Initial

!Loop through each body to find GPE of each
P = 0.
Vabs_Squared = 0.
!For each body
Do i = 0,n-1
        !For each other body
        Do j = 0,n-1
        !Ignore self
        If (i /= j) then
            !Find absolute separation from position vectors
            s(i,j) = ((r(i,0,0)-r(j,0,0))**2 + (r(i,1,0) -r(j,1,0))**2 + (r(i,2,0) - r(j,2,0))**2)**0.5
            !Sum calculated GPE
            P(i) = P(i) - G*M(i)*M(j)/s(i,j)
        End if
        End Do
        !Sum kinetic energy in each dimension (k)
        Do k = 0,2
            Vabs_Squared(i) = Vabs_Squared(i) + v(i,k,0)**2
        End Do
End Do
!Half to account for doubling up of counting in loop
P = P/2.
E1 = 0.
!Sum total energy
Do i = 0,n-1
    E1 = E1 +   0.5*M(i)*Vabs_Squared(i) + P(i)
End Do

!Report final energy and accuracy test
Write(6,*) ""
write(6,*) "Final Energy of System (J)"
write(6,*) E1
Write(6,*)
Write(6,*) "Percentage of Initial Energy Retained (Indicative of Sim Accuracy)"
write(6,*) E1/E0 * 100
Write(6,*) Step
Write(6,*) t
!Stops program from closing automatically
Write(6,*) ""
Write(6,*) "Clock Time (s) ~", (Time() - systime)
call cpu_time(finish)
Write(6,*) "CPU time (s)", finish - start
Write(6,*) ""
End Do
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





















end program NBody

