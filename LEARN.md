Notes of learning process
-------------------------

- We have created several harmonics-<something>.hpp files in order to
  test the effect that each Jn has in the orbit.

- The "artificial" one helps us to intensify the effect of different
  zonal harmonics.

- Using artificial we have learnt that J2 has mainly the effect of
  precessing the orbit without significant changes in eccentricity and
  orbital size.

- The effect of J3 is changing eccentricity and semimajor axis almost
  continuously until collision.

- Combining J2 and J3 the orbit precesses and change its size but
  never collide with Earth.

- We have calculated the total energy of the satellite and the test
  program save the residula in "conservation.dat".  You can plot this
  value using gnuplot:

  > plot 'conservation.dat'

- A complete test can be performed using: 

  + Choosing the proper harmonics file in propagator.hpp
    (initPropagator routine).

  + Choosing the total integration time (duration in test.cpp).

  + Running the integration:

    make test.out

  + Plotting the evolution of the elements:
  
    python test.py

- Once an orbit is integrated it can be visualized in 3D using visual
  python and the script plot.py
