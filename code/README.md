
## Molecular dynamics simulation of a driven synchronized particle
### American Journal of Physics

Danielle McDermott
Tiare Guerrero

The supplementary codes were written in Python3.  To run these codes, install a Python 3 interpreter.

https://wiki.python.org/moin/BeginnersGuide/Download
https://www.anaconda.com/

The programs use common Python packages including numpy and matplotlib.  These come preinstalled in the Anaconda interpreter, but may need to be installed with another interpreter.

### Getting started

To create Figure 1 (fig1.pdf), simply run the following command

```python3 plot_landscape.py```

To create Figures 2-4, set line 578 in MD_colloid.py

```make_fig = 4```

This flag will choose simulation settings.  Then simply run 

```python3 MD_colloid.py```

#-------------------------------------------------------------------
plot_landscape.py
creates Figure 1 (fig1.pdf), which is a schematic birds eye view of the system.

two panels, side by side
(a) birds eye view of Sy by Sx for a substrate period of Sy/3
(b) U(y) for the same period, calculation of force Fp = -grad(U)

to fix - could downgrade the "D" to "d" in F^D(t) (and all in text)

the magnitude of U0 / Ap is scaled for viewing rather than simulation.

#-------------------------------------------------------------------
MD_colloid.py performs the simulations discussed in Section IIIA of the paper.  the A single colloid on a corrugated or washboard potential/substrate is subject to a periodic applied driving force.  An animation of the motion and the system properties are plotted.

#-----------------------------------------------------------------
Figure 2 - basic dynamics
(a) driving force versus time
(b) y/lambda versus time

parameters
AP  = 0.1
FAC = 0.05
FDC = 0.1
f = 0.01

dt = 0.1
maxtime = 4000

It this frequency and driving regime the steps are fractional.

#-----------------------------------------------------------------
Figure 3 - sweep FDC for fixed FAC

parameters
AP  = 0.1
(a) FAC = 0.0
(b) FAC = 0.05
FDC = 0.0 - 0.3
f = 0.01

#-----------------------------------------------------------------
Figures 4 - phase variables

parameters
(a)
AP  = 0.1
FAC = 0.05
FDC = 0.1
f = 0.01

(b) ...

dt = 0.1
maxtime = 400

#-----------------------------------------------------------------

Solutions to Exercises

Exercise 1

refers to the complete MD code contained within MD_colloid.py

#-----------------------------------------------------------------
Fig 5 - Brownian

AP  = 0.1
FAC = 0.0
FDC = 0.0
T/Ap = 3, 3.5, 4

maxtime = 
writemovietime =

#-----------------------------------------------------------------
Figures 6 - variable frequency

parameters
AP  = 0.1
(a) FAC = 0.0
(b) FAC = 0.05
FDC = 0.0 - 0.3
f = 0.1, 0.05, 0.015, 0.005, 0.001

dt = 0.1
maxtime = 400

#-----------------------------------------------------------------
Figure 7 - variable FAC

#-----------------------------------------------------------------
Figure 8 - variable FAC vs. FDC - as in juniper?

#-----------------------------------------------------------------
Figure 9 - step width?

