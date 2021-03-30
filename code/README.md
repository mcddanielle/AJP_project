
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

To see the authors solutions to exercises 2 and 5, set line 43 in solutions_exercises.py

```solve_exercise = "2a"```

This flag will choose simulation settings.  Then simply run 

```python3 solve_exercise.py```

### About figures in manuscript (and codes)
----------
```plot_landscape.py```

#### Figure 1

(fig1.pdf) which is a schematic birds eye view of the system.

two panels, side by side
(a) birds eye view of Sy by Sx for a substrate period of Sy/3
(b) U(y) for the same period, calculation of force Fp = -grad(U)

the magnitude of U0 / Ap is scaled for viewing rather than simulation.

-------------

```MD_colloid.py```

performs the simulations discussed in the paper.  A single colloid on a corrugated or washboard potential/substrate is subject to a periodic applied driving force.  System properties are plotted.

-------------

```solutions_exercises.py```

performs the simulations discussed in the problems section using calls to the ```MD_colloid``` module

-------------
#### Figure 2

(fig2.pdf) basic dynamics
(a) driving force versus time
(b) y/lambda versus time, where y is the particle position and lambda is the substrate period

parameters
AP  = 0.1
FAC = 0.05
FDC = 0.1
f = 0.01

dt = 0.1
maxtime = 4000 (integer timesteps)

It this frequency and driving regime the steps are fractional.

---------------
#### Figure 3 

(fig3.pdf) sweep FDC for fixed FAC

parameters
AP  = 0.1
(a) FAC = 0.0
(b) FAC = 0.05
FDC = 0.0 - 0.3
f = 0.01

---------------
#### Figures 4 

(fig4.pdf) phase variables

parameters

AP  = 0.1
FAC = 0.05
f = 0.01
dt = 0.1
maxtime = 400

(a) FDC = 0.04
(b) FDC = 0.07
(c) FDC = 0.1
(d) FDC = 0.125

----------------

### About solutions to problems

#### Problem 1 - Write your own MD code

refers to the complete MD code contained within ```MD_colloid.py```

--------------
#### Problem 2 - Exploring model parameters

refers to Problem 2a-c (fig5.pdf, fig6.pdf, fig7.pdf) which may be generated with line 43 in solutions_problems,py

---------------
##### Problem 2a - variable FAC

parameters['maxtime']=3000        
F_AC = [0.2, 0.3, 0.4]

---------------
##### Problem 2b - variable frequency

parameters['maxtime']=4500        #total time steps in simulation
frequency = [0.1, 0.015, 0.01, 0.005, 0.002, 0.001]
------------------
##### Problem 2c - variable FAC - sweep FDC - produce Shapiro steps

parameters['freq'] = 0.1
delta_Fdc = 0.01
F_AC = [0.1,0.2,0.3,0.4]
Fdc_max=0.6+delta_Fdc

---------------
#### Problem 3 - Drag models and Reynolds numbers

Simple analytical calculation.  See ```calc_Reynolds.py``` for solution.
No figure produced.

---------------
#### Problem 4 - Equation of motion

Simple analytical calculation.  No figure produced.

------------------
#### Problem 5 - Brownian motion

AP  = 0.1
FAC = 0.0
FDC = 0.0
T/Ap = 4.0, 5.0, 6.0

maxtime = 300000
writemovietime = 1000 

-------------------


