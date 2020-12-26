
## Molecular dynamics simulation of synchronization in driven particles
### American Journal of Physics

Danielle McDermott
Tiare Guerrero

The supplementary codes were written in Python3.  To run these codes, install a Python 3 interpreter.

https://wiki.python.org/moin/BeginnersGuide/Download
https://www.anaconda.com/

The programs use common Python packages numpy and matplotlib.  These come preinstalled in the Anaconda interpreter, but may need to be installed with another interpreter.

MD_colloid.py performs the simulations discussed in Section IIIA of the paper.  A single colloid on a corrugated or washboard potential/substrate is subject to a periodic applied driving force.  An animation of the motion and the system properties are plotted.

MD_parameter.dat contains the presets to run MD_colloid.py as discussed in Section IIIA.  The values stored in this file are read into a Python dict to be passed to the various subroutines of the MD algorithm.