This python codes implement the homeostatic model in 

______________________________________________________________________________________________________________________
"Gating of homeostatic regulation of intrinsic excitability produces cryptic long-term storage of prior perturbations"
Leandro M, Alonso*, Mara C.P. Rue*, and Eve Marder
______________________________________________________________________________________________________________________

The model definition is in 

- singlecell_liu_bound_dynTauG.py

To demonstrate how the model works, it is convenient to start a simulation from random initial conductances and wait until it self-assembles into a periodic bursting pattern. 

The script "integrateContinue_bound_dynTauG_from_random.py" starts the model from random and integrates it for 10 minutes. The solution will be stored at ./simulation/init-<randomname>/
The script "plotSimulation_masters_with_alpha_and_sf_2022_V2.py" loads the simulation data and produces a plot. 

These two scripts are sufficient to reproduce the results in Fig. 4 as follows

> python integrateContinue_bound_dynTauG_from_random.py
> python plotSimulation_masters_with_alpha_and_sf_2022_V2.py

________________________________________________________________________________________________________________________
**********************************
* Description of the other files *
**********************************

- singlecell_liu_bound_dynTauG.py
Model definition 

- integrateContinue_bound_dynTauG_from_random.py
Script to integrate the model from random initial conditions 

- plotSimulation_masters_with_alpha_and_sf_2022_V2.py
Script to load simulations and produce a plot

- analizeSolution.py
Auxiliary functions to detect spikes, etc. 

- auxfunctions.py
Auxiliary functions

- currents_visualization.py
Scripts to plot currentscapes

- saveCurrents_bound_dynTauG_V2.py
Auxiliary scripts to save the currents needed to plot currentscapes

- plot.currentscapes.2022.py
This script shows how to plot currentscapes

