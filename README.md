# NMDA-paper-2024
This repository contains all the code needed to reproduce all analyses and modeling in Drotos et al. 2024.

## Neuron model
The model neuron used here was based, with minimal changes, on the fast-spiking interneuron model described previously by Wang and Buzsáki (Wang and Buzsáki, 1996). The code found here can be used to run the model and generate the figures from the paper. 
* The main neuron model is contained in the function `NMDAmodel.m`. The inputs to this function are the modulation frequency of the stimulus along with the AMPA conductance, NMDA conductance, and leak conductance of the cell you would like to model. The `alpha`, `beta`, `h`, and `m` functions are required for the Hodgkin-Huxley channels. 
* The `plotMTFsmodel.m` will run the NMDA model at varying modulation frequencies to produce a temporal and rate MTF. 
* The `runModelSimulations.m` function will run a given model cell with set parameters 5 times and then average the MTFs.
* The `NMDAmodelTrains.m` function will produce the result from a train of inputs given AMPA and NMDA conductances with identical timing to that used in the slice electrophysiology experiments in the paper (5 stimuli at 30 Hz). 

## In vivo electrophysiology analyses
This in vivo data in this paper was analysed using custom MATLAB scripts. The functions here were used to analyze and plot the data for the paper. An explanation of each individual function is available by request. 
