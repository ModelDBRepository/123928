This readme and the code were contributed by Timothee Masquelier
timothee.masquelier@alum.mit.edu
Sept 10th 2009

This code was used in:

Masquelier T, Hugues E, Deco G, Thorpe SJ (2009) Oscillations,
Phase-of-Firing Coding, and Spike Timing-Dependent Plasticity: An
Efficient Learning Scheme. J Neurosci 29(43):13484-13493
http://dx.doi.org/10.1523/JNEUROSCI.2207-09.2009

Feel free to use/modify but please cite us.

We use the Brian simulator described in: Goodman D, Brette R (2008)
Brian: a simulator for spiking neural networks in python. Front
Neuroinformatics 2:5.  and available from:
http://www.briansimulator.org/

This code has been tested with:
   - Brian 1.1.2
   - Python 2.5 and 2.6
   - Windows XP and Linux

Note that we did not use the STDP class provided in Brian. For
flexibility issues we preferred to code STDP manually (using embedded
C code for faster computations)

The base implementation corresponds to the 'all-to-all' additive STDP
of: Song S, Miller K, Abbott L (2000) Competitive hebbian learning
through spike-timing-dependent synaptic plasticity. Nat Neurosci

A 'nearest-spike' implementation is also provided.

Use mu > 0 to interpolate with multiplicative STDP as in: Gutig R,
Aharonov R, Rotter S, Sompolinsky H (2003) Learning input correlations
through nonlinear temporally asymmetric hebbian plasticity. J Neurosci

Rate-based homeostatic terms w_in and w_out can also be added as in:
Kempter R, Gerstner W, van Hemmen JL (1999) Hebbian learning and
spiking neurons. Phys Rev E

Main file (should be called like that "python -i main.py")
Calls init.py to set the parameters. The current configuration
corresponds to the first 3 seconds of the oscillation case (see Fig. 5
in the paper).

Instantiate the all the needed neurons (input layer, mirror layer
(copy of the input layer useful for implementation issues), output
layer), connect them, and finally runs the Brian simulator.

All the data files are read and dumped in ../data

There are two main modes:

   recomputeSpikeList = True (used in the paper)
       input layer is computed from activation levels in files
       inputValues.###.txt (the number is the random seed) format:
       '%010.5f ' (time) '%07.5f ' (value for neuron 0) '%07.5f '
       (value for neuron 1) ...  The files used in the paper are
       available upon requests (timothee.masquelier@alum.mit.edu)

   recomputeSpikeList = False
       spikes from the input layer are read from spikeList.###.###.mat
       files (first number is the random seed, second number is the
       period number).  each file contain a (n,2) array called sl :
       neuron (start from 0), and spike time in s This mode thus
       allows you to test STDP with input spike trains generated
       somewhere else (for eg with Matlab)

To take advantage of Python vectorization we can simulate multiple
output neuron in parallel, with different gmax and LTD/LTP ratios
(this is useful for parameter exploration)

nG is the number of gmax values
nR is number of ratio LTD/LTP
M = nG*nR is the number of output neurons, numbered like that [
(r_0,g_0)...(r_0,g_nG),(r_1,g_0)...(r_1,g_nG),...,(r_nR,g_0)...(r_nR,g_nG)]

Note that instead of LTD/LTP ratios you can also explore various w_in/w_out ratios.

If you want to use resets set useReset = True in init.py and provide a
reset.###.mat file with a (1,n) array containing reset times in s and
called resetTimes

For plotting it may be useful to provide a realValuedPattern.XXX.mat
file with the repeating pattern activation levels in a (1,n) array
called realValuedPattern
For plotting and mutual info computation it may be useful to provide a
patternPeriod.XXX.mat file with a (n,2) array (start,end) in s called
patternPeriod

other Python files:
init.py contains all the parameters
analyze.py graphical plotting
mutualInfo.py computes, plots and dumps the mutual info between
       presence of the stimulus and postsynaptic responses
saveWeight.py dump final weigth in weight.###.mat
customrefractoriness.py Brian file to handle both a refractory period
       and a user-defined reset function
toMatlab.py (under development): exports data in a mat file for later
       use with Matlab
pickleAll.py (underdevelopment): dump all variables
unpickleAll.py (underdevelopment): load all (previously dumped)
       variables
restore.py (underdevelopment): restore the final state of a previous
       simulation

Final note: I was new to Python when I started this project. I'm sure
many things could be coded more efficiently.
