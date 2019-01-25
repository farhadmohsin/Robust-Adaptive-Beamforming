# Robust-Adaptive-Beamforming

A robust and optimized system architecture has been developed and designed for adaptive beamformer with a Microphone Array. The system includes following subsystems -
MMSE STSA Estimator, DOI (Direction of Interest) Estimator and an Adaptive Beamformer. 
The system architecture had been implemented and tested for Xtensa Processor which was configured for HiFi-2 DSP Standard for audio
processing.

This work was done as participants of the Cadence Design Contest. A variant of this work was submitted and accepted in the 2015 18th International Conference on Computer and Information Technology (ICCIT). A raw report written before the submission in the conference has been included in the repo.
https://ieeexplore.ieee.org/abstract/document/7488064

The repository consists of three items.

1. beamformer_singlewindow.c: This code was done specifically for the HiFi-2 DSP Processor. We provide with inputs for this code as well.
2. The input for the code are text files, defining the audio input. 
- window_in: eight arrays of double precision point variables of size 512. These are the current window's input from the eight channels.
- window_prev: This carries eight arrays of 307 double precision point variables. This arrays contain values provided from the run of the previous window. For the first window's run, these arrays need to be set at 0.
This file next carries an array, Y of (512-307) = 205 elements, that are also carried from the output of the last run. This, again is initialized as 0 for the first window.
- noise_in:  This is the noise spectrum estimated from the silent zone for the use of MMSE. This contains two arrays of 257 elements, Noise and Lambda.
These noise could be generated easily from code defined in Beamformer.cpp, lines 70 onwards.

3. Beamformer.cpp: This contains a .cpp code, executable and easily convertible to Xtensa processor's need. The code is full of comments, to ensure understanding.
This code runs for a number of windows at once and gives a complete audio output.
The I/O for this file is defined in the code comment section.
