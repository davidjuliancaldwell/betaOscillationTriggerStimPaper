### Code to analyze evoked potentials from beta oscillation triggered direct electrical stimulation in humans.

This repository contains MATLAB and R code to analyze the data from the beta oscillation triggered stimulation experiments.

The main script is ***master_script_betaStim.m***, which calls other sub scripts and analysis functions.

---

### Finding stimuli timing

The ***find_stims*** folder has code to figure out when the stimuli were delivered and build tables to help further data extraction. 

---

### Phase R_analysis_scripts

The function ***B_phaseCalc_allChans_processed.m*** is the key function for the nonlinear sinusoid fits to estimate the phase of delivery across all the channels.

---

### Peak Extraction

The ***peak_extraction*** folder has scripts illustrating sweeping through the subjects and extracting peak to peak voltages.


---

### R analysis

The ***R_analysis_scripts*** folder contains the R scripts required to fit linear mixed models and generate statistical plots after the data structure generated from ***master_script_betaStim*** has been run

---

### Filter analysis

Code to analyze the performance of the real time filter on the TDT is in the ***test_real_time_filter*** folder
