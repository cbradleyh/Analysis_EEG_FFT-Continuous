# Manuscript_MT-tACS-PAS_2017
Analysis of continuous EEG data to extract FFT information (before and after tACS)

Step by step: open these files and follow instructions in this order:

### A01_Stimulation_Offset.docx
Determine the exact offset of stimulation and the electrodes that show strong artefacts

### A02_Pre_Processing.m
Downsample, label with ‘active’ and ‘sham’ (if relevant), segment, filter (low-pass 100Hz. Notch 50Hz, demean).

### A03_Interpolation_ICA.m
Remove the ‘contaminated’ channels, interpolate the flat and noisy channels (automatic detection), run the ICA on filtered data.

### A04_ICA_Rejection.docx
Reject only eye movements (vertical and horizontal). To do this, you have to choose an electrode with which to correlate the signal, open each file, reject components and save.

### A05_ICA_Application.m
Applies the ICA correction.

### A06_Artefact_Rejection.docx
Reject manually all periods that have large amounts of speech, muscle artefact, etc...

### A07_Segment.m
Cuts into individual segments.

### A08_ReRef_FFT_Calculation.m
Decide if you want to change reference, specify FFT input, run FFT

### A09_SubSelect_Reformat.m
Define frequency band values, define subjects to exclude, reorganise data. This will not work if you don’t have the same organisation (4 segments with a baseline segment).

### A10_DataVis.m
Plot various things (spectra, group values, individual values, etc…). Not the most flexible/practical at the moment. If you can come up with a better script that works for you, do it!

### A11_FieldTrip_Permutation.m
Select the data you want to compare, select electrodes, setup an appropriate layout, set the statistical thresholds and values and perform a permutation analysis.

