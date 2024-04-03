# PFF-Analysis: Peak Finding and Fitting for System Identification

These MATLAB routines are used to identify nonlinear frequency and damping trends from transient time series data.

This code is provided free of charge to aid in research. No guarantee is made about the results.

If you use this code, please cite the journal paper:
```
@article{jinIdentificationInstantaneousFrequency2020,
	title = {Identification of Instantaneous Frequency and Damping From Transient Decay Data},
	author = {Jin, Mengshi and Chen, Wei and Brake, Matthew R. W. and Song, Hanwen},
	date = {2020-06-26},
	journaltitle = {Journal of Vibration and Acoustics},
	shortjournal = {Journal of Vibration and Acoustics},
	volume = {142},
	number = {051111},
	issn = {1048-9002},
	doi = {10.1115/1.4047416}
}
```

## Usage 

### Graphical User Interface

1. Clone the repository to your local directory. 
2. Run the script PFF_Analyze.m in MATLAB. 
3. Follow the instructions in the GUI to processes the data. Example LMS data is included in the `data` folder here.
4. After running you will be interested in the variables `Amp`, `Freq`, `Damp`, `Time`, which are amplitude (same units as input), frequency [Hz], damping ratio [fraction of critical], and time [s].

### Function 

There is now a clear function that can be used in `pff_function.m`, but this is not the cleanest. 
Run the file `tests/test_pff_function.m` to verify that you get the expected results with this function.

## Future Work

1. Write additional verification tests
2. Clean-up code
