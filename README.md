# PFF-Analysis: Peak Finding and Fitting for System Identification

These MATLAB routines are used to identify nonlinear frequency and damping trends from transient time series data.

To Do: 
1. Insert Paper referece
3. Write tests to verify that code behaves as expected.
2. Refactor GUI code so that there is a clear API that one can use if they do want to use the GUI. 

There is now a clear function that can be used in `pff_function.m`, but this is not the cleanest. 
Run the file `tests/test_pff_function.m` to verify that you get the expected results with this function.

## Usage 

1. Clone the repository to your local directory. 
2. Run the script PFF_Analyze.m in MATLAB. 
3. Follow the instructions in the GUI to processes the data. Example LMS data is included in the `data` folder here.
4. After running you will be interested in the variables `Amp`, `Freq`, `Damp`, `Time`, which are amplitude (same units as input), frequency [Hz], damping ratio [fraction of critical], and time [s].
