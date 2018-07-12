# ae4d
Used for up to 4 dimensional processing of acoustoelectric data




%%%%% GETTING SENSITIVITY %%%%%
1) Load raw data
2) Enhance if desired
3) Choose the LF channel to compare to
4) Input into box 5 the pressure of the US transducer if known
5) Adjust gains for LF and HF channels
6) Press sensitivity and adjust cursor to the x,z value desired. Note that the sensitivity is read as a function of t at that point.
7) Readouts: 1 - Slope in microvolts per mA, 
             2 - Mean AE value $$ Not really useful
             3 - STD of AE values $$ Not really useful
             4 - Adjusted slope calculated by slope and pressure
             5 - Pressure (user input)
             6 - 2nd slope calculation, generally less accurate
             7 - The two signals (LF - Black, AE - Red) are shown in top right with given R^2 value as title
8) If slope is negative, check Invert on bottom left and press Modify to change the sign of the AE signal
9) If there is a time delay between LF and AE, use the TimeShift box in bottom left
Note: Is is best to use non-enveloped and non-rectified data for both LF and HF



%%%%% USING 2D OVERLAY %%%%%%%%
1) Load PE
2) Press Use
3) Filter PE data - Make sure there is only 1 time point in chopped version
4) Press To PEData
5) Load AE
6) Filter AE data - Make sure there is only 1 time point in chopped version
7) Press Overlay
8) Adjust height of AE signal compared to PE using Depth Shift



%%%%%%% USING IRADON %%%%%%
1) Create file 
2) Set Range to perform (more x and z points results in more accurate image)
3) Set Enhancesments
