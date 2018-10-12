# ae4d
Used for up to 4 dimensional processing of acoustoelectric data

%%%%% Basic Usage %%%%%
***Creating filtered 4D AE Matrix***
Top left ==> Choosing filters
In order to use a fast or slow time filter, check their corresponding On boxes.
  To Match filter check the match box
  To Bandpass filter, choose range in the boxes and leave match unchecked.
  The transducer option dropdown menu only changes the default values for bandpass filter right now
Check the OneMHz box if using single element US transducer
Check TC box and then press TC button to skew frequency spectrum to lower frequencies
  In TC GUI choose the thickness, speed of sound in bone and US frequency range then press Set
In order to filter multiple channels of AE data, enter the channel numbers in HF chans box
  Leave this box empty if only one channel. Leaving it empty with multiple channels will take the 1st channel
  Channels will be remapped such that if you recorded from 0, 2 and 7, they will be considered 1, 2 and 3
If you are saving very large file (>1.5GB), check the Large box, otherwise the 4D AE variable will not be saveable 
Press Create when ready 

***Loading filtered AE data***
Press load and choose 4d_data set
If using multiple channels, enter then number of channels you want to load into the box to the left of MERGE and press Merge
  Toggle Max Range box to set the max ranges for all dimensions

***Chopping and merging***
After loading an AE file, the 4D data will be in a matrix called Xfilt
In order to process this matrix without changing any of its contents, choose a range of X, Y, Z and T values and press Chop
  This will create a new variable called X_c
  Insofar as Use Chopped is checked, any filtering peformed will be done on X_c and not Xfilt
  To overwrite X_c from Xfilt just press Chop again
    This is an effective way to change ranges and enhancements without having to reload Xfilt
A similar process can be used for dealing with multiple channels of AE data
  After loading multiple channels using Merge, a cell array called Xmerged is created that contains the Xfilt for each channel
  Enter the channel you want to further process in the box to the right of Merge and press Use to extract the Xfilt related to that channel
  That particular Xfilt can then be chopped and modified like normal
  Press return to overwrite the copy of xfilt for that channel in the Xmerged cell array.
  
  ***Looking at Pulse Echos***
Pulse echos can be loaded from either verasonics beamformed or raw RF data
  To look at verasonics beamformed data, check the BSQ box and press Create PE
  After processing, press Load PE and select the PE file
  After loading press Use to convert the PE data to an Xfilt variable to view and modify as if it were AE data
  To return the PE data to its original PE matrix press ToPEData
  Overlay was designed to overlay PE with AE data, but this option got replaced by Stitch (see below)
  
  ***Enhancing data***
To intepolate or mean/median filter data check their associated boxes at the top.
  Median and Mean filters need to be odd numbers
  Make sure Use Chopped is checked to modify X_c rather than Xfilt
    This is very important if interpolating since Xfilt tends to already be a very large variable
  Check Square to make the mm/pixel ratio for all dimensions the same during interpolation
    Currently the script takes the average mm/pixel ratio between all axes
    Note that the number in each interpolate box represents how much to interpolate the variable by
      In order to set the mm/pixel ratio to a particular value, do the math.
  
***Inverting polarity***
Sometimes you think a circle should be red after basebanding but its actually blue, to fix this:
  Check Invert in the bottom left panel and press Modify
  Make sure Time Shift and BB Freq are set to 0 to prevent them from also acting
  
***Basebanding***
Speaking of which, in order to baseband your data 
  Enter the demodulation frequency in the BB Freq box
  Check the S_env and dB boxes
  Press Modify
  To view:
    Press dB, check hc (hotcold), and set the Rng to something like -10 10 before plotting
  Taking envelope before basebanding will screw things up
  Using Mean filter before basebanding will screw things up
  Interpolating before basebanding will screw things up
  
  
***Advanced basebanding***
  Just talk to me in person
  
***Time or Shifting***
Sometimes you plot sensitivity and the LF and AE waveforms are a sample or two misaligned
  To fix this:
    Enter the number of samples in the Time Shift box you want to shift the AE data by
    Press plot to view new image
    Note that this will shift it from where it currently is and not from its original position.

  
***Regular plotting***
To simply plot a 2D image, select the dimensions you want in the Image box and press Plot
  Enter the range (min max) for each dimension in the box on the left underneath Movie axes menu
  The box to the right of the range box is the value to use for that dimension when plotting other dimensions
  The Rng box sets the C range for all plots.
    Leaving Rng blank will result in basic linear scale of all data
    Press dB to change everything to dB scale
      Set range to -10 0 or something like that
To plot all dimensions either press Plot4 or click on the image after using Plot
  Clicking on an image will update the used coordinate for all axes
Check the Show FFT box to see the spectral data of the AE image
  Use the ZT plot to get slowtime and fasttime FFT


***Movies***
Much like using Plot, select the dimenions you want to use in the Movie menu and press Movie
  Check the All box to view all dimensions simulatneously
  Check the Show LF box to view the LF signal in relation to the movie frames.
  
  
***Cool things you can do with plots
    


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
 1) Press Stitch button.
 2) The rest is pretty self-explanitory

 %%%%%%% USING IRADON %%%%%%
1) Create file 
2) Set Range to perform (more x and z points results in more accurate image)
3) Set Enhancements
