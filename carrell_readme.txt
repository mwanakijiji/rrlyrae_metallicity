Notes on the bkgrnd.cc program:

1. to compile use the command:
     g++ -o bkgrnd bkgrnd.cc

2. to see the list of options run the program with no input:
     ./bkgrnd

3. the input file needs to have the input and output directories 
   in the first line followed by the filenames one each per line,
   for example:
    indata outdata
    file1.dat
    file2.dat
       .
       .
    fileN.dat

 !!IMPORTANT!! the *same* filename is used, so if you give the 
  same input and output directories it will WRITE OVER the 
  original data, make sure the input and output directories 
  are different to keep the original data.

4. the runtime options:
    --boxcar - smooth the data before estimating background
               (yes if the option is given, no if not)
    --smooth # - integer used for smoothing. 33 is default
    --sismoo # - boxcar smoothing during background estimation
                 (this is given as +/- pixels, so a value of 
                  1 will smooth to 3 pixels)
    --no-plot - suppresses plotting of the output
                (default is to plot results)
    --delay # - how many seconds to show results before moving on
                (a value of -1 will require you to hit Enter in 
                 the command window to continue, default is 3)

5. Examples.

to run with:

the default options:
./bkgrnd in.list

a smoothing parameter of 25:
./bkgrnd --smooth 25 in.list

a plotting delay of 5 seconds:
./bkgrnd --delay 5 in.list

a smoothing parameter of 25, no simultaneous smoothing and no plotting:
./bkgrnd --smooth 25 --sismoo 0 --no-plot in.list
