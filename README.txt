To run this script, you need the following python packages : skyfield, shapely, numpy and time. You can simply install them with pip. 

You can call the script using python predictsat.py -i <inputfile> -o <outputfolder>. 
The input file contains the data needed for the script to run (in tsv format). To use the script on your own data, 
just change the values of the observatory location, start and stop time as well as min/max Ra/dec of the region of the sky observed.
Be careful not to modify the first entry of each line (observatory location, start time (YMDHMin), stop time (YMDHMin), Ra min/max, dec min/max).

You will also need to provide a .txt file containing the TLEs of the satellites you want to check. The file needs to be called '3le.txt' and be 
placed in the same folder as the script itself. 
The format required is the 3le format containing the name of the satellite in the first line, and then the two line elements of the satellite 
(see example just below).

SWISSCUBE               
1 35932U 09051B   20116.45284354  .00000113  00000-0  36666-4 0  9998
2 35932  98.6092 297.0369 0006399 293.2581  66.7942 14.56380908562114

Please note that some websites call this format TLE, while others call it 3LE (since it contains a third line with just the satellite name).
You can use a website like space-track.org, and use the query builder to get TLEs for the specific time of your observation. The script 
should remove any duplicates in the TLEs, as long as they are ordered by NORAD ID (you can specify how to order them with the space-track.org
query builder).

The script creates two output files, one called 'tracks.reg' which is a region file you can open with SAOImageds9 and that should plot the
tracks of the satellites on top of your FITS file. The other file, called 'output.txt', stores data like the satellite name, norad ID and 
mean altitude during the time it crosses the observed region. 

Note that you might get a warning in the terminal telling you that an invalid value has been encountered