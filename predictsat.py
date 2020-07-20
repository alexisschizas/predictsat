# uses the skyfield library from : https://rhodesmill.org/skyfield/
from skyfield.api import EarthSatellite, Topos, load, Distance
import math
import numpy as np
from shapely.geometry import LineString, Polygon
import time
import os
import argparse, sys

start_time = time.time()    # used to compute the script total runtime

text = 'this script predicts the tracks of satellites on an observed region of the sky'

parser = argparse.ArgumentParser(description=text)
parser.add_argument("-i", "--input")
parser.add_argument("-o", "--output")
args = parser.parse_args()

if args.input:
   inputfile = args.input
if not args.input:
   print('predictsat.py -i <inputfile> -o <outputfolder>')
   sys.exit(2)
if args.output:
   outputfolder = args.output
if not args.output:
   print('predictsat.py -i <inputfile> -o <outputfolder>')
   sys.exit(2)

def initregionfile():   #initialize region file
    f = open(regionfilename, "w")
    f.write('# Region file format: DS9 version 4.1 \nfk5\n')
    f.close()

def writepointtoregion(ra, dec, objectname = '', color = 'Cyan'):   #creates a point in a ds9 region file
    f = open(regionfilename, 'a')
    f.write('point(' + str(ra) + ',' + str(dec) + ') # color='+color+' point = x text={'+objectname+'}\n')

def writelinetoregion(ra_start, dec_start,ra_stop,dec_stop, objectname = '', color = 'cyan'):   #creates a line in a ds9 region file
    f = open(regionfilename, 'a')
    f.write('line(' + str(ra_start) + ',' + str(dec_start) + ',' + str(ra_stop) + ',' +str(dec_stop) + ')# line=0 0 color='+color+' text={'+objectname+'}\n')

def initoutput():   #initialize output file
    f = open(str(outputfolder)+'/output.txt', 'w')
    f.write('')

def outputdata(sat, elevation):     #writes in the output file
    f = open(str(outputfolder)+'/output.txt', 'a')
    f.write(sat.name + '\t' + 'norad ID: ' + str(sat.model.satnum) +'\t' + 'mean altitude: ' + str(elevation) + '\n')

def processTLE(tlefilename):    #removes duplicates from a 3le file. The satellites need to be ordered by norad id number
    list_of_satellites = load.tle_file(tlefilename)
    print('Loaded', len(list_of_satellites), 'satellites')
    processed_list = []
    for i in range(len(list_of_satellites)):
        if i + 1 < len(list_of_satellites) and list_of_satellites[i + 1].model.satnum == list_of_satellites[i].model.satnum:
            continue
        else:
            processed_list.append(list_of_satellites[i])

    print('Removed duplicates,', len(processed_list), 'satellites remaining')
    return processed_list

def detectsat(sat, i, y, m, d, h, minstop, minstart):
    #creates a number i of segments that fit a satellite's trajectory.
    #returns the intersection points between the segments and the area of the sky that is observed,
    #the smallest distance between the segments and the area observed (0 if the intersect)
    #and an array of altitudes of the satellite
    minutes = np.arange(minstart, minstop, (minstop - minstart) / i).tolist()
    t_array = ts.utc(y, m, d, h, minutes)   #creates a time array of length i

    pts = []
    diff = sat - obs    #computes the position of the satellite in the sky
    topos = diff.at(t_array)    #computes the position of the satellite in the sky (array of positions since t is an array)
    ra, dec, dist = topos.radec()

    geo = sat.at(t_array)
    subpt = geo.subpoint()
    alt = subpt.elevation       #computes the satellite's altitude (array)

    ra, dec = ra._degrees, dec._degrees     #convert them to degrees (may create warning)

    for x in range(len(ra)):    #create a list of (ra, dec) points that can be used to create a Shapely LineString
        pts.append((ra[x],dec[x]))

    if np.isnan(ra[0]) or np.isnan(dec[0]) or np.isnan(ra[-1]) or np.isnan(dec[-1]):    #avoid nan values
        print('NAN value found')
        inter = area.intersection(area)
        dist = 10000000
        alt = 0
        return inter, dist, alt
    else:
        #create segments using the list of points, compute the intersection and the smallest distance between
        #the segments and the observed area
        satline = LineString(pts)
        inter = satline.intersection(area)
        dist = satline.distance(area)
        return inter, dist, alt

# ------------------------------------------- read the input file ------------------------------------------------------

f = open(inputfile, 'r')

path = str(outputfolder)
if not os.path.isdir(path):
    os.mkdir(path)

for line in f.readlines():
    #print(line)
    if 'observatory location' in line:
        strline = line.split('\t')
        lat, lon, obsalt = strline[1], strline[2], float(strline[3])

    if 'start time' in line:
        strline = line.split('\t')
        y,m,d,h,minstart = float(strline[1]),float(strline[2]),float(strline[3]),float(strline[4]),float(strline[5])

    if 'stop time' in line:
        strline = line.split('\t')
        y,m,d,h,minstop = float(strline[1]),float(strline[2]),float(strline[3]),float(strline[4]),float(strline[5])

    if 'Ra min/max, dec min/max' in line:
        strline = line.split('\t')
        ra_min, ra_max, dec_min, dec_max = float(strline[1]),float(strline[2]),float(strline[3]),float(strline[4])
# ---------------------------------------- end of reading the input file -------------------------------------------

ts = load.timescale()

# Sets the position of the observer. Is later used to compute the satellite's position in the sky
obs = Topos(lat,lon,elevation_m=obsalt)
print('observer is located at :', obs)

#process the TLE/3LE
tlefilename = '3le.txt'
list_of_satellites = processTLE(tlefilename)
regionfilename = str(outputfolder)+'/tracks.reg'  # ds9 region output
initregionfile()    # create region file, fill region file header
initoutput()        # create output file


#create some geometrical elements with Shapely
listofpts = [(ra_max, dec_min), (ra_max, dec_max),
                (ra_min, dec_max), (ra_min, dec_min)]
area = Polygon(listofpts)   #   creates a rectangle for the region observed

diag1 = LineString([listofpts[3],listofpts[1]])
diag2 = LineString([listofpts[0],listofpts[2]])
center = diag1.intersection(diag2)  #compute the center of the area, not used later in the code but could be useful

# main loop, iterate on all satellites
for s in range(len(list_of_satellites)):

    if 'DEB' in list_of_satellites[s].name:     #skip all the debris. You can remove this line if you want to include them
        continue

    inter, dist, alt = detectsat(list_of_satellites[s], 2, y, m, d, h, minstop, minstart)

    if dist < 30:
        # only consider satellites that are close to the region. This part can definitely be improved
        ra1, dec1,ra2,dec2 = np.nan, np.nan, np.nan, np.nan
        newinter, newdist, newalt = detectsat(list_of_satellites[s], 2, y, m, d, h, minstop, minstart)
        i = 2
        while i < 256 and newdist < 100:
            # keep fitting with more and more segments as long as you're not going too far (stops at 256)
            # this could be improved, maybe by observing the variation between dist and newdist, which tells you
            # if you're getting closer to the observed region or not
            i = 2*i
            newinter, newdist, newalt = detectsat(list_of_satellites[s], i, y, m, d, h, minstop, minstart)

            if len(newinter.coords) >= 2:   #if the segments intersect the observed region rectangle
                # get the 2 intersection points
                ra1 = newinter.coords[0][0]
                dec1 = newinter.coords[0][1]
                ra2 = newinter.coords[-1][0]
                dec2 = newinter.coords[-1][1]
                print('satellite :', list_of_satellites[s].name, '\nintersection points :',
                      '(',ra1,',',dec1,') ','(', ra2, ',', dec2, ')')
                #print(list_of_satellites[s].name, newdist)         #used to check dist of satellite to observation area
                print('satellite altitude :', (newalt.km[0]+newalt.km[-1])/2, 'km')

        if np.isnan(ra1) or np.isnan(dec1) or np.isnan(ra2) or np.isnan(dec2):  # avoid nan values
            continue
        else:   #fill the output and region file
            if 'DEB' in list_of_satellites[s].name:
                writelinetoregion(ra1, dec1, ra2, dec2, '', 'black')
                continue
            else:
                outputdata(list_of_satellites[s], (newalt.km[0]+newalt.km[-1])/2)
                writepointtoregion(ra1, dec1, '','green')  # start point
                writepointtoregion(ra2, dec2, '', 'red')  # end point -> deduce direction
                writelinetoregion(ra1, dec1, ra2, dec2, '', 'Cyan') # track
                writepointtoregion((ra1+ra2)/2, (dec1 + dec2)/2,list_of_satellites[s].name,'Cyan') # center point of the line for visual clarity on ds9

#remove temp files
if os.path.exists("deltat.data"):
    os.remove("deltat.data")

if os.path.exists("deltat.preds"):
    os.remove("deltat.preds")

if os.path.exists("Leap_Second.dat"):
    os.remove("Leap_Second.dat")

print('runtime :', (time.time() - start_time)/60, 'minutes')    #prints the total runtime