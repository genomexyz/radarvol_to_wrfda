'''
 -------------------------------------------------------------------------------------
   vol2wrfdainput.py

   Purpose: This script is generating ob.radar as wrfdainput
			from specified radar using wradlib module

   Status: Development 

   History:
	  --/2016 - Donaldi Permana - Initial development in radar data extraction
	  05/2018 - Muhammad Ryan - Writing radar.stat, radar.temp, and ob radar
	  12/2018 - Abdullah Ali - Add regriding in radar data extraction
							   Joint between radar data extraction and wrfdainput

   Copyright 2018. BMKG. All rights reserved. 
 ------------------------------------------------------------------------------------
 '''

from netCDF4 import Dataset
from time import gmtime, strftime
from datetime import datetime
import wradlib as wrl
import numpy as np
import csv,warnings,os,glob
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=RuntimeWarning)

def rungema(radartime,radarpath):
	listpcdfpath=list()
	runtime=radartime.strftime("%Y%m%d%H%M")
	print(radarpath)
	print(os.listdir(radarpath))
	for i in os.listdir(radarpath):
		if os.path.isfile(os.path.join(radarpath,i)) and runtime in i:
			fpath=radarpath+'/'+i
			listpcdfpath.append(fpath)
	return listpcdfpath

def GEMA(listpcdfpath,radarpath,res,rmax,sweep):
	f0=wrl.util.get_wradlib_data_file(listpcdfpath[0])
	raw=wrl.io.read_Rainbow(f0)
	try :	
		llon = float(raw['volume']['sensorinfo']['lon'])
		llat = float(raw['volume']['sensorinfo']['lat'])
		lalt = float(raw['volume']['sensorinfo']['alt'])
	except :
		llon = float(raw['volume']['radarinfo']['@lon'])
		llat = float(raw['volume']['radarinfo']['@lat'])
		lalt = float(raw['volume']['radarinfo']['@alt'])
	res_coords=res/111229.
	xmax,xmin=llon+(rmax/111229),llon-(rmax/111229)
	ymax,ymin=llat+(rmax/111229),llat-(rmax/111229)
	n_grid=np.floor(((xmax-xmin)/res_coords)+1)
	x_grid=np.linspace(xmax,xmin,n_grid)
	y_grid=np.linspace(ymax,ymin,n_grid)			

	for fpath in listpcdfpath:
		filename=fpath[len(radarpath)+1:]
		f=wrl.util.get_wradlib_data_file(fpath)
		raw=wrl.io.read_Rainbow(f)
		
		# load site radar attribute
		try :	
			llon = float(raw['volume']['sensorinfo']['lon'])
			llat = float(raw['volume']['sensorinfo']['lat'])
			lalt = float(raw['volume']['sensorinfo']['alt'])
		except :
			llon = float(raw['volume']['radarinfo']['@lon'])
			llat = float(raw['volume']['radarinfo']['@lat'])
			lalt = float(raw['volume']['radarinfo']['@alt'])
		sitecoords=(llon,llat,lalt)

		i=sweep
		try :	
			elevation = float(raw['volume']['scan']['slice'][i]['posangle'])
		except :
			elevation = float(raw['volume']['scan']['slice'][0]['posangle'])
		strsweep=str(i+1) 

		# load azimuth data
		try:
			azi = raw['volume']['scan']['slice'][i]['slicedata']['rayinfo']['data']
			azidepth = float(raw['volume']['scan']['slice'][i]['slicedata']['rayinfo']['@depth'])
			azirange = float(raw['volume']['scan']['slice'][i]['slicedata']['rayinfo']['@rays'])  
		except:
			azi0 = raw['volume']['scan']['slice'][i]['slicedata']['rayinfo'][0]['data']
			azi1 = raw['volume']['scan']['slice'][i]['slicedata']['rayinfo'][1]['data']
			azi = (azi0/2) + (azi1/2)
			del azi0, azi1
			azidepth = float(raw['volume']['scan']['slice'][i]['slicedata']['rayinfo'][0]['@depth'])
			azirange = float(raw['volume']['scan']['slice'][i]['slicedata']['rayinfo'][0]['@rays'])			
		try:
			azires = float(raw['volume']['scan']['slice'][i]['anglestep'])
		except:
			azires = float(raw['volume']['scan']['slice'][0]['anglestep'])
		azi = (azi * azirange / 2**azidepth) * azires

		# load range data
		try:
			stoprange = float(raw['volume']['scan']['slice'][i]['stoprange'])
			rangestep = float(raw['volume']['scan']['slice'][i]['rangestep'])
		except:
			stoprange = float(raw['volume']['scan']['slice'][0]['stoprange'])
			rangestep = float(raw['volume']['scan']['slice'][0]['rangestep'])
		r = np.arange(0, stoprange, rangestep)*1000
			
		# projection for lat and lon value
		polargrid=np.meshgrid(r,azi)
		x,y,z=wrl.georef.polar2lonlatalt_n(polargrid[0],polargrid[1],elevation,sitecoords)

		# regriding
		grid_xy=np.meshgrid(x_grid,y_grid)
		xgrid=grid_xy[0]
		ygrid=grid_xy[1]
		grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()
		xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
		radius=r[np.size(r)-1]
		center=[x.mean(),y.mean()]
		
		# load radar data (depend on what file you load)
		if fpath[-7:]=='dBZ.vol':
			print 'Extracting data '+filename+' : SWEEP-'+strsweep +' at Elevation Angle '+str(elevation)+'  deg ...'
			data = raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['data']
			datadepth = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@depth'])
			datamin = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@min'])
			datamax = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@max'])
			data = datamin + data * (datamax - datamin) / 2 ** datadepth
			# data masking and preprocessing
			clutter=wrl.clutter.filter_gabella(data, tr1=6, n_p=6, tr2=1.3, rm_nans=False)
			try:
				data_noclutter=wrl.ipol.interpolate_polar(data, clutter, Interpolator = wrl.ipol.Linear)
			except:
				data_noclutter=data
			data_dbz=data_noclutter
			# if the shape of data and georeferencing is not the same	
			a,b=np.shape(x)
			c,d=np.shape(data_dbz)
			selisih1=a-c
			selisih2=b-d
			if selisih2>0:
				print ('Matching data and coordinate shape...')
				data_=np.zeros((a,b))
				for k in range(c):
					for j in range(d):
						data_[k,j]=data_dbz[k,j]
					for ii in range(selisih2):
						data_[c-1,d+ii]=np.nan
					data_dbz=data_
			if selisih1>0:
				print ('Matching data and coordinate shape...')
				for ii in range(selisih1):
					data_[c+ii,d-1]=np.nan
				data_dbz=data_
			# interpolate polar coordinate data to cartesian cordinate
			# option : Idw, Linear, Nearest
			gridded_dbz = wrl.comp.togrid(xy, grid_xy, radius, center, data_dbz.ravel(), wrl.ipol.Linear)	
			gridded_datadbz = np.ma.masked_invalid(gridded_dbz).reshape((len(x_grid), len(y_grid)))
							
		elif fpath[-5:]=='V.vol':
			print 'Extracting data '+filename+'   : SWEEP-'+strsweep +' at Elevation Angle '+str(elevation)+'  deg ...'
			data = raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['data']
			datadepth = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@depth'])
			datamin = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@min'])
			datamax = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@max'])
			data = datamin + data * (datamax - datamin) / 2 ** datadepth
			data_v = data
			# interpolate polar coordinate data to cartesian cordinate
			# option : Idw, Linear, Nearest
			gridded_v = wrl.comp.togrid(xy, grid_xy, radius, center, data_v.ravel(), wrl.ipol.Linear)	
			gridded_datav = np.ma.masked_invalid(gridded_v).reshape((len(x_grid), len(y_grid)))
				
	return gridded_datadbz,gridded_datav,sitecoords,datetime,elevation,xgrid,ygrid,xmin,xmax,ymin,ymax,z

def write_radartemp(wrffile,radarex,radarstatdata,rest,scantotal,radartime,radarpath):
	# extract wrf domain
	wrfread = Dataset(wrffile, 'r')
	lonwrf = wrfread.variables['XLONG'][0]
	latwrf = wrfread.variables['XLAT'][0]
	minlonwrf, maxlonwrf = lonwrf[0,0], lonwrf[0,-1]
	minlatwrf, maxlatwrf = latwrf[0,0], latwrf[-1,0]

	# radar.stat
	print 'Writting radar.stat'
	listpcdfpath=rungema(radartime,radarpath)
	radaropen = GEMA(listpcdfpath,radarpath,res_radar,rmax,0)
	radstat = open(radarstatdata, 'w')
	radstat.write(str(round(radaropen[2][0],3))+','+str(round(radaropen[2][1],3))+','+str(round(radaropen[2][2],3)))

	# radar.temp
	tulisfile = open(radarex, 'w')
	datawriter = csv.writer(tulisfile, delimiter=',')

	'''
	radaroutput[0]  : gridded_datadbz
	radaroutput[1]  : gridded_datav
	radaroutput[2]  : sitecoords
	radaroutput[3]  : datetime
	radaroutput[4]  : elevation
	radaroutput[5]  : xgrid
	radaroutput[6]  : ygrid
	radaroutput[7]  : xmin
	radaroutput[8]  : xmax
	radaroutput[9]  : ymin
	radaroutput[10] : ymax
	radaroutput[11] : z
	radaroutput[12] : scantotal

	'''
	print '\nWritting radar.temp'
	for i in xrange(scantotal):
		radaroutput=GEMA(listpcdfpath,radarpath,res_radar,rmax,i)
		for j in xrange(len(radaroutput[0])):
			for k in xrange(len(radaroutput[0][1])):
				if radaroutput[1][j,k] and radaroutput[0][j,k] and radaroutput[5][j,k] > minlonwrf and radaroutput[5][j,k] < maxlonwrf and \
				   radaroutput[6][j,k] > minlatwrf and radaroutput[6][j,k] < maxlatwrf:
					datawriter.writerow([round(radaroutput[5][j,k],3), round(radaroutput[6][j,k],3), round(radaroutput[0][j,k],1), round(radaroutput[0][j,k],3), round(radaroutput[1][j,k],3)])
				elif not(radaroutput[0][j,k]) and radaroutput[1][j,k] and radaroutput[5][j,k] > minlonwrf and radaroutput[5][j,k] < maxlonwrf and \
					 radaroutput[6][j,k] > minlatwrf and radaroutput[6][j,k] < maxlatwrf and radaroutput[0][j,k] < 0:
					datawriter.writerow([round(radaroutput[5][j,k],3), round(radaroutput[6][j,k],3), round(radaroutput[0][j,k],1), -888888.000, round(radaroutput[1][j,k],3)])
				elif radaroutput[0][j,k] and not(radaroutput[1][j,k]) and radaroutput[5][j,k] > minlonwrf and radaroutput[5][j,k] < maxlonwrf and \
					 radaroutput[6][j,k] > minlatwrf and radaroutput[6][j,k] < maxlatwrf:
					datawriter.writerow([round(radaroutput[5][j,k],3), round(radaroutput[6][j,k],3), round(radaroutput[0][j,k],1), round(radaroutput[0][j,k],3), -888888.000])

def write_obradar(wrfdainput,radarstatdata,radartempfile,banneddata,dbzthres,site,date):
	print '\nWritting ob.radar'
	alldata = np.genfromtxt(radartempfile,delimiter=',')
	datacnt = []
	for i in xrange(len(alldata)):
		if not(alldata[i,4] in banneddata) and alldata[i,4] > dbzthres:
			datacnt.append(alldata[i])
	datacnt = np.asarray(datacnt)

	# calculate total sounding
	totsound = len(datacnt)

	# calculate total height
	heightcnt = []
	for i in xrange(len(datacnt)):
		if not(alldata[i,2] in heightcnt):
			heightcnt.append(alldata[i,2])

	heightcnt = np.asarray(heightcnt)
	totheight = len(heightcnt)

	# get radar stat
	radarstat = np.genfromtxt(radarstatdata, delimiter=',')

	#write header
	datawriter = open(wrfdainput, 'w')
	datawriter.write('TOTAL NUMBER =  1\n')
	datawriter.write('#-----------------#\n')
	datawriter.write('\n')

	datawriter.write('RADAR	  '+site+' '*(13-len(str(radarstat[0])))+str(radarstat[0])+' '*(10-len(str(radarstat[1])))+\
	str(radarstat[1])+' '*(10-len(str(radarstat[2])))+str(radarstat[2])+'  '+\
	date+' '*(6-len(str(totsound)))+str(totsound)+' '*(6-len(str(totheight)))+str(totheight)+'\n')

	datawriter.write('#-------------------------------------------------------------------------------#\n')
	datawriter.write('\n')
	for i in xrange(totsound):
		datawriter.write('FM-128 RADAR   '+date+' '*(14-len(str(datacnt[i,1])))+str(datacnt[i,1])+' '*(14-len(str(datacnt[i,0])))+\
		str(datacnt[i,0])+' '*(10-len(str(radarstat[2])))+str(radarstat[2])+'	   1\n')
		
		datawriter.write(' '*(15-len(str(datacnt[i,2])))+str(datacnt[i,2])+' '*(12-len(str(datacnt[i,4])))+\
		str(datacnt[i,4])+'   0		   1'+' '*(14-len(str(datacnt[i,3])))+str(datacnt[i,3])+'   0		   1'+'\n')

	
''' Main program '''
# Set parameter
# Scan total can be defined or default of total radar elevation
# For gematronik radar, sometimes very high elevation gives an error
os.environ['WRADLIB_DATA']='/home/genomexyz/radarvol_to_wrfda'
wrffile = 'WRF/wrfinput_d01'
radarpath = 'VOL'
radarstatdata = 'out/radar.stat'
radartempfile = 'out/radar.temp'
wrfdainput = 'out/gemaob.radar'
banneddata = np.asarray([-32.118,-888888.0]);dbzthres = 0
rest = 0.001;res_radar = 1000.;rmax=250000.
radartime=datetime(2019,1,20,0,0)
listpcdfpath=rungema(radartime,radarpath)
print('total adalah', listpcdfpath)
f=wrl.util.get_wradlib_data_file(listpcdfpath[0])
raw=wrl.io.read_Rainbow(f)
#scantotal=int(np.size(raw['volume']['scan']['slice']))
scantotal=9
site='jaktes'
date = '2017-06-15_00:12:00'

write_radartemp(wrffile,radartempfile,radarstatdata,rest,scantotal,radartime,radarpath)
write_obradar(wrfdainput,radarstatdata,radartempfile,banneddata,dbzthres,site,date)
