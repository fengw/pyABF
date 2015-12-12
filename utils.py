#!/usr/bin/env python 
""" 
Utilities for CyberShake related

""" 

# General libraries
import os,sys, shutil
import glob, time
import numpy as np 

from scipy import stats 

import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D 

from time import sleep 
import MySQLdb as mdb

# =============
# NGA package with Directivity
# =============
from pynga import *
from pynga.CB08 import *
from pynga.BA08 import *
from pynga.CY08 import *
from pynga.AS08 import *
from pynga.SC08 import *
from pynga.utils import *

# ====================
# use my_util tools
# ====================
from my_util.geopy import *
from my_util.signal import *
from my_util.numer import *
from my_util.indpy import *
from my_util.functools import *
from my_util.image import * 


def GetKey4(group_i,index_j,rup_var_id, Ti):
    """
    Generate key for rupture variations (related to hypocenter group)
    """
    return '(%s,%s,%s,%s)'%(group_i,index_j,rup_var_id,Ti)

def GetKey2(ih, Ti):
    """
    Generate key for ihypo and T_i (related to hypocenter group)
    """
    return '(%s,%s)'%(ih,Ti)

# communicate with HPC to download files 
# This can be used for other platforms
def hpc_scp(srcfile, dst, machine='hpc', opt='-r', inverse=False):
    # transfer file between local machine and hpcc
    # inverse = False: Download from hpc 
    # inverse = True:  Upload to hpc 

    usrnam = 'fengw@hpc-login1.usc.edu:'
    scp = 'scp'
    
    if inverse == False:
	if not os.path.exists( dst ):
	    os.mkdir(dst)
	prefix = ' '.join( (scp,opt,usrnam) )
	cmd = ' '.join( (prefix+srcfile, dst+os.sep) )
	#print cmd
	os.system( cmd )
    else:
	prefix1 = ' '.join((scp,opt))
	prefix2 = usrnam
	cmd = ' '.join( (prefix1, srcfile, prefix2 + dst + os.sep) )
	#print cmd
	os.system( cmd )


# ================================================
# waveform extraction 
# ================================================
def ExtractRups(cursor,sids,rids,erf_id=35,rup_scenario_id=3):
    # can be modified
    rups_info = {}   

    if rup_scenario_id == 4:
	rup_scenario_id = 3

    for irup in xrange( len(sids) ):
	sid = sids[irup]
	rid = rids[irup]
	sr_key = '%s_%s'%(sid,rid)
	
	# Ruptures
	query = "select * from %s where %s = %s and %s = %s and %s = %s"%\
		('Ruptures','ERF_ID',erf_id,'Source_ID',sid,'Rupture_ID',rid)
	cursor.execute( query )       # run query
	row_rup = cursor.fetchone()
	nvar = len(row_rup)   
	if nvar == 0:
	    print 'There is no rupture for (erf_id,rup_var_id,sid,rid) = (%s,%s,%s,%s)\n'%(erf_id,rup_scenario_id,sid,rid)
	    return rups_info
	else:
	    pass

	# Rupture Variances (hypocenter and slip)
	query = "select * from %s where ERF_ID = %s and Rup_Var_Scenario_ID = %s and Source_ID = %s and Rupture_ID = %s"%('Rupture_Variations',erf_id,rup_scenario_id,sid,rid)
	cursor.execute( query )       # run query
	row_rup_var = cursor.fetchall()  
	nvar = len(row_rup_var)   
	if nvar == 0:
	    print 'There is no rupture variations for (erf_id,rup_var_id,sid,rid) = (%s,%s,%s,%s)\n'%(erf_id,rup_scenario_id,sid,rid)
	    return rups_info
	else:
	    pass
	

	# regroup the hypocenter and slip variation
	for k in xrange( len(row_rup_var) ):
	    slip_hypo = row_rup_var[k][5].strip().split('-')
	    tmp = int(slip_hypo[2][1:])
	    if tmp != k:
		break
	    Nh = tmp + 1 
	Ns = nvar/Nh
	rups_info[sr_key] = [Nh,Ns]
	
    return rups_info



def ExtractListGen(cursor, Sites, Sources, Ruptures, fid, rup_model_ids=(35,5,3,1), ih=None, islip=None, bb=None ):
    
    erf_id, sgt_id, rup_scenario_id, vel_id = rup_model_ids
    rups_info = ExtractRups( cursor, Sources, Ruptures, erf_id=erf_id, rup_scenario_id = rup_scenario_id )
    
    SourceRups = []
    for stanam in Sites:

	# Get site id for further selection
	query = "select * from %s where %s = '%s'"%('CyberShake_Sites','CS_Short_Name',stanam)
	cursor.execute( query )       # run query
	row_sites = cursor.fetchone()  

	# fetch run id (select verified runs)
	# and also Max_Frequency and Low_Frequency_Cutoff to get deterministic and broadband 
	if bb == None:
	    query = "select * from %s where ERF_ID = %s and \
		     SGT_Variation_ID = %s and Rup_Var_Scenario_ID = %s and Velocity_Model_ID= %s and \
		     Site_ID = %s and Status = '%s'"%('CyberShake_Runs',erf_id,sgt_id, rup_scenario_id, vel_id, row_sites[0],'Verified')
	else:
	    query = "select * from %s where ERF_ID = %s and \
		     SGT_Variation_ID = %s and Rup_Var_Scenario_ID = %s and Velocity_Model_ID= %s and \
		     Site_ID = %s and Status = '%s' and Max_Frequency = %s"%('CyberShake_Runs',erf_id,sgt_id, rup_scenario_id, vel_id, row_sites[0],'Verified', bb)
	cursor.execute( query )       # run query
	try:
	    row_run = cursor.fetchone()
	    Nrun = len( row_run )
	    site_run_id = row_run[0]
	except:
	    print 'There is no verified run for parameter set (erfid,sgtid,rup_var_scenario_id,vel_id,site_name)=(%s,%s,%s,%s,%s)\n'%(erf_id,sgt_id,rup_scenario_id,vel_id,stanam)
            continue

	RunID = '%s'%site_run_id
	sids = []; rids = []
	for irup in xrange( len(Sources) ):
	    sid = Sources[irup] 
	    rid = Ruptures[irup]

	    #  CyberShake_Site_Ruptures selection
	    query = "select * from %s where CS_Site_ID = %s and Source_ID = %s and Rupture_ID = %s"%('CyberShake_Site_Ruptures',row_sites[0],sid,rid)
	    cursor.execute( query )       # run query
	    row_site_rup = cursor.fetchall()  
	    nvar = len(row_site_rup)
	    if nvar == 0:
		print 'There is no ruptures for (Site_name, Site_ID, run_id, sid,rid) = (%s, %s,%s,%s,%s)\n'%(row_sites[2],row_sites[0],site_run_id,sid,rid)
		continue
            
	    sr_key = '%s_%s'%(sid,rid)
	    
	    sids.append(sid)
	    rids.append(rid)

	    Nh,Ns = rups_info[sr_key] 
	    Nvar = Nh*Ns 

            if ih == None and islip == None:
		# choose the center hypocenter location and slip = 1, 2 for this hypocenter? 
		ih = int(Nh/2)    # centroid biased 
		
		# since for one hypocenter, there are two more slip distributions generated (Ns = 2*Nh for each hypocenter, but we only take the two generated for the hypocenter)
		islip1 = 2*ih 
		ivar1 = HypoSlip2RupVar((ih,islip1), Nh, Ns)
		islip2 = 2*ih + 1
		ivar2 = HypoSlip2RupVar((ih,islip2), Nh, Ns)
		#print sr_key, Nh, Ns, ih, ivar1, ivar2

		fid.write( '%s   %s   %s   %s   %s\n'%(stanam,RunID,sid,rid,ivar1) )
		fid.write( '%s   %s   %s   %s   %s\n'%(stanam,RunID,sid,rid,ivar2) )
	    else: 
		ivar = HypoSlip2RupVar((ih,islip), Nh, Ns)
		fid.write( '%s   %s   %s   %s   %s\n'%(stanam,RunID,sid,rid,ivar) )

        
	SourceRups.append( [sids,rids] ) 

    return Sites, SourceRups


# Correlation analysis (dataset of CyberShake)
def PlotDataSet(Inpth, fname, Plotpth, plotname, savetype):
    
    # read in flatdata
    filename = Inpth + fname
    data = np.loadtxt( filename )
    Mws = data[:,2]
    rakes = data[:,3]
    Rjbs = data[:,4]
    Vs30s = data[:,5]
    
    Nrecord = len( Mws )

    # plot 
    fig = init_fig( num=1, figsize=(14,10), dpi=100 )
    fig.text( 0.2, 0.95, 'Total number of seismograms: %s'%Nrecord )
    axs = init_subaxes( fig, subs=(2,2), basic=(0.7,0.7,0.7,0.7) )
    
    ax = fig.add_axes( axs[0] )
    ax.hist( Mws, bins = 30 )
    ax.set_xlabel( '$M_w$' )
    #ax.set_ylabel( '# of records' )

    ax = fig.add_axes( axs[1] )
    group,group_names = RakeBin(rakes)
    Ng = len(group.keys())
    ind = np.arange( Ng )+1; width = 0.5
    Nrake = []; clr = []; Names = []; rects = []
    for ikey, key in enumerate( group.keys() ):
	Nrake = len( group[key] ) 
	#Names.append( group_names[key][1] )
	Names.append( key )
	rect0 = ax.bar(ind[ikey], Nrake, width, color = group_names[key][0])
	rects.append( rect0[0] )
    ax.set_xticks(ind+width/2)
    ax.set_xticklabels( tuple( Names ) )
    #ax.set_ylabel( '# of records' )
    ax.set_title( 'Faulting' )

    ax = fig.add_axes( axs[2] )
    ax.hist( Rjbs, bins = 30 )
    ax.set_xlabel( '$R_{jb}$' )
    #ax.set_ylabel( '# of records' )

    ax = fig.add_axes( axs[3] )
    group,group_names = Vs30Bin(Vs30s)
    Ng = len(group.keys())
    ind = np.arange( Ng )+1; width = 0.5
    Nrake = []; clr = []; Names = []; rects = []
    for ikey, key in enumerate( group.keys() ):
	Nrake = len( group[key] ) 
	#Names.append( group_names[key][1] )
	Names.append( key )
	rect0 = ax.bar(ind[ikey], Nrake, width, color = group_names[key][0])
	rects.append( rect0[0] )
    ax.set_xticks(ind+width/2)
    ax.set_xticklabels( tuple( Names ) )
    #ax.set_ylabel( '# of records' )
    ax.set_title( '$Vs_{30}$' )

    fig.savefig( Plotpth + plotname, format = savetype, dpi = 300 )
    #fig.savefig( Plotpth + plotname, format = savetype )


# ================================================
# Calculate the spectral intensities (for correlation analysis)
# For validation
# ================================================
def rspectra( a, dt, period, dr = 0.05, chan='a', pseudo=True ):
    # use irc to compute
    n = len(a)  # length of the input series
    pr = period
    wn = 2*np.pi/pr

    if chan == 'v':
	u,a = v2ua( a, dt )    # velocity to acceleration
    elif chan == 'u':
	v,a = u2va( a, dt )    # displacement ot acceleration 
    
    # impulse response convolution (irc)
    wd = wn*np.sqrt(1-dr*dr)
    t = dt*np.arange(n)

    im0 = -1./wd * np.exp( -dr*wn*t ) * np.sin( wd*t )   # implusive response of oscillator
    Nfft,ntmp = n2pow(n)

    # convolution in frequency domain (using FFT)                                
    x = np.fft.ifft( np.fft.fft( im0,Nfft )*np.fft.fft( a,Nfft ) ) * dt                    
    
    t0 = 0.0
    nt = n
    
    # header info
    hd = [pr, dr, dt, t0, nt]
    
    if not pseudo:
	im1 = np.exp( -dr*wn*t ) * ( dr*wn/wd * np.sin( wd*t ) - np.cos( wd*t ) )
	im2 = np.exp( -dr*wn*t ) * ( (wd**2-(dr*wn)**2)/wd * np.sin( wd*t ) + 2.*dr*wn*np.cos( wd*t ) )
	x1 = np.fft.ifft( np.fft.fft( im1,Nfft )*np.fft.fft( a,Nfft ) ) * dt                    
	x2 = np.fft.ifft( np.fft.fft( im2,Nfft )*np.fft.fft( a,Nfft ) ) * dt                    
        return hd, x, x1, x2
    else: 
	return hd, x


def calc_SI(input,dt,h1,h2,dr=0.05,chan='a'):
    # Compute spectra intensity (integral over time window, or duration)
    if chan == 'v':
        v = input
	u,a = v2ua( input, dt )    # velocity to acceleration
    elif chan == 'u':
	u = input
	v,a = u2va( input, dt )    # displacement to acceleration 
    elif chan == 'a':
	a = input
	u,v = a2uv( input, dt )    # acceleration to displacement and velocity
    else:
	print 'chan should be u, v, or a'
	raise ValueError
    chan = 'a'

    if h1 <= h2: 
	print 'h1 should be much larger than h2'
	raise ValueError

    periods = np.arange( 0.1, 2.5+h1, h1 )
    PSV_T = []; SA_T = []
    for it in xrange( len(periods) ):
	period = periods[it]
        wn = 2.*np.pi / period

	hd, x = rspectra( a, dt, period, dr=dr, chan=chan )
	Sd,ind_d = maxi(abs(x))
	
	# we use pseudo PSV and PSA
	Sv = wn*Sd
	Sa = wn*Sv
        
	SA_T.append(Sa)
	PSV_T.append(Sv)
    
    periods_new = np.arange( 0.1, 2.5+h2, h2 )
    PSV = interp1d( periods, PSV_T, periods_new, 0, 0 )
    SI = sum( PSV ) * h1 

    periods_new = np.arange( 0.1, 0.5+h2, h2 )
    SA = interp1d( periods, SA_T, periods_new, 0, 0 )
    ASI = sum( SA ) * h1
    return SI, ASI 



def calc_IMs( input, dt, periods, SI_ASI=None, dr=0.05, chan='a'):
    """
    # Compute IMs for given accelation
    periods: at which SA will be computed
    SI_ASI = [flag,h1,h2], if flag != 0:
    h1 = 0.05; h2 = 0.01
    h1, h2 : for SI and ASI computation ( h1 > h2 ), interpolation
    """
    input = np.array( input, 'f' )
    if chan == 'v':
        v = input
	u,a = v2ua( input, dt )    # velocity to acceleration
    elif chan == 'u':
	u = input
	v,a = u2va( input, dt )    # displacement to acceleration 
    elif chan == 'a':
	a = input
	u,v = a2uv( input, dt )    # acceleration to displacement and velocity
    else:
	print 'chan should be u, v, or a'
	raise ValueError
    chan = 'a'
    output = {}

    # PGD, PGV, PGA, SI, and ASI (Now a is acceleration)
    pgd = max( abs(u) )
    pgv = max( abs(v) )
    pga = max( abs(a) )

    output['PGD'] = pgd
    output['PGV'] = pgv
    output['PGA'] = pga

    try: 
	Np = len(periods)
    except: 
	periods = [periods]
	Np = len(periods)

    # SD, SV and SA
    output['SA'] = []
    for it in xrange( len(periods) ):
	period = periods[it]
        wn = 2.*np.pi / period
       
        hd, x = rspectra( a, dt, period, dr=dr, chan=chan )
        t0 = hd[3]

	Sd,ind_d = maxi(abs(x))
	Sv = wn*Sd
	Sa = wn*Sv
	output['SA'].append(Sa)    # SA in cm/s^2
    
    # SI and ASI calculating
    if SI_ASI != None:
	h1,h2 = SI_ASI
	SI, ASI = calc_SI(a, dt, h1, h2)

	output['SI'] = SI
	output['ASI'] = ASI

    return output


def CalcGMRotIpp(input1,input2,osc1,osc2,periods,pp=-1, GeoMean=True, RotD=False):
    
    input1 = np.array( input1 )
    input2 = np.array( input2 )
    osc1 = np.array( osc1 )
    osc2 = np.array( osc2 )
    Np = len(osc1) 
    
    # outputs:
    output = {}; cita_min = {}
    
    # compute as-record PGA (what about the rotation-independent PGA?)
    if GeoMean:
	output['PGA'] = np.sqrt( max(abs(input1))*max(abs(input2)) )
    else: 
	output['PGA'] = [max(abs(input1)), max(abs(input2))]

    cita_min['PGA'] = -1

    if pp < 0:
	# compute as-recorded SA
	GMsa = []
	for ip in xrange( Np ):
	    SA1p = max(abs(osc1[ip,:]))
	    SA2p = max(abs(osc2[ip,:]))
	    if GeoMean: 
		GMsa.append(np.sqrt(SA1p*SA2p))
	    else: 
		GMsa.append( [SA1p,SA2p] ) # x and y (determined by the computing box of CyberShake)
	
	output['SA'] = GMsa
	cita_min['SA'] = -1
     
    else:
	# do the rotation
	citas = range( 0, 91 )   # rotation angle
	Nc = len(citas)
	
	GM0 = []; GM1 = [] 
	for ip in xrange( Np ):
	    
	    GMsa = []
	    for ic in xrange( Nc ):
		cita = citas[ic] * np.pi / 180.  # degree to radius
		c = np.cos(cita); s = np.sin(cita)
		osc1p = c * osc1[ip,:] + s * osc2[ip,:]
		osc2p = -s * osc1[ip,:] + c * osc2[ip,:]
		SA1p = max(abs(osc1p))    # get SA value
		SA2p = max(abs(osc2p))    
		GMsa.append(np.sqrt(SA1p*SA2p))
	    GM0.append( GMsa )  # GM(cita,T)
	    
	    # find GMRotDpp(T)
	    GM_rank,index = sorti( GMsa )
	    GM1.append( np.percentile(GM_rank, pp) ) 

        if not RotD: 
	    
	    GM0 = np.array( GM0 )
	    GM1 = np.array( GM1 )
	    
	    # Compute penalty(cita)
	    penalty = 10000.
	    for ic in xrange( Nc ):
		tmp = 1./Np * sum( (GM0[:,ic]/GM1[:]-1)**2 )
		if tmp <= penalty:
		    penalty = tmp
		    index_cita = ic
	    cita_min0 = citas[index_cita]   # find the minimum cita

	    output['SA'] = list( GM0[:,index_cita] )
	    cita_min['SA'] = citas[index_cita]

	else: 
	    output['SA'] = GM1
	    cita_min['SA'] = -1

    return output, cita_min


def CalcRotIppIMs( inputx, inputy, dt, periods, dr=0.05, chan='a', pp=-1, GeoMean = True, RotD=False ):
    """
    Compute Rot-Independent IM (Geometrical mean) 
    Following Boore 2006

    periods: at which SA will be computed
    default gives the geometrical mean or the Rot independent geometrical mean 
    """
    inputx = np.array( inputx, 'f' )
    inputy = np.array( inputy, 'f' )

    if chan == 'v':
        vx = inputx
	vy = inputy
	ux,acx = v2ua( inputx, dt )    # velocity to acceleration
	uy,acy = v2ua( inputy, dt )    # velocity to acceleration
    elif chan == 'u':
	ux = inputx
	vx,acx = u2va( inputx, dt )    # displacement to acceleration 
	uy = inputy
	vy,acy = u2va( inputy, dt )    # displacement to acceleration 
    elif chan == 'a':
	acx = inputx
	ux,vx = a2uv( inputx, dt )    # acceleration to displacement and velocity
	acy = inputy
	uy,vy = a2uv( inputy, dt )    # acceleration to displacement and velocity
    else:
	print 'chan should be u, v, or a'
	raise ValueError

    chan = 'a'
    
    try: 
	Np = len(periods)
    except: 
	periods = [periods]
	Np = len(periods)
    
    # response(w0,t) as calculated
    osc1 = []; osc2 = []
    for it in xrange( len(periods) ):
	period = periods[it]
	hd,x,x1,osc10 = rspectra( acx, dt, period, dr=dr, chan=chan, pseudo=False )
	hd,x,x1,osc20 = rspectra( acy, dt, period, dr=dr, chan=chan, pseudo=False )
	osc1.append( osc10 )
	osc2.append( osc20 )
    
    output, cita_min = CalcGMRotIpp(acx, acy, osc1, osc2, periods, pp=pp, GeoMean=GeoMean, RotD=RotD)
    
    # the unit of SA is depending on the input series !!!

    # output has key: PGA, SA
    return output



# correlation coefficient (for validation)
# follow Bradley 2011
def corrcoef_sd(rho,Nsample,alpha):
    # http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
    ci = []

    # table (just for 95%, and 90% confidence interval)
    mu2 = {'0.05': 1.96, '0.10':1.6449} 
    SE = 1./np.sqrt( Nsample-3 )
    sigma = mu2['%s'%('%3.2f'%alpha)]*SE
    for ir in xrange( len(rho) ):
	a = np.arctanh( rho[ir] )
	ci.append( [np.tanh(a-sigma), np.tanh(a+sigma)] )
    return ci


def rho_model( T, model='BB', Tref = None ):
    if model == 'BB':
	# BB 2011 model
	if T < 0.01 or T > 10.0:
	    print 'period exceeds the limit [0.01,10.0]'
	    raise ValueError
	else:
	    if 0.01<=T<0.2:
		an = 1.00
		bn = 0.895
		cn = 0.06
		dn = 1.6
	    else:
		an = 0.97
		bn = 0.25
		cn = 0.8
		dn = 0.8
	    return 0.5*(an+bn) - 0.5*(an-bn) * np.tanh( dn*np.log(T/cn) )
    
    if model == 'JB01':
	if T <0.05 and T > 5:
	    print 'period exceeds the limit [0.01,10.0]'
	    raise ValueError
	else:
	    return (0.5-0.127*np.log(T))*(0.05<=T<0.11) + \
		    (0.968+0.085*np.log(T))*(0.11<=T<0.25) + \
		    (0.568-0.204*np.log(T))*(0.25<=T<=5)
    
    if model == 'JJ08':
	# Jaker and Jayaram 2008 Earthquake Spectra
	# This is for the total residual (not for the inter- or intra- residual)
	if Tref == None:
	    print 'Tref should have valid value...'
	    raise ValueError
	else:
	    Tmax = max( T, Tref )
	    Tmin = min( T, Tref )
	    C1 = 1 - np.cos( np.pi/2. - 0.366*np.log(Tmax/max(Tmin,0.109)) )
	    if Tmax >= 0.2: 
		C2 = 0.0 
	    else: 
		C2 = (1-0.105*(1-1/(1+np.exp(100*Tmax-5)))*(Tmax-Tmin)/(Tmax-0.0099) )
	    C3 = C2* (Tmax<0.109) + C1 * (Tmax>=0.109)
	    C4 = C1 + 0.5*(np.sqrt(C3)-C3) * (1+np.cos(np.pi*Tmin/0.109))
	    if Tmax <= 0.109:
		return C2
	    elif Tmin > 0.109:
		return C1
	    elif Tmax < 0.2:
		return min(C2,C4)
	    else:
		return C4




# ==================================
# CyberShake NGA related
# ==================================
def cybershk_sites(sitedata):
    """
    All CyberShake Sites information (Vs30, Z2.5, Z1.0 etc.)
    sitedata give the full path of the site file
    """
    sites = {}
    
    # required header of the sitedata file
    keys = 'lat','lon', 'Vs30_Wills_2006', 'Vs30_Wald_2007', 'Z2.5', 'Z1.0'  # Z2.5,Z1.0 all in km !!
    lines = open( sitedata, 'r' ).readlines()
    for i in range( 1, len(lines) ):
	spl = lines[i].strip().split()
	name = spl[0]
	sites[name] = {}
	for ikey in xrange( len(keys) ):
	    sites[name][keys[ikey]] = spl[ikey+1]
    return sites


# Source related 
def RupSelect(sids,cybershk_database_password, erf_id=35):
    """
    Select from UCERF and CyberShake database (by Region and Magnitude)
    Input: 
        Select Dictionary (contain all information need to select ruptures)
	erf_id = 35 (CyberShake database)
    Output: 
	sids
	rids
    """
    rids = []; sids0 = []
    
    hostname = 'focal.usc.edu'  # host                            
    username = 'cybershk_ro'    # usrname                         
    password = cybershk_database_password    # password                        
    database = 'CyberShake'     # Database name 
    db = mdb.connect( host = hostname, user = username, \
		 passwd = password, db = database )
    cursor = db.cursor()
    
    # constains on Ruptures Table in CyberShake to Select Ruptures
    Nsid = len(sids)

    for isid in xrange( Nsid ):
	sid = sids[isid]
	query = 'select * from Ruptures where ERF_ID=%s and Source_ID=%s'%(erf_id,sid)
	cursor.execute( query )
	rows = cursor.fetchall()
	Nrup = len(rows)
	rid0 = []; sid0 = []; Mws = []
	Area = rows[-1][10]   # preserve the rupture area!
	for irup in xrange( Nrup ):
	    # deal with different cases
	    if rows[irup][10] == Area:
		if not rows[irup][5] in Mws:
		    rid0.append( rows[irup][2] )
		    Mws.append( rows[irup][5] )
		    sid0.append( sid ) 
		else: 
		    continue 
	sids0.append(sid0)
	rids.append( rid0 )
    
    cursor.close()
    db.close()

    return sids0, rids

# =============================
# Extract Rupture Info from CyberShake Database 
# =============================
def rup_gen(cursor, sid, rups_info, fid_stdout, erf_id=35, rup_scenario_id=3, hypo_info=True):
    """
    Extract src detail info from UCERF 
    Given one source, get all source infomation: 
	stress drop variability (Mw)
	hypocenter location  (Nh)
	slip distribution  (Nvar)
    """
    # correction (due to database error)
    # for example: select * from Rupture_Variations where ERF_ID=35 and Source_ID=10 and Rupture_ID=3 and Rup_Var_Scenario_ID=4
    # which gives all the same hypocenter locations (just this)
    if rup_scenario_id == 4:
	rup_scenario_id = 3

    # Ruptures
    query = "select * from Ruptures where ERF_ID = %s and Source_ID = %s"%(erf_id,sid)
    cursor.execute( query )       # run query
    row_rup = cursor.fetchall()
    Nrup = len(row_rup)   
    if Nrup == 0:
	fid_stdout.write('There is no SourceID=%s for erf_id = %s\n'%(sid, erf_id))
	return rups_info
    else:
	pass
    
    # fault information (use the last one)
    SourceName = row_rup[-1][3]   # section name
    Nrow = row_rup[-1][8]  # vertical points on the fault plane
    Ncol = row_rup[-1][9]  # horizontal poinst on the fault plane 
    Mws = []; rids = []
    Area = row_rup[-1][10]   # preserve the rupture area!
    for irup in xrange( Nrup ):
	# deal with different cases
	if row_rup[irup][10] == Area:
	    if not row_rup[irup][5] in Mws:
		rids.append( row_rup[irup][2] )
		Mws.append( row_rup[irup][5] )
	else: 
	    continue

    # Fault surface and trace line
    # use the last Rupture ID to get the Full length of fault segment (Type A)
    query = "select * from Points where ERF_ID= %s and Source_ID= %s and Rupture_ID=%s"%(erf_id,sid, Nrup-1)   
    cursor.execute( query )      
    row_point = cursor.fetchall()
    Npoints = len(row_point)   
    if Npoints == 0:
	fid_stdout.write('There is no faults surface points for (erf_id,sid) = (%s,%s)\n'%(erf_id,sid))
	return rups_info
    else:
	pass
    row_point = np.array( row_point ) 

    Ztor = row_point[0,6]
    Zbom = row_point[-1,6] 
    rake = np.mean( row_point[:,7] )
    dip = np.mean( row_point[:,8] )
    
    # surface points
    index_surf = (row_point[:,6] == Ztor).nonzero()[0]
    lon = row_point[index_surf,5]
    lat = row_point[index_surf,4]
    strike = row_point[index_surf,9]
    
    # decimate the surface points
    dec = 5 
    lon = lon[::dec].tolist()
    lat = lat[::dec].tolist()
    strike = strike[::dec].tolist() 

    # all fault points
    lonf = row_point[:,5].tolist() 
    latf = row_point[:,4].tolist()
    depf = row_point[:,6].tolist() 

    if hypo_info == True:
	# Rupture Variances (hypocenter and slip)
	query = "select * from %s where ERF_ID = %s and Rup_Var_Scenario_ID = %s and Source_ID = %s and Rupture_ID = %s"%('Rupture_Variations',erf_id,rup_scenario_id,sid,Nrup-1)
	cursor.execute( query )       # run query
	row_rup_var = cursor.fetchall()  
	nvar = len(row_rup_var)   
	if nvar == 0:
	    fid_stdout.write('There is no rupture variations for (erf_id,rup_var_id,sid) = (%s,%s,%s)\n'%(erf_id,rup_scenario_id,sid))
	    return rups_info
	else:
	    pass
	
	# regroup the hypocenter and slip variation 
	tmp = row_rup_var[-1][5].strip().split('-')
	Nh = int(tmp[2][1:]) + 1
	Nf = int(tmp[1][1:]) + 1
	hypo_loc = {}
	
	# get hypocenter locations
	for ih in xrange( Nh ):
	    # attention: for the new rupture generator, the rupture variations are the combination of 
	    # s and h. for example, for the same h id, but different s id, the hypocenter locations are not the same due to 
	    # the randomization process, for the old rupture generator, for the same h id, the location and
	    # the depth are the same regardless the s id. 
	    rup_var_id = HypoSlip2RupVar( (ih,0), Nh, Nf )
	    hypo_key = 'h%4.4i'%ih
	    hlon = row_rup_var[rup_var_id][7]
	    hlat = row_rup_var[rup_var_id][6]
	    hdep = row_rup_var[rup_var_id][8]
	    hypo_loc[hypo_key] = [hlon,hlat,hdep]
	rups_info['hypo_slip'] = [Nh, Nf, hypo_loc]
    else:
	rups_info['hypo_slip'] = []

    # Save fault geometry info to dictionary for further use
    rups_info['name'] = SourceName
    rups_info['MR'] = [Mws,rids,Nrow,Ncol,rake,dip,Ztor,Zbom]  # magnitude and points on the fault, averege rake, dip, Ztor, and Zbom
    rups_info['surface'] = [lon,lat,strike]    # fault plane surface projection
    rups_info['fault'] = [lonf, latf, depf]    # fault plane
    
    return rups_info


def HypoSlip2RupVar(index, Nh, Nf, inverse=False):
    """
    Relation between rupture variation ID and hypocenter and slip distribution ID 
    """
    if not inverse: 
	ih, islip = index 
	ivar = islip * Nh + ih
    else: 
	for islip in xrange( Nf ):
	    for ih in xrange( Nh ):
		if index == islip * Nh + ih: 
		    ivar = (ih, islip)
		    break
    return ivar 


# ========================================
# Extract IM data from CyberShake databas
# ========================================
def im_gen(cursor, sid, rid, stanam, \
	   sites_info, Ek,\
	   fid_stdout, rup_model_ids=(35,5,3,1), bb=None, site_only=False):
    """
    Give CyberShake site (stanam), UCERF2.0 Source (sid)

    rup_model_ids = (erf_id, sgt_id, rup_var_scenario_id, velocity_id)
    if site_only == True, then just save the site info, not other
    bb = None: get the deterministic runs
    bb = 0: get all verified simulations
    otherwise Max_Frequency = bb
    """
    
    # CyberShake study
    erf_id, sgt_id, rup_scenario_id, vel_id = rup_model_ids

    # Get site id for further selection
    query = "select * from %s where %s = '%s'"%('CyberShake_Sites','CS_Short_Name',stanam)
    cursor.execute( query )       # run query
    row_sites = cursor.fetchone()  

    # fetch run id (select verified runs)
    # and also Max_Frequency and Low_Frequency_Cutoff to get deterministic and broadband 
    if bb == None:
	# just use the low frequency runs
	query = "select * from %s where ERF_ID = %s and \
		 SGT_Variation_ID = %s and Rup_Var_Scenario_ID = %s and Velocity_Model_ID= %s and \
		 Site_ID = %s and Status = '%s' and Max_Frequency IS %s"%('CyberShake_Runs',erf_id,sgt_id, rup_scenario_id, vel_id, row_sites[0],'Verified', 'NULL')
    else:
	if bb != 0:
	    # specify the max frequency 
	    query = "select * from %s where ERF_ID = %s and \
		     SGT_Variation_ID = %s and Rup_Var_Scenario_ID = %s and Velocity_Model_ID= %s and \
		     Site_ID = %s and Status = '%s' and Max_Frequency = %s"%('CyberShake_Runs',erf_id,sgt_id, rup_scenario_id, vel_id, row_sites[0],'Verified', bb)
	else:
	    # include all (there might be repeated sites, this extraction will overwrite the existing sites)
	    query = "select * from %s where ERF_ID = %s and \
		     SGT_Variation_ID = %s and Rup_Var_Scenario_ID = %s and Velocity_Model_ID= %s and \
		     Site_ID = %s and Status = '%s'"%('CyberShake_Runs',erf_id,sgt_id, rup_scenario_id, vel_id, row_sites[0],'Verified')

    cursor.execute( query )       # run query
    try:
	row_run = cursor.fetchone()
	Nrun = len( row_run )
	site_run_id = row_run[0]
    except:
	fid_stdout.write('There is no verified run for parameter set (erfid,sgtid,rup_var_scenario_id,vel_id,site_name)=(%s,%s,%s,%s,%s)\n'%(erf_id,sgt_id,rup_scenario_id,vel_id,stanam))
	return sites_info, Ek

    # PeakAmplitude selection (rupture_ID, some sources donot use part of their rupture))
    query = "select Source_ID, Rupture_ID from %s where Run_ID = %s and Source_ID = %s and Rup_Var_ID = %s"%('PeakAmplitudes',site_run_id,sid,0)
    cursor.execute( query )       # run query
    row_Sa = cursor.fetchall()  
    Nrup = len(row_Sa)
    if Nrup == 0:
	fid_stdout.write('There is SA for (Site_name, Site_ID, run_id, sid) = (%s, %s,%s,%s)\n'%(row_sites[2],row_sites[0],site_run_id,sid))
	return sites_info, Ek
    else:
	pass

    # deal with the problem that not all ruptures for one source are used in the simulation
    
    # if all selection rules above are satisfied, then leave the sites and do the following stuff.
    # site info dict (create if run_id for the site exists)
    # From now on, stanam is used as the combination of (SiteName, RunID, SourceID)
    sites_info[stanam] = {}
    sites_info[stanam]['id'] = row_sites[0]  
    sites_info[stanam]['name'] = row_sites[2]    # short name of Sites
    sites_info[stanam]['lat'] = row_sites[3]
    sites_info[stanam]['lon'] = row_sites[4]
    sites_info[stanam]['type'] = row_sites[5]
    sites_info[stanam]['run_id'] = site_run_id

    if site_only == False:
	
	# get hypo slip index 
	query = "select * from %s where ERF_ID = %s and Rup_Var_Scenario_ID = %s and Source_ID = %s and Rupture_ID = %s"%('Rupture_Variations',erf_id,rup_scenario_id,sid,rid)
	cursor.execute( query )       # run query
	row_rup_var = cursor.fetchall()  
	tmp = row_rup_var[-1][5].strip().split('-')
	Nh = int(tmp[2][1:]) + 1
	Nf = int(tmp[1][1:]) + 1
	
	# Get all periods (from IM_Types table based on IM_Type_ID)
	query = "select * from %s"%('IM_Types')
	cursor.execute( query )       # run query
	row_imtype = cursor.fetchall()  
	periods = {}; 
	for ir in xrange( len(row_imtype) ):
	    periods[str(ir+1)] = row_imtype[ir][2]     # unit: sec; Sa unit: cm/s^2
						       # ir+1 is the IM_Type_ID as shown in the table

	# intensity Measure dictionary
	Ek[stanam] = {}

	# PeakAmplitude selection (rupture_ID, some sources donot use part of their rupture))
	query = "select * from %s where Run_ID = %s and Source_ID = %s and Rupture_ID = %s"%('PeakAmplitudes',site_run_id,sid,rid)
	cursor.execute( query )       # run query
	row_Sa = cursor.fetchall()  
	nvar = len(row_Sa)
	if nvar == 0:
	    fid_stdout.write('There is no SA for (Site_name, Site_ID, run_id, sid,rid) = (%s, %s,%s,%s,%s)\n'%(row_sites[2],row_sites[0],site_run_id,sid,rid))
	    return sites_info, Ek
	else:
	    pass
	
	# get all available periods
	Nt = 0; Ts = []
	for ik in xrange( len(row_Sa) ):
	    IMtype = row_Sa[ik][4]
	    if row_Sa[ik][2] == 0:  # just deal with the first rupture variation
		if 1<= IMtype <= 26 or 82<= IMtype <=99:   
		    # For geometric mean and 1-26: 2 to 10 sec, 82-99: 0.1-1.666 sec (broadband)
		    Nt = Nt+1
		    Ts.append( periods[str(IMtype)] )
	    else:
		break

	tmpSa = np.zeros( (Nh, Nf, Nt) )
	for ih in xrange( Nh ):
	    for islip in xrange( Nf ):
		ivar = islip*Nh + ih 
		for it in xrange( Nt ):
		    ik = ivar * Nt + it
		    tmpSa[ih,islip,it] = row_Sa[ik][5]/980.   # in gravity
	Ek[stanam] = tmpSa.tolist()
	
	# all periods information (corresponding to it)
	Ek['periods'] = Ts
  
    return sites_info, Ek 


#Utilities to use OpenSHA produce NGA site flatfiles 
#Could be extended to generate NGA site flatfiles using pynga.utils
def cpt_OpenSHA_nga(cwd, NGAmeta, sids, rids, Ti, SiteName=None, erf_id=35):
    """
    run NGA of OpenSHA
    """
    try: 
	Nt = len(Ti)
    except:
	Ti = [Ti,]
	Nt = len(Ti)

    nga_input = NGAmeta + 'nga_inputs'
    fid = open( nga_input, 'w' )
    for irup in xrange( len(sids) ):
	sid = sids[irup]
	rid = rids[irup]
	for it in xrange( Nt ):
	    fid.write( '%s %s %s\n'%(sid,rid,'%.3f'%Ti[it]) )
    fid.close()
    os.chdir( NGAmeta )
    if SiteName == None:
	os.system( 'nga_comparison_calc nga_inputs' )
    else:
	os.system( 'nga_comparison_calc --site %s nga_inputs'%SiteName )
    os.chdir( cwd ) 

	    
def OpenSHA_nga_files(NGAmeta, sid, rid, Ti, SiteName=None, RunID=None, erf_id=35):
    """
    Return OpenSHA output file names (select using sid, rid, Ti, and SiteName)
    """
    if SiteName == None:
	ftmp = 'ERF%s_src%s_rup%s_SA%s.csv'%(erf_id, sid, rid, '%2.1f'%(Ti))
    else:
	if RunID != None:
	    ftmp = 'ERF%s_%s.csv'%(erf_id,SiteName)
	else:
	    print 'When SiteName is not None, RunID cannot be None'
	    raise ValueError
    
    OpenSHA_output = NGAmeta + ftmp
    if not os.path.exists( OpenSHA_output ):
	print 'OpenSHA computed NGA file %s are not found...'%OpenSHA_output
	raise ValueError
    else:
	if RunID != None and SiteName != None:
	    ftmp0 = 'ERF%s_%s_RunID%s.csv'%(erf_id,SiteName,RunID)
	    shutil.copyfile( NGAmeta+ftmp, NGAmeta+ftmp0 )
            OpenSHA_output = NGAmeta+ftmp0

    return OpenSHA_output


# this will be used by other classes (not used)
# original output from OpenSHA (two ways to use it)
def OpenSHA_nga(NGAmeta,NGAmodel,sid,rid,Ti,SiteName=None,erf_id=35):
    """
    Computed OpenSHA (information rearangement)
    Using dictionary (attention to the structure)
    """
    OpenSHA_output = OpenSHA_nga_files(NGAmeta, sid, rid, Ti, SiteName=SiteName )

    ngaO = {}
    for inga,nga in enumerate( NGAmodel ):
	ngaO[nga] = {}

    if SiteName == None:
	lines = open(OpenSHA_output, 'r').readlines()
	for il in range( 1, len(lines) ):
	    spl = lines[il].strip().split(',')                                                                         
	    stanam = spl[1]   # short name
	    for inga,nga in enumerate( NGAmodel ):
		index1 = 2*inga
		ngaO[nga][stanam] = [float(spl[index1-8]),float(spl[index1-7])]   # return (median,std) in natural log
    else:
	lines = open(OpenSHA_output, 'r').readlines()
	for il in range( 3, len(lines) ):
	    spl = lines[il].strip().split(',')
	    sid = spl[0]
	    rid = spl[1]
	    sr_key = '%s_%s'%(sid,rid)
	    for inga,nga in enumerate( NGAmodel ):
		index1 = 2*inga
		ngaO[nga][sr_key] = [float(spl[index1-8]),float(spl[index1-7])]   # return (median,std) in natural log
    return ngaO

    
def CSmetaRup_gen(metafile, cybershk_database_password):
    """
    Extract rupture info and save as metadata
    """
    rups_info = {}
    
    hostname = 'focal.usc.edu'  # host                            
    username = 'cybershk_ro'    # usrname                         
    password = cybershk_database_password    # password                        
    database = 'CyberShake'     # Database name 
    db = mdb.connect( host = hostname, user = username, \
		 passwd = password, db = database )
    suffix = metafile.strip().split('/')[-1].split('_')
    erf_id,rup_scenario_id = suffix[2:4]
    sid = suffix[4].split('.')[0]
    
    cursor = db.cursor()
    fid_stdout = open( './stdout', 'w' )
    
    rups_info = rup_gen( cursor, sid, rups_info, fid_stdout, rup_scenario_id=rup_scenario_id )
    
    fid_stdout.close()
    cursor.close()
    db.close()

    # save rupture info into file
    meta1 = dict( 
	    rups_info = rups_info, 
	    )
    header = '# Rupture meta file for rupture %s\n'%sid
    save(metafile,meta1,header=header)
    
    return 1


def CSmetafile_gen(metafile, cybershk_sitedata, cybershk_database_password ):

    sites_info = {}
    Ek = {}
    hostname = 'focal.usc.edu'  # host                            
    username = 'cybershk_ro'    # usrname                         
    password = cybershk_database_password    # password                        
    database = 'cybershake'     # database name 
    db = mdb.connect( host = hostname, user = username, \
		 passwd = password, db = database )

    suffix = metafile.strip().split('/')[-1].split('_')
    erf_id,sgt_id, rup_scenario_id, vel_id = suffix[1:5]
    sid,rid = suffix[5], suffix[6].split('.')[0]

    rup_model_ids = (erf_id, sgt_id, rup_scenario_id, vel_id)

    # stdout file
    fid_stdout = open( './stdout', 'w' )
    rup_info = {}
    sites = cybershk_sites( cybershk_sitedata )
    for stanam in sites.keys():
	cursor = db.cursor()
	sites_info, Ek = \
	im_gen( cursor,sid,stanam,\
		sites_info,Ek,fid_stdout, rup_model_ids = rup_model_ids )
	cursor.close()

    meta1 = dict(
		sites_info = sites_info,
		im11 = Ek
		)
    header = '# Meta file for src %s\n'%sid
    save( metafile, meta1, header = header )

    # close what left
    fid_stdout.close()
    db.close()


def SC08_input_NGAinfo(fD_input,NGAmodel,Tlist):
    # NGA_info file
    fid = open(fD_input+'NGA_info','w')
    for inga in xrange( len(NGAmodel) ):
	if inga == len(NGAmodel)-1:
	    fid.write( '%s6\n'%NGAmodel[inga] )
	else:
	    fid.write('%s6, '%NGAmodel[inga] )

    for it in xrange( len(Tlist) ):
	if it != 0:
	    fid.write(', %s'%('%3.2f'%Tlist[it]))
	else:
	    fid.write('%s'%('%3.2f'%Tlist[it]))
    fid.close()
    return 1


def SC08_input_SiteInfo(CSmetafile, RupMetafile, cybershk_sitedata, fD_input, projection, kwds, dimS ):
    if not os.path.exists( CSmetafile ):
	CSmetafile_gen(CSmetafile, cybershk_sitedata)
    meta = load( CSmetafile )
    sites_info = meta.sites_info
    sites_run = meta.sites_info.keys()

    suffix = CSmetafile.strip().split('/')[-1].split('_')
    erf_id,sgt_id,rup_scenario_id, vel_id = suffix[1:5]
    sid = suffix[5].split('.')[0]
    
    meta_rup = load( RupMetafile )
    rups_info = meta_rup.rups_info
    
    # just compute the last one (for geometry parameters)
    Mw = rups_info['MR'][0][-1]
    rid = rups_info['MR'][1][-1]
    
    fD_input1 = fD_input + '/ERF%s_RupVar%s/'%(erf_id, rup_scenario_id)
    if not os.path.exists( fD_input1 ):
	os.mkdir( fD_input1 )
    fD_input2 = fD_input1 + '/SourceID%s_RuptureID%s/'%(sid,rid)
    if not os.path.exists( fD_input2 ):
	os.mkdir( fD_input2 )

    # Site_info file (fix the origion for all cases, or fix the coordinate)
    fid = open(fD_input2+'sites_info','w')
    Nsta = len(sites_run)
    fid.write('%s %s %s %s\n'%(Nsta,dimS[0],dimS[1],dimS[2]))    # header:total # of run sites
    #print kwds
    for ista in xrange( Nsta ):
	slon = sites_info[sites_run[ista]]['lon']
	slat = sites_info[sites_run[ista]]['lat']
	siteID = sites_info[sites_run[ista]]['id']
	
	# Site ID   X (east)   Y(north)  Z(up, here=0)  [loop over sites_run, get siteID from sites_info]
	# X,Y are computed based on the pre-defined project
	Xsta,Ysta = projection( slon, slat,**kwds )
	Xsta,Ysta = Xsta/1000., Ysta/1000.   # should be in km
	fid.write( '%s %s %s %s\n'%(siteID,Xsta,Ysta,0) )
    fid.close()
    return 1


def SC08_input_FaultInfo(RupMetafile,fD_input,projection,kwds):
    # Some bugs related to the strike and hypocenter loation 
    # Strikes cannot change very rapidly (<60) and hypocenter should be 
    # within the fault surface

    if not os.path.exists( RupMetafile):
	CSmetaRup_gen( RupMetafile )

    suffix = RupMetafile.strip().split('/')[-1].split('_')
    erf_id,rup_scenario_id = suffix[2:4]
    sid = suffix[4].split('.')[0]
    
    meta_rup = load( RupMetafile )
    rups_info = meta_rup.rups_info
    
    # just compute the last one (for geometry parameters)
    Mw = rups_info['MR'][0][-1]
    rid = rups_info['MR'][1][-1]

    Nrow = rups_info['MR'][2]
    Ncol = rups_info['MR'][3]
    rrake = rups_info['MR'][3]
    Nh,Nf,hypo_loc = rups_info['hypo_slip']

    # ===================
    # Fault surface info
    # ===================
    # Option2: based on differential strike obtaind from Points and select connection points to get fault segments dstrike(colpoints)
    lon1d = rups_info['fault'][0]
    lat1d = rups_info['fault'][1]
    dep1d = rups_info['fault'][2]

    lon2d = np.array(lon1d).reshape( (Nrow,Ncol) )
    lat2d = np.array(lat1d).reshape( (Nrow,Ncol) )
    dep2d = np.array(dep1d).reshape( (Nrow,Ncol) )
    
    # surface trace of the fault
    lon = rups_info['surface'][0]
    lat = rups_info['surface'][1]
    strike1 = rups_info['surface'][2]
	    
    # Option1: based on Ncol and Nh to devide the entire fault (slower)
    if 0:
	c_point = (np.mod(np.arange( Ncol ), Ncol/Nh)==1).nonzero()[0]  # connect point
	c_point = list( c_point )

	if c_point[-1] != Ncol-1:
	    c_point.append( Ncol-1 )
    else:
	# Option2: use the first and the last point to Construct one segment (faster)
	# doesn't make much difference
	c_point = [0, Ncol-1]

    # set up segment points (lon,lat)
    Nseg = len(c_point)-1    # number of Segments based on connect points
    segs = {}
    segs0 = np.zeros( (4,Nseg,3) )
    for iseg in xrange( Nseg ):
	start_index = c_point[iseg]
	end_index = c_point[iseg+1]
	# four points determine one segment
	Irow = [0,0,Nrow-1,Nrow-1]
	Icol = [start_index,end_index,end_index,start_index]
	for ip in xrange( len(Irow) ):
	    lon_ip = lon2d[Irow[ip],Icol[ip]]
	    lat_ip = lat2d[Irow[ip],Icol[ip]]
	    dep_ip = -dep2d[Irow[ip],Icol[ip]]
	    Xtmp,Ytmp = projection(lon_ip,lat_ip,**kwds)
	    Xtmp,Ytmp = Xtmp/1000., Ytmp/1000.
	    segs0[ip,iseg,0] = Xtmp
	    segs0[ip,iseg,1] = Ytmp
	    segs0[ip,iseg,2] = dep_ip

    # correct the self-strike (distortion due to the UTM projection)
    # the third points of one segment
    for iseg in xrange( Nseg ):
	for icmp in xrange(2):
	    if iseg != 0:
		segs0[3,iseg,icmp] = segs0[2,iseg-1,icmp]
	points = []
	for ip in xrange( 3 ):
	    if ip == 2:
		points.append([segs0[ip+1,iseg,0],segs0[ip+1,iseg,1],segs0[ip+1,iseg,2]])
	    else:
		points.append([segs0[ip,iseg,0],segs0[ip,iseg,1],segs0[ip,iseg,2]])
	dist, point1 = point_plane(segs0[2,iseg,:],points)
	segs0[2,iseg,:] = point1
    
    for iseg in xrange( Nseg ):
	for ip in xrange( 4 ):
	    key = '(%s,%s)'%(iseg,ip)
	    segs[key] = '%s %s %s'%(segs0[ip,iseg,0],segs0[ip,iseg,1],segs0[ip,iseg,2])

    # generate fault segment file
    fD_input1 = fD_input + '/ERF%s_RupVar%s/'%(erf_id, rup_scenario_id)
    if not os.path.exists( fD_input1 ):
	os.mkdir( fD_input1 )
    fD_input2 = fD_input1 + '/SourceID%s_RuptureID%s/'%(sid,rid)
    if not os.path.exists( fD_input2 ):
	os.mkdir( fD_input2 )
    fid_seg = open( fD_input2+'rupture_info', 'w' )
    lines = {}
    for ip in xrange( 4 ):
	key_p = '%s'%ip
	key0 = '(%s,%s)'%(0,ip)
	lines[key_p] = segs[key0]
	for iseg in xrange( Nseg ):
	    if iseg != 0:
		key = '(%s,%s)'%(iseg,ip)
		lines[key_p] = lines[key_p]+' '+segs[key]
	fid_seg.write('%s\n'%lines[key_p])
    fid_seg.close()
    
    # ===================
    # Hypocenter info
    # Hypocenter should be within the segment plane obtained above
    # ===================
    rlon=[]; rlat=[]; rdep = []
    for ih in xrange( Nh ):
	hypo_key = 'h%4.4i'%ih
	rlon.append( hypo_loc[hypo_key][0] )
	rlat.append( hypo_loc[hypo_key][1] )
	rdep.append( hypo_loc[hypo_key][2] )
    rlon = np.array(rlon)
    rlat = np.array(rlat)
    rdep = np.array(rdep)

    for ih in xrange( Nh ):
	# generate hypo info file
	fid_hypo = open( fD_input2+'hypo_%s'%(ih+1),'w' )
	
	# original location of the hypocenters
	srcX,srcY = projection( rlon[ih], rlat[ih], **kwds )    
	srcZ = -rdep[ih]   # in km (Z-up as positive direction)
	srcX, srcY = srcX/1000., srcY/1000.
        
	if 0:
	    # compare with Matlab
	    # this section is done in Matlab isodirect_util/msegcl5.m
	    # search which fault segment the hypocenter is located and do the correction due to projection
	    for iseg in xrange( Nseg ):
		start_index = c_point[iseg]
		end_index = c_point[iseg+1]
		
		lon_se = lon2d[0,start_index],lon2d[0,end_index],lon2d[-1,end_index],lon2d[-1,start_index]
		lat_se = lat2d[0,start_index],lat2d[0,end_index],lat2d[-1,end_index],lat2d[-1,start_index]
		
		# find the segment where hypocenter is located
		if min(lon_se)<= rlon[ih] <= max(lon_se) and min(lat_se) <= rlat[ih] <= max(lat_se):
		    points = []   # three points of the segment (will be used to determine this segment)
		    for ip in xrange( 3 ):
			if ip == 2:
			    points.append([segs0[ip+1,iseg,0],segs0[ip+1,iseg,1],segs0[ip+1,iseg,2]])
			else:
			    points.append([segs0[ip,iseg,0],segs0[ip,iseg,1],segs0[ip,iseg,2]])
		    flag = 1
		    print 'Hypocenter %s is located on Segment %s'%(ih+1,iseg+1)
		    break
		else:
		    flag = 0
		    continue
	    
	    if flag == 0:
		print 'The hypocenter %s cannot be located on fault plane, excluding the rupture'%ih
		raise ValueError

	    # compute the distance and projected point in the plane
	    dist, point1 = point_plane([srcX,srcY,srcZ],points)
	    print dist
	    
	    # use the closest point in the fault plane (projection on the fault plane)
	    #srcX,srcY,srcZ =  point1   
	    #dist, point11 = point_plane( point1, points )
	    #print point1
	    #print dist   # this should be very small

	fid_hypo.write('%s\n'%(srcX))
	fid_hypo.write('%s\n'%(srcY))
	fid_hypo.write('%s\n'%(srcZ))
	fid_hypo.write('%s\n'%rrake)
	fid_hypo.write('%s\n'%Mw)
	fid_hypo.close()

    return Nh


def SC08_input(sidstr,ridstr,Nhstr,fDcpt_pth, util, erf_id=35, rup_scenario_id=3,mflag=1):
    # generate input files for matlab package (all in east-nroth-up coords)
    # inputfile for run.m in scripts/fD_compute/
    fid = open(fDcpt_pth+'/input_file','w')
    fid.write('%s/\n'%util)    # directory where isodirect_util and fD_compute are located
    fid.write('%s\n'%erf_id)   # for ERF id
    fid.write('%s\n'%rup_scenario_id)   # for RupVar
    fid.write('%s\n'%sidstr)   # in 86,89,115 
    fid.write('%s\n'%ridstr)   # in 2,4,2
    fid.write('%s\n'%Nhstr)      # in 17,28,11
    
    if mflag == '0':
	fid.write('%s,%s,%s,%s,%s\n'%(0,0,0,'%.1f'%3.0,'BA6') )   # visually test
    else:
	fid.write('%s,%s,%s,%s,%s\n'%(0,0,0,'%.1f'%3.0, 'BA6') )  # don't plot

    fid.close()
    return 1


def SC08_matlab_run(fDcpt_pth,cwd):
    # run matlab from python
    os.chdir( fDcpt_pth )
    os.system( 'matlab -nojvm < %s/run.m'%fDcpt_pth )
    os.chdir( cwd )
    return 1


# prepare the input file for SC computing of directivity
def cpt_SC08(sids,rids,CSmeta, RupMeta, cybershk_sitedata,\
	     cwd, fDcpt_pth,util, projkwds, dimS, \
	     Tlist = [0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0],\
	     NGAmodel = ['CB','BA','CY','AS'], \
	     rup_model_ids=(35,5,3,1), mflag=1):
    """
    Compute all hypocenters
    sids: [128,]
    rids: [1263,]
    cwd: current work directory
    """
    fD_input = fDcpt_pth+'/inputs/'
    fD_output = fDcpt_pth+'/outputs'
    
    # Set up projection (same for all ruptures)
    # Attention, you should use the same projection for all application (the way to get self.slon1d for example)
    kwds = projkwds
    
    erf_id, sgt_id, rup_scenario_id, vel_id = rup_model_ids

    sidstr = str(sids[0])
    ridstr = str(rids[0])
    Nhs = np.zeros(len(sids),int)
    for irup in xrange( len(sids) ):
	sid = sids[irup]
	rid = rids[irup]
        if irup != 0:
	    sidstr = sidstr + ',' + str(sids[irup])
	    ridstr = ridstr + ',' + str(rids[irup])

	# NGA_info file
	SC08_input_NGAinfo(fD_input,NGAmodel,Tlist)

	# Site_info file (compute all sites flatfile, not just for given source and rupture)
	CSmetafile = CSmeta + '%s/'%sid + \
	                'meta_%s_%s_%s_%s_%s_%s.py'%(erf_id, sgt_id, rup_scenario_id, vel_id,sid,rid)
	# Rupture and Hypocenters
	RupMetafile = RupMeta + '%s/'%sid + \
		   'meta_rup_%s_%s_%s.py'%(erf_id,rup_scenario_id,sid)
	
	SC08_input_SiteInfo( CSmetafile, RupMetafile, cybershk_sitedata, fD_input, projection, kwds, dimS)
	Nhs[irup] = SC08_input_FaultInfo( RupMetafile, fD_input, projection, kwds )
    
    Nhstr = str( Nhs[0] )
    for irup in xrange( len(Nhs) ):
	if irup != 0:
	    Nhstr = Nhstr + ',' + str( Nhs[irup] )
    
    SC08_input(sidstr, ridstr, Nhstr, fDcpt_pth, util, erf_id=erf_id, rup_scenario_id=rup_scenario_id, mflag=mflag)
    SC08_matlab_run( fDcpt_pth, cwd )
    return 1


def SC08_files( sid, rid, fD_output, model_name, erf_id=35, rup_scenario_id = 3):
    
    # Select files using sid,rid
    ftmp = 'ERF%s_RupVar%s/SourceID%s_RuptureID%s/hypo*_%s6.txt'%(erf_id,rup_scenario_id,sid,rid,model_name)
    fnamtmp = glob.glob( fD_output + ftmp )
   
    if len( fnamtmp ) == 0:
	print 'You should first generate directivity files'
	raise ValueError

    # get correct file order (orderred by hypocenter)
    SC08_output = {}
    for ifile in xrange( len( fnamtmp ) ):
	name = fnamtmp[ifile].strip().split('/')[-1].strip().split('_')[0]   # hypo*
	SC08_output[str(int(name[4:])-1)] = fnamtmp[ifile]

    return SC08_output



# =================
# ABF block
# =================
def im_decom1(Gkxmfs,pdf_km, pdf_kx, pdf_kf=None, RupVar=True):
    
    # Gkxmfs could have dim as (Nk,Nr,Nh,Nsta) or (Nk,Nr,Nh,Nf,Nsta)
    # for the case (Nk,Nr,Nh,Nsta), have to set pdf_kf = None (and no return for Nf-related term)

    Nk = len( Gkxmfs )

    if RupVar:
	Fkxmfs = []
    
    Ekxms = []
    Hkxs = []
    
    for ik in xrange( Nk ):

	f_xmfs = np.array( Gkxmfs[ik] )
	
	pdf_m = pdf_km[ik]
	pdf_x = pdf_kx[ik]
	if pdf_kf != None: 
	    pdf_f = pdf_kf[ik]
	    Nm,Nh,Nf,Nsta = f_xmfs.shape
	    f_xms = np.average( f_xmfs, axis=2, weights = pdf_f )   # average over slip variations
	else: 
	    pdf_f = None 
	    Nm,Nh,Nsta = f_xmfs.shape
	    f_xms = f_xmfs 
	
	# average over magnitudes
	f_xs = np.average( f_xms, axis=0, weights=pdf_m )
	Hkxs.append( f_xs )      
        
	if pdf_kf != None:
	    tmp_xmfs = np.zeros((Nm,Nh,Nf,Nsta))
	    for islip in xrange( Nf ):
		tmp_xmfs[:,:,islip,:] = f_xmfs[:,:,islip,:] - f_xms
	    
	    if RupVar:
		Fkxmfs.append( tmp_xmfs )
		#tmp_sd_xms = np.average( tmp_xmfs**2, axis=2, weights=pdf_f ) 
		#SFkxms.append( np.sqrt( tmp_sd_xms ) )
		
		# average the standard deviation over hypocenter and magnitude (weakly dependence)
		#tmp_sd = np.average( np.average( tmp_sd_xms, axis=0, weights=pdf_m), axis = 0, weights = pdf_x )
		#SFks.append( tmp_sd ) 
	
	tmp_xms = np.zeros((Nm,Nh,Nsta))
	for ir in xrange( Nm ):
	    tmp_xms[ir,:,:] = f_xms[ir,:,:] - f_xs
	Ekxms.append( tmp_xms )
	#tmp_sd_xs = np.average( tmp_xms**2, axis=0, weights=pdf_m ) 
	#SEkxs.append( np.sqrt( tmp_sd_xs ) )
	#tmp_sd = np.average( tmp_sd_xs, axis=0, weights=pdf_x )
	#SEks.append( tmp_sd )
	
    if pdf_kf != None and RupVar:
	return Fkxmfs, Ekxms, Hkxs 
    else: 
	return Ekxms, Hkxs 



def im_decom2(Hkxs, Nloc, pdf_kx, pdf_k, pdf_s):
    # This requires each source has the same site sets (pre-selected) for all models 

    Nk = len(Hkxs)
    Eks = np.zeros( (Nk,Nloc) )
    Es = np.zeros( Nloc )
    for ik in xrange( Nk ):
	pdf_x = pdf_kx[ik]
	Nh, Nsta = Hkxs[ik].shape
	tmp = np.average( Hkxs[ik], axis=0, weights=pdf_x )
	Eks[ik] = tmp[:Nloc]
    Es = np.average( Eks, axis=0, weights=pdf_k )
    A = np.average( Es, axis=0, weights=pdf_s )

    # compute dk_xs and its sigma
    Dkxs = []; SDks = []
    for ik in xrange( Nk ):
	Nh,Nsta = Hkxs[ik].shape
	Dkxs.append(np.zeros((Nh,Nloc)))
	for ih in xrange( Nh ):
	    Dkxs[ik][ih,:] = Hkxs[ik][ih,:Nloc] - Eks[ik]
        #SDks.append( np.sqrt( np.average( Dkxs[ik]**2, axis=0, weights=pdf_kx[ik] ) ) ) 

    # compute ck_s and its sigma
    Cks = np.zeros( (Nk,Nloc) )
    for ik in xrange( Nk ):
	Cks[ik,:] = Eks[ik,:] - Es 
    #SCs = np.sqrt( np.average( Cks**2, axis=0, weights = pdf_k ) ) 

    # compute b_s and its sigma 
    Bs = Es - A
    
    return Dkxs, Cks, Bs, A 



# ================================================================
# Scatter plot residuals
# use un-interpolated factors
# ================================================================
# Plot Vs30 CyberShake-Ref (without site effect) and BA-Ref  
# put in the thesis
def CompareVs30(rup_model_ids, mflag, sigma, T):
    
    erf_id, sgt_id, rup_id, vel_id = rup_model_ids
    Bpth = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs0/Bs/Sigma%s'%\
	    (mflag, erf_id, sgt_id, rup_id, vel_id, sigma)   # as test 
    Ppth = './plots/SiteSpecific'
    if not os.path.exists( Ppth ): 
	os.mkdir(Ppth) 

    pfmt = 'eps' 
    Bfile = Bpth + '/CyberShake.NGAs.%3.2f.bs'%T 

    Bs = np.loadtxt(Bfile) 
    Vs30 = Bs[:,2] 
    CS_B = Bs[:,9]   # CS
    CS_Ref_B = Bs[:,10]     #CS-Ref
   
    CS_BA_B = Bs[:,6]    # CS-BA
    
    BA_B = CS_B - CS_BA_B 
    Ref_B = CS_B - CS_Ref_B 
    BA_Ref_B = BA_B - Ref_B 

    # linear fit 
    if 1:
	slope, intercept, r, p, std_err = stats.linregress(Vs30,CS_Ref_B) 
	CSy = slope * Vs30 + intercept  
	
	slope, intercept, r, p, std_err = stats.linregress(Vs30,BA_Ref_B) 
	BAy = slope * Vs30 + intercept  

    # plot CS_Ref_B and BA_Ref_B .vs. Vs30  (group) 
    fig = plt.figure(1, (14,8)) 
    plotnam = '/CyberShake.BA.SiteEffect.Vs30.scatter.'
    
    ax = fig.add_subplot(111) 
    ax.scatter( Vs30, CS_Ref_B, s=80, facecolors='none', edgecolors='r', label='CyberShake' )
    ax.scatter( Vs30, BA_Ref_B, s=80, facecolors='none', edgecolors='b', label='BA08' )
    lg = ax.legend(loc=(0.6,0.65), fontsize=24,scatterpoints=1)
    #lg.draw_frame(False) 
    ax.plot( Vs30, CSy, 'r', lw = 2 )
    ax.plot( Vs30, BAy, 'b', lw = 2 )
    ax.set_ylim([-1.5,1.75])
    ax.set_ylabel( '$b(\it{r})$ (ln)' )
    ax.set_xlabel( '$V_{s30}$ (m/s)' )
    fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300, transparent=True )
    #plt.show()


# Plot Basin depth CB-BA, CY-BA, AS-BA (with basin depth lines), and CyberShake-BA (without basin depth line)
def CompareBasin(rup_model_ids, mflag, sigma, T):
    
    erf_id, sgt_id, rup_id, vel_id = rup_model_ids
    Bpth = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs0/Bs/Sigma%s'%\
	    (mflag, erf_id, sgt_id, rup_id, vel_id, sigma)   # as test 
    Ppth = './plots/SiteSpecific'
    pfmt = 'eps' 
    plt.rc('font',family='Arial')
    Bfile = Bpth + '/CyberShake.NGAs.%3.2f.bs'%T 

    Bs = np.loadtxt(Bfile) 
    Vs30 = Bs[:,2] 
    Z10 = Bs[:,3] 
    Z25 = Bs[:,4]
    CS_B = Bs[:,9:10] # CS 
    CS_NGA = Bs[:,5:9]   # CS-NGA 
    
    NGA_B = -(CS_NGA - CS_B)
    NGA_B = NGA_B - NGA_B[:,1:2]    # NGA - BA
    NGA_B[:,1] = CS_NGA[:,1]        # CyberShake - BA 

    # compute correlations 
    print [ np.corrcoef(NGA_B[:,0],Z25)[0,1],
	    np.corrcoef(NGA_B[:,2],Z10)[0,1], 
	    np.corrcoef(NGA_B[:,3],Z10)[0,1], 
	    np.corrcoef(CS_B[:,0],Z10)[0,1], 
	    np.corrcoef(CS_B[:,0],Z25)[0,1], ] 

    # calculate the basin term (normalized)
    # CB 
    NGA1 = CB08.CB08_nga()
    kwds = {'Tother':T, 'Z25':Z25}
    CB_B = mapfunc(NGA1.basin_function,**kwds) 
    CB_B = CB_B - sum(CB_B)/len(CB_B) 
    
    # CY
    NGA1 = CY08.CY08_nga()
    kwds = {'Tother':T, 'Z10':Z10}
    CY_B = mapfunc(NGA1.basin_function,**kwds) 
    CY_B = CY_B - sum(CY_B)/len(CY_B) 

    # AS
    NGA1 = AS08.AS08_nga()
    kwds = {'Tother':T, 'Z10':Z10,'Vs30':Vs30}   # Vs30 will affect AS basin effects
    AS_B = mapfunc(NGA1.soil_function,**kwds) 
    AS_B = AS_B - sum(AS_B)/len(AS_B) 
    
    # plot CS_Ref_B and BA_Ref_B .vs. Vs30  (group) 
    fig = plt.figure(1, (12,8)) 
    fig.clf() 
    clr1='r'
    clr2='b'
    plt.rc('font',family='Arial')
    fsize=14
    plt.rcParams['xtick.labelsize']=fsize
    plt.rcParams['ytick.labelsize']=fsize
    plt.rcParams['axes.labelsize']=fsize

    plotnam = '/NGAs.BA.BasinsEffect.BasinDepth.scatter.'
    ss = 80    # symbol size
    ax = fig.add_subplot(221)

    ax.scatter( Z25, NGA_B[:,0], s=ss,marker='o',facecolors='none', edgecolors=clr1, label='ABF' )
    ax.scatter( Z25, CB_B, s=ss, marker='^', facecolors='none', edgecolors=clr2, label='MBF' )
    ax.legend(loc=4).draw_frame(False)

    ax.set_ylim([-1.5,1.75])
    ax.set_ylim([-1,1])
    #ax.xaxis.set_ticks_position('top')
    #ax.xaxis.set_label_position('top')
    ax.set_xlabel( '$Z_{2.5}$ (km)' ) 
    ax.set_ylabel( '$b(\it {r})$' )
    ax.text( -0.13, 1.07, '(a)', fontsize=fsize, transform=ax.transAxes )
    ax.grid(True)
    
    ax = fig.add_subplot(222)
    ax.scatter( Z10, NGA_B[:,2], s=ss, marker='o', facecolors='none', edgecolors='r' )
    ax.scatter( Z10, CY_B, s=ss, marker='^',facecolors='none', edgecolors='b' )
    ax.set_ylim([-1.5,1.75])
    ax.set_ylim([-1,1])
    #ax.xaxis.set_ticks_position('top')
    #ax.xaxis.set_label_position('top')
    ax.set_xlabel( '$Z_{1.0}$ (m)' )
    ax.text( -0.13, 1.07, '(b)', fontsize=fsize, transform=ax.transAxes )
    ax.grid(True)
    
    ax = fig.add_subplot(223)
    ax.scatter( Z10, NGA_B[:,3], s=ss, marker='o', facecolors='none', edgecolors='r' )
    ax.scatter( Z10, AS_B, s=ss, marker='^',facecolors='none', edgecolors='b' )
    ax.set_ylim([-1.5,1.75])
    ax.set_ylim([-1,1])
    #ax.xaxis.set_ticks_position('top')
    #ax.xaxis.set_label_position('top')
    ax.set_xlabel( '$Z_{1.0}$ (m)' )
    ax.set_ylabel( '$b(\it {r})$' )
    ax.text( -0.13, 1.07, '(c)', fontsize=fsize, transform=ax.transAxes )
    ax.grid(True)
    
    ax = fig.add_subplot(224)
    ax.scatter( Vs30, NGA_B[:,3], s=ss, marker='o',facecolors='none', edgecolors='r' )
    ax.scatter( Vs30, AS_B, s=ss, marker='^', facecolors='none', edgecolors='b' )
    ax.set_ylim([-1.5,1.75])
    ax.set_ylim([-1,1])
    #ax.xaxis.set_ticks_position('top')
    #ax.xaxis.set_label_position('top')
    ax.set_xlabel( '$V_{S30}$ (m/s)' )
    #ax.set_ylabel( '$b(\it {r})$' )
    ax.text( -0.13, 1.07, '(d)', fontsize=fsize, transform=ax.transAxes )
    ax.grid(True)
    fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300, transparent=True )

    fig = plt.figure(2,(14,6)) 
    plotnam = '/CyberShake.BA.BasinsEffect.BasinDepth.scatter.'
    ax = fig.add_subplot(121)
    ax11 = fig.add_subplot(122)
    ax.scatter( Z25, NGA_B[:,1], s=ss, facecolors='none', edgecolors='r')
    ax.set_xlabel( '$Z_{2.5}$ (km)',color='r' ) 
    ax11.scatter( Z10, NGA_B[:,1], s=ss, facecolors='none', edgecolors='b')
    ax.set_ylim([-1.5,1.75])
    ax11.set_xlabel( '$Z_{1.0}$ (m)', color='b')
    for tks in ax11.get_xticklabels():
	tks.set_color('b')
    ax.set_xlabel( '$Z_{2.5}$ (km)', color='r' ) 
    for tks in ax.get_xticklabels():
	tks.set_color('r')
    ax.set_ylim([-1.5,1.75])

    ax.set_ylabel( '$b(\it {r})$' )
    ax.grid(True)
    fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300, transparent=True )
    
    fig = plt.figure(3,(14,8)) 
    plotnam = '/AS.BA.BasinsEffect.BasinDepth.scatter.'
    ax = fig.add_subplot( 111, projection='3d' ) 
    ax.scatter( Z10, Vs30, NGA_B[:,3], c='k', marker='o' )
    ax.plot( Z10, np.ones(len(Z10))*1000, NGA_B[:,3], 'o',mfc='none',mec='r' )
    ax.plot( 1000*np.ones(len(Z10)), Vs30, NGA_B[:,3], 'o',mfc='none',mec='b' )
    ax.set_xlabel( '$Z_{1.0}$ (m)')
    ax.set_ylabel( '$V_{s30}$ (m/s)') 
    ax.set_zlabel( '$b(\it{r})$' )
    ax.grid(True)
    
    if 0:
	ax = fig.add_subplot(223)
	ax.scatter( Z10, NGA_B[:,3], s=ss, facecolors='none', edgecolors='k' )
	#ax.scatter( Z10, AS_B, s=ss, facecolors='none', edgecolors='b' )
	ax.set_ylim([-1.5,1.75])
	ax.set_xlabel( '$Z_{1.0}$ (m)' )
	ax.set_ylabel( '$b(\it {r})$' )
	ax.text( 0.1, 0.9, '(c)', fontsize=16, transform=ax.transAxes )
	ax.grid(True)
	
	ax = fig.add_subplot(224)
	ax.scatter( Vs30, NGA_B[:,3], s=ss, facecolors='none', edgecolors='k' )
	#ax.scatter( Vs30, AS_B, s=ss, facecolors='none', edgecolors='b' )
	ax.set_ylim([-1.5,1.75])
	ax.set_xlabel( '$V_{S30}$ (m/s)' )
	#ax.set_ylabel( '$b(\it {r})$' )
	ax.text( 0.1, 0.9, '(d)', fontsize=16, transform=ax.transAxes )
	ax.grid(True)
    fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300, transparent=True )
    
    #plt.show()


def BmapCorrelation(rup_model_ids0, rup_model_ids1, ref, mflag, sigma, T):
    
    # opt: Bmap (between bmaps of different models)
    #      BmapFT (between bmap and FTmax)

    erf_id0, sgt_id0, rup_id0, vel_id0 = rup_model_ids0
    erf_id1, sgt_id1, rup_id1, vel_id1 = rup_model_ids1
    
    Bpth0 = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs/Bs/Sigma%s'%\
	    (mflag, erf_id0, sgt_id0, rup_id0, vel_id0, sigma)   # as test 
    Bfile0 = Bpth0 + '/CyberShake.NGAs.%3.2f.bs'%T 
    
    Bpth1 = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs/Bs/Sigma%s'%\
	    (mflag, erf_id1, sgt_id1, rup_id1, vel_id1, sigma)   # as test 
    Bfile1 = Bpth1 + '/CyberShake.NGAs.%3.2f.bs'%T 

    bmap_corr = [] 
    bmapFT_corr = []

    NGAmodels = ['CB','BA','CY','AS'] 
    iref = 0
    for inga in xrange( len(NGAmodels) ): 
	if ref == NGAmodels[inga]: 
	    iref = inga
	    break 

    Bs0 = np.loadtxt(Bfile0) 
    CS_B0 = Bs0[:,9:10] # CS 
    CS_NGA0 = Bs0[:,5:9] # CS-NGA
    
    NGA_B0 = -(CS_NGA0 - CS_B0)     # NGA
    NGA_B0 = NGA_B0 - NGA_B0[:,1:2]    # NGA - BA
    NGA_B0[:,1] = CS_NGA0[:,1]        # CyberShake - BA 
    
    print 'model0 and Ref, model1 and Ref, model0 and model1'
    bmap_corr.append( np.corrcoef(NGA_B0[:,1],NGA_B0[:,iref])[0,1] )   # Between CyberShake(0)-BA and REF-BA
    
    
    Bs1 = np.loadtxt(Bfile1) 
    CS_B1 = Bs1[:,9:10] # CS 
    CS_NGA1 = Bs1[:,5:9]   # CS-NGA 
    
    NGA_B1 = -(CS_NGA1 - CS_B1)
    NGA_B1 = NGA_B1 - NGA_B1[:,1:2]    # NGA - BA
    NGA_B1[:,1] = CS_NGA1[:,1]        # CyberShake - BA 
    
    bmap_corr.append( np.corrcoef(NGA_B1[:,1],NGA_B1[:,iref])[0,1] )   # Between CyberShake(1)-BA and REF-BA
    bmap_corr.append( np.corrcoef(NGA_B0[:,1],NGA_B1[:,1])[0,1] )   # Between CyberHSake(0)-BA and CyberShake(1)-BA 

    # From Wavelet analysis
    if 0: 
	bmap_corr0 = []
	P0s = []; T0s = []
	P1s = []; T1s = []
	for icmp in [1,2]:
	    FTpth0 = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Excitation%s/ERF%s_SGT%s_RupVar%s_Vel%s'%\
		    (mflag, erf_id0, sgt_id0, rup_id0, vel_id0)   
	    FTfile0 = FTpth0 + '/FTmax_average%3.2f_cmp%s.txt'%(T,icmp) 
	    
	    FTpth1 = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Excitation%s/ERF%s_SGT%s_RupVar%s_Vel%s'%\
		    (mflag, erf_id1, sgt_id1, rup_id1, vel_id1)  
	    FTfile1 = FTpth1 + '/FTmax_average%3.2f_cmp%s.txt'%(T,icmp) 
	    FT0 = np.loadtxt( FTfile0 )
	    P0 = FT0[:,2] 
	    T0 = FT0[:,3] 
	    bmap_corr0.append( np.corrcoef( NGA_B0[:,1], P0 )[0,1] )
	    bmap_corr0.append( np.corrcoef( NGA_B0[:,1], T0 )[0,1] )

	    FT1 = np.loadtxt( FTfile1 )
	    P1 = FT1[:,2] 
	    T1 = FT1[:,3] 
	    
	    bmap_corr0.append( np.corrcoef( NGA_B1[:,1], P1 )[0,1] )
	    bmap_corr0.append( np.corrcoef( NGA_B1[:,1], T1 )[0,1] )
	    P0s.append( P0 ) 
	    T0s.append( T0 ) 
	    P1s.append( P1 )
	    T1s.append( T1 )
	
	bmap_corr.append( bmap_corr0 )
    
	# geometrical mean of two components
	bmap_corr.append( np.corrcoef( NGA_B1[:,1], np.sqrt( P0s[0]**2 + P0s[1]**2) )[0,1] )
	bmap_corr.append( np.corrcoef( NGA_B1[:,1], np.sqrt( T0s[0]**2 + T0s[1]**2) )[0,1] )
	bmap_corr.append( np.corrcoef( NGA_B1[:,1], np.sqrt( P1s[0]**2 + P1s[1]**2) )[0,1] )
	bmap_corr.append( np.corrcoef( NGA_B1[:,1], np.sqrt( T1s[0]**2 + T1s[1]**2) )[0,1] )
    

    return bmap_corr   



def CompareDistance(rup_model_ids, mflag, sigma, T): 
    erf_id, sgt_id, rup_id, vel_id = rup_model_ids
    Cpth = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs0/Cks/Sigma%s'%\
	    (mflag, erf_id, sgt_id, rup_id, vel_id, sigma)   # as test 
    Ppth = './plots/Distance'
    if not os.path.exists( Ppth ): 
	os.mkdir( Ppth )
    
    pfmt = 'eps' 
    ss = 5    # symbol size
    lgs = ['CB08','BA08','CY08','AS08','CyberShake'] 
    lgsR = ['Rjb','Rrup','Rx'] 
   
    T = [2.0, 3.0, 5.0, 10.0] 

    tmp1 = []
    for it in xrange( len(T) ): 
	period = T[it]
	Cfiles = glob.glob( Cpth + '/CyberShake.NGAs.%3.2f.Source*.cks'%(period)  )
    
	#pdf_k = {'87': 0.031381083105959572, '10': 0.056710843160550251, '89': 0.019702661875046653, '15': 0.020430098754702634, '93': 0.018678811463310765, '218': 0.03234306739048444, '219': 0.031757038801596513, '271': 0.034166418522314709, '88': 0.018571395550180343, '273': 0.039605198930116908, '64': 0.055258008944139396, '254': 0.054441512035385956, '112': 0.076403987169896714, '267': 0.022970009205114782, '8': 0.019837271690223278, '255': 0.02576486248717012, '86': 0.13627204235727591, '231': 0.050274046547844217, '232': 0.062252280548632764, '85': 0.19317936146005402}
	#sum_test = []
	#Nsrc = 0
	#weights = []
        tmp = []; sids = []
	R = []
	for Cfile in Cfiles:
	    sid = Cfile.strip().split('/')[-1].split('.')[-2][6:]
	    sids.append( sid )  

	    #print 'Source %s'%sid
	    Cks = np.loadtxt(Cfile) 

	    Rjb = Cks[:,2] 
	    Rrup = Cks[:,3] 
	    Rx = Cks[:,4] 
	    #Dist = [Rjb, Rrup, Rx]  # km
	    R.append( Rjb + 0.0001 )
	    Nsta = len(Rx) 

	    CS_C = Cks[:,9:10]     # CS 
	    
	    CS_NGA = Cks[:,5:9]   # CS-NGA 
	    CS_Ref = Cks[:,10:11]   # CS-Ref
	    
	    NGA_C = -(CS_NGA - CS_C)    # NGA
	    Ref_C = -(CS_Ref - CS_C)    # Ref  
	    
	    tmp.append([NGA_C[:,0],NGA_C[:,1],NGA_C[:,2],NGA_C[:,3],CS_C[:,0]])   # CB,BA,CY,AS,CS
            #print np.array( tmp ).shape
	    NGA_C = NGA_C - Ref_C[:,0:1]   # NGA - Ref (referece) 
        tmp1.append( tmp )  # [Nperiod, Nsrc, Nmodel]
    tmp1 = np.array(tmp1)
    Nt, Ns, Nm, Nsta = tmp1.shape 
    
    from scipy import polyfit, polyval 

    # CyberShake
    fig = plt.figure(1, (10,40)) 
    fig.clf() 
    plotnam = '/CyberShake.NGAmean.Crk.'
    p = []
    for i in xrange( Ns ): 
	for j in xrange( Nt ): 
	    isub = j + Nt*i  
	    ax = fig.add_subplot(Ns,Nt,isub+1)
	    
	    # mean of CyberShake
	    mu = sum(tmp1[j,i,-1])/Nsta
	    mu = 0.0
	    CS_D = tmp1[j,i,-1] - mu 
	    NGA_D = (tmp1[j,i,0]+tmp1[j,i,1]+tmp1[j,i,2]+tmp1[j,i,3])/4.-mu
	    if 1:
		p0 = ax.scatter( R[i], CS_D, s=ss, facecolors='none', edgecolors='r' )
		p1 = ax.scatter( R[i], NGA_D, s=ss, facecolors='none', edgecolors='b' )
		ax.plot( R[i], np.zeros(len(R[i])), 'k--' )
		ax.set_ylim([-2.3,2.3])
	    else:
		#ax.set_xscale( 'log' ) 
		npoly = 1 
		w = polyfit( np.log10(R),CS_D, npoly) 
		line = polyval( w, np.log10(R)) 
		p0 = ax.semilogx( R, line, 'r-' ) 
		w = polyfit( np.log10(R),NGA_D, npoly) 
		line = polyval( w, np.log10(R)) 
		p1 = ax.semilogx( R, line, 'b-' ) 
		p0 = p0[0]
		p1 = p1[0]
	    if i == 0 and j == 0: 
		p.append( p0 ) 
		p.append( p1 ) 
	    if i == 0 and j == 1:
		ax.legend( p, tuple( ['CS11','NGAmean'] ), scatterpoints=1, loc=3, bbox_to_anchor=(0, 1.5, 1.0, .102), title='C(r,k)')

	    if i != Ns-1: 
		ax.set_xticks([])
		if i == 0: 
		    ax.set_title( 'SA-%s'%('%3.1fs'%T[j]) )

	    else: 
		ax.set_xticks([-20,20,60,100])
		ax.set_xlabel('Distance (km)')
	    
	    if j != 0: 
		ax.set_yticks([])
	    else: 
		ax.set_yticks([-2.0,0.0,2.0])
		ax.set_ylabel(sids[i]) 


    if pfmt == 'png':
	fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300, transparent=False )
    else: 
	fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300, transparent=True )
    


def CompareDirectivity(rup_model_ids, mflag, sigma, T, sid): 
    erf_id, sgt_id, rup_id, vel_id = rup_model_ids
    Dpth = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs0/Dkxs/Sigma%s'%\
	    (mflag, erf_id, sgt_id, rup_id, vel_id, sigma)   # as test 
    Ppth = './plots/Directivity'
    if not os.path.exists( Ppth ): 
	os.mkdir( Ppth )
    Ppth = Ppth + '/%s'%sid 
    if not os.path.exists( Ppth ): 
	os.mkdir( Ppth )

    pfmt = 'png' 
    ss = 40    # symbol size
    fig = plt.figure(1, (12,6)) 
    fig.clf() 
    Dfiles = glob.glob( Dpth + '/%s/CyberShake.NGAs.%3.2f.Source%s.Ih*.dkxs'%(sid,T,sid)  )
    for Dfile in Dfiles:
	ihypo = Dfile.strip().split('/')[-1].split('.')[-2][2:]
	print 'Hypocenter %s'%ihypo
	Dkxs = np.loadtxt(Dfile) 
    
	# use BA coefficients to compute the directivity parameters
	fD = Dkxs[:,3]   # already normalized
	IDP = Dkxs[:,2] 
	
	CS_D = Dkxs[:,8:9]     # CS 
	CS_NGA = Dkxs[:,4:8]   # CS-NGA 
	CS_Ref = Dkxs[:,9:10]   # CS-Ref
	
	NGA_D = -(CS_NGA - CS_D)    # NGA
	Ref_D = -(CS_Ref - CS_D)    # Ref

	NGA_D = NGA_D - Ref_D[:,0:1]   # NGA - Ref (referece) To get residual directivity effect in NGA models
	print np.corrcoef(NGA_D[:,1],fD)[0,1] 

	# Plot CS_Ref, and NGA_D vs IDP, and plot fD superposed
	plotnam = '/CyberShake.BA08.Ref.IDP.Ih%s.scatter.'%(ihypo)
	fig.clf()
	ax = fig.add_subplot(121)
	ax.scatter( IDP, NGA_D[:,1], s=ss, facecolors='r', edgecolors='r', label='BA08d' )
	ax.scatter( IDP, fD, s=ss, facecolors='none', edgecolors='k', label='$f_D$' )
	ax.set_ylim([-1.5,1.75])
	ax.set_xlim([0,5])
	ax.xaxis.set_ticks_position('top')
	ax.xaxis.set_label_position('top')
	ax.set_xlabel( 'IDP' ) 
	ax.set_ylabel( '$d(\mathbf {r},\it{k},\it{x})$' )
	ax.text( 0.1, 0.9, '(a)', fontsize=16, transform=ax.transAxes )
	ax.grid(True)
	ax.legend(loc=0)
	
	ax = fig.add_subplot(122)
	ax.scatter( IDP, CS_Ref, s=ss, facecolors='r', edgecolors='r', label='CyberShake' )
	ax.scatter( IDP, fD, s=ss, facecolors='none', edgecolors='k', label='$f_D$' )
	ax.set_ylim([-1.5,1.75])
	ax.set_xlim([0,5])
	ax.xaxis.set_ticks_position('top')
	ax.xaxis.set_label_position('top')
	ax.set_xlabel( 'IDP' )
	ax.text( 0.1, 0.9, '(b)', fontsize=16, transform=ax.transAxes )
	ax.grid(True)
	ax.legend(loc=0)
	if pfmt == 'png':
	    fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300, transparent=False )
	else: 
	    fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300, transparent=True )


def directivitySD(sigmas, slabs1, slabs2, target1=(35,5,3,1), target2=None, reference='BA', mflag=54, Ts=[3.00,]):
    """
    Parameter space study for the directivity effects
    1. simple way: 
    # just use SDks and read the pdf_k and pdf_s to do the average to get sigma value (so far use this one) 

    2. complicated way:
    Get Dkxs for different models (specified by target and reference) 
    then read in the hypocenter location information to compute the original sigma_d(k,s) for a given source at a site 
    then average over k and s to get a number 
    """ 
    
    NGAs = ['CB','BA','CY','AS'] 
    for i in xrange( len(NGAs) ): 
	if reference == NGAs[i]: 
	    iref = i 
	    break 
    
    # File Path
    erf_id, sgt_id, rup_id, vel_id = target1
    Gpth1 = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs0'%(mflag, erf_id, sgt_id, rup_id, vel_id)
    SourceInfo = Gpth1 + '/SourceInfo' 
    lines = open(SourceInfo).readlines() 
    sids = []; pdf_k = []
    for il in xrange( len(lines) ): 
	spl = lines[il].strip().split()
	sids.append(spl[0]) 
	pdf_k.append( float(spl[2] ) )

    sigmaD = {} 
    key1 = str(target1)
    key3 = '%s-%s'%(str(target1), reference)
    sigmaD[key1] = {}
    sigmaD[key3] = {}
    for isigma in xrange( len(sigmas) ): 
	sigmaKey = sigmas[isigma] 
	sigma_tmp1 = []; sigma_tmp3 = []
	for it in xrange( len(Ts) ):
	    T = Ts[it]
	    tmp1 = []
	    tmp3 = []
	    for ik in xrange( len(sids) ): 
		sid = sids[ik]
		SdFile1 = Gpth1 + '/SDks/Sigma%s/CyberShake.NGAs.%s.Source%s.Sdks'%(sigmaKey, '%.2f'%T, sid) 
		output1 = np.loadtxt( SdFile1, usecols=(iref+2,4+2) ) 
		tmp1.append( output1[:,1] )
		tmp3.append( output1[:,0] )
	    tmp1 = np.array(tmp1)
	    tmp3 = np.array(tmp3)
	    sigma_tmp1.append( np.average( np.average(tmp1,axis=1,weights=np.ones(tmp1.shape[1])/tmp1.shape[1]), axis=0, weights=pdf_k ) )
	    sigma_tmp3.append( np.average( np.average(tmp3,axis=1,weights=np.ones(tmp1.shape[1])/tmp1.shape[1]), axis=0, weights=pdf_k ) )
	sigmaD[key1][sigmaKey] = sigma_tmp1
	sigmaD[key3][sigmaKey] = sigma_tmp3
    
    if target2 != None:
	erf_id, sgt_id, rup_id, vel_id = target2
	Gpth2 = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs0'%(mflag, erf_id, sgt_id, rup_id, vel_id)
	key2 = str(target2)
	key4 = '%s-%s'%(str(target2), reference)
	sigmaD[key2] = {}
	sigmaD[key4] = {}
	for isigma in xrange( len(sigmas) ): 
	    sigmaKey = sigmas[isigma] 
	    sigma_tmp2 = []; sigma_tmp4 = []
	    for it in xrange( len(Ts) ):
		T = Ts[it]
		tmp2 = []
		tmp4 = []
		for ik in xrange( len(sids) ): 
		    sid = sids[ik]
		    SdFile2 = Gpth2 + '/SDks/Sigma%s/CyberShake.NGAs.%s.Source%s.Sdks'%(sigmaKey, '%.2f'%T, sid) 
		    output2 = np.loadtxt( SdFile2, usecols=(iref+2,4+2) ) 
		    tmp2.append( output2[:,1] )
		    tmp4.append( output2[:,0] )
		tmp2 = np.array(tmp2)
		tmp4 = np.array(tmp4)
		sigma_tmp2.append( np.average( np.average(tmp2,axis=1,weights=np.ones(tmp2.shape[1])/tmp2.shape[1]), axis=0, weights=pdf_k ) ) 
		sigma_tmp4.append( np.average( np.average(tmp4,axis=1,weights=np.ones(tmp4.shape[1])/tmp4.shape[1]), axis=0, weights=pdf_k ) )
	    sigmaD[key2][sigmaKey] = sigma_tmp2
	    sigmaD[key4][sigmaKey] = sigma_tmp4
    
    # Plot 
    if target2 != None:
	Keys = [key1,key2,key3,key4] 
	clrs = ['r','b','k','g']
    else: 
	Keys = [key1, key3] 
	clrs = ['r','k'] 

    pfmt = 'eps' 
    pfmt = 'png'
    plt.rc('font',family='Arial')
    Ppth = './plots/DirectivitySD'
    if not os.path.exists( Ppth ): 
	os.mkdir( Ppth )
    
    # plot CHD vs. sigma for given T (different models)
    xt = np.arange( len(sigmas) ) 
    fig = plt.figure(1, (12,8)) 
    for it in xrange(len(Ts)):
	T = Ts[it]
	fig.clf() 
	ax = fig.add_subplot(111) 
	ikey = 0
	for key in Keys: 
	    values = [] 
	    for i in xrange( len(sigmas) ): 
		sigmaKey = sigmas[i]
		values.append( sigmaD[key][sigmaKey][it] )
	    ax.plot(xt, values,clrs[ikey]+'o',label=key)
	    ikey += 1

	ax.set_xticks(xt)
	xticknames = plt.setp(ax, xticklabels=slabs1)
	plt.setp( xticknames, rotation=30,fontsize=14)
	ax.set_ylabel(r'$\langle\sigma_{d}(r,k)\rangle_{r,k}$')
	ax.set_xlabel('Conditional Hypocenter Distributions')
	ax.set_xlim([-1,7])
	ax.set_ylim([0.0, 0.7])
	ax.set_title( 'SA at T = %s s'%(slabs2[it]) )
	ax.legend(loc=0)
	plotnam='/SingleSigma_T%s_Ref%s.'%('%.2f'%T,reference)
	fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300 )

    # plot T vs sigma for given CHD (different models)
    xt = np.arange( len(Ts) ) 
    fig = plt.figure(1, (12,8)) 
    for i in xrange(len(sigmas)): 
	sigmaKey = sigmas[i] 
	fig.clf() 
	ax = fig.add_subplot(111) 
	ikey = 0
	for key in Keys: 
	    values = sigmaD[key][sigmaKey]
	    print sigmaKey, key, values
	    ax.plot(xt, values, clrs[ikey]+'o',label=key)
	    ikey += 1
	
	ax.set_xticks(xt)
	xticknames = plt.setp(ax, xticklabels=slabs2)
	plt.setp( xticknames, rotation=30,fontsize=14)
	ax.set_ylabel(r'$\langle\sigma_{d}(r,k)\rangle_{r,k}$')
	ax.set_xlabel('period [s]')
	ax.set_xlim([-1.0, 4.0])
	ax.set_ylim([0.0, 0.7])
	ax.set_title( 'CHD: %s'%slabs1[i])
	ax.legend(loc=0)
	plotnam='/SingleSigma_CHD%s_Ref%s.'%(sigmaKey,reference)
	fig.savefig( Ppth + plotnam + pfmt, format=pfmt, dpi=300 )

    

if __name__ == '__main__': 
    import sys
    opt = sys.argv[1] 

    mflag =54
    if opt == 'Regression':
	rup_model_ids = (35,5,3,1)
	sigma = '1.00_1.00'  # uniform
	period = float(sys.argv[2])
	#CompareVs30(rup_model_ids, mflag, sigma, period) 
	CompareBasin(rup_model_ids,mflag,sigma,period)
	#CompareDistance(rup_model_ids, mflag, sigma, period)
	#CompareDirectivity(rup_model_ids, mflag, sigma, period, 255)
    
    if opt == 'BasinCorr': 
	rup_model_ids0 = (35,7,4,1)
	rup_model_ids1 = (35,7,4,4)
	sigma = '1.00_1.00'
	Ts = [2.0, 3.0, 5.0, 10.0]
	for T in Ts:
	    basin_corr = BmapCorrelation(rup_model_ids0, rup_model_ids1,'AS',mflag,sigma,T ) 
	    print T, basin_corr

    if opt == 'CHDs': 
	sigmas = ['1.00_1.00','2.00_2.00','0.50_0.50','2.00_3.50','0.50_1.50','3.50_2.00','1.50_0.50']
	slabs1 = ['UN','CB','PB','A1-S','A2-S','A1-N','A2-N']
	
	sigmas = ['1.00_1.00',]
	slabs1 = ['UN',]
	#Ts = [3.0, 5.0, 10.0]
	#slabs2 = ['3.0','5.0','10.0']
	#Ts = [3.0,] 
	#slabs2 = ['3.0',]
        #directivitySD(sigmas, slabs1, slabs2,  target1=(35,5,3,1), target2=None, reference='BA', mflag=mflag, Ts=Ts)

	Ts = [2.0, 3.0, 5.0, 10.0]
	slabs2 = ['2.0', '3.0','5.0','10.0']
	directivitySD(sigmas, slabs1, slabs2, target1=(35,7,4,1), target2=(35,7,4,4), reference='BA', mflag=mflag, Ts=Ts)
	#directivitySD(sigmas, slabs1, slabs2, target1=(35,5,3,1), target2=(35,7,4,1), reference='BA', mflag=mflag, Ts=Ts)
