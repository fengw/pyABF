#!/usr/bin/env python
"""
Class of doing analysis using CyberShake 
including: 
    ABF (compare different hazard models) averaging-based factorization
"""

from utils import * 

class cybershk_nga:
    """
    Class of CyberShake Project
    Combined with NGA
    """
    def __init__( self, wkd, cybershk_database_Password,  \
                  rup_model_ids=(35,5,3,1), \
                  sids=[79,], periods=[3.0,], mflag='0', bb=None, \

                  NGAs = {'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
                          'BA':{'NewCoefs':None,'terms':(1,1,1)},\
                          'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
                          'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)},\
                          'SC':{'NewCoefs':{'CB':None,'BA':None,'CY':None,'AS':None},\
                                'predictors':{'cuteps':{'vr/vs':0.8,'c_cut':2.45,'s_cut':75,'r_cut':0.2,'d_cut':[40,70],'m_cut':[5.6,6.0]},}} \
                          }, \
                  Reference = { \
                      'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
                      'BA':{'NewCoefs':None,'terms':(1,1,1)},\
                      'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
                      'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)},\
                      'SC':{'NewCoefs':{'CB':None,'BA':None,'CY':None,'AS':None},\
                            'predictors':{'cuteps':{'vr/vs':0.8,'c_cut':2.45,'s_cut':75,'r_cut':0.2,'d_cut':[200,250],'m_cut':[5.6,6.0]},}} \
                      },  \
                  ngaModelVersion='2014'
		  ):

	# Reference models (two flags, 4 digits)
	self.rup_model_ids = rup_model_ids
	self.erf_id, self.sgt_id, self.rup_scenario_id, self.vel_id = rup_model_ids
	self.bb = bb   # broadband

	self.prefix_id = 'ERF%s_SGT%s_RupVar%s_Vel%s'%(self.erf_id, self.sgt_id, self.rup_scenario_id, self.vel_id)
	self.prefix_id1 = 'ERF%s_RupVar%s'%(self.erf_id, self.rup_scenario_id)

	self.sids = sids
	self.periods = periods

	# set outside this class (make notes for all possible models)
	self.mflag = mflag      # Model flag (ABF)  (related with NGAs and Reference dictionary)

	# dictionary with key: CB,BA,CY,AS, and with subdictionary contains updated (if no, then set None) coefficients
	self.NGAs = NGAs     
	self.Reference = Reference
	
	self.ngaModelVersion = ngaModelVersion 

	# Those paths should be ready before running anything
	self.wkd = wkd
	self.util = os.path.join( self.wkd, 'utils' )
	self.scripts = os.path.join( self.wkd, 'scripts' )     # changable
	self.datapth = os.path.join( self.wkd, 'data' )    # data 

	# ===================== The following folders will be generated 
	self.plots = os.path.join( self.wkd, 'plots' )     
	#self.isodirect_pth = os.path.join( self.util, 'isodirect_util')
	#self.rNGA_pth = os.path.join( self.util, 'nga_util' )
	
	# Directivity flatfile utility (compute directivity)
	self.fDcpt_pth = os.path.join( self.util, 'fD_compute' )
	self.fD_input = self.fDcpt_pth+'/inputs/' 
	self.fD_output = self.fDcpt_pth + '/outputs/' 

	self.metapth = os.path.join( self.wkd, 'metadata/' )   # preprocess
	self.CSmeta0 = self.metapth+'CyberShake/'   # CyberShake
	self.RupsMeta0 = self.metapth + 'Ruptures/' # Ruptures
	self.CSmeta = self.CSmeta0 + self.prefix_id + '/'
	self.RupsMeta = self.RupsMeta0 + self.prefix_id1 + '/'
	self.NGAmeta = self.metapth + 'NGA08/'     # OpenSHA NGA outputs

	# major outputs for analysis (revise the directory distribution)
	self.mapin = os.path.join( self.wkd, 'scripts/map_input/' )
	self.modelrup_pth = self.mapin + 'Model_Rups%s/'%(self.mflag)
	self.output_pth = os.path.join( self.modelrup_pth, self.prefix_id )

	# GMT plot related
	self.mapin_sites = self.mapin + 'sites_map'
	self.mapin_rups = self.mapin + 'rups_map'

	self.modelrup_plot_pth = self.plots + '/Model_plots/'
	self.modelrup_plot_pth0 = self.modelrup_plot_pth + 'Model%s/'%self.mflag
	self.modelrup_plot_pth00 = self.modelrup_plot_pth0 + self.prefix_id 

	# scatter plot related
	self.flatfile_plot_pth = self.plots + '/FlatFile_plots'
	self.flatfile_plot_pth0 = self.flatfile_plot_pth + '/Model%s/'%self.mflag
	self.flatfile_plot_pth00 = self.flatfile_plot_pth0 + self.prefix_id

	# create required paths
	newpth = [ self.plots, \
	           self.metapth, self.mapin, \
	           self.modelrup_pth, self.output_pth, \
	           self.fD_input, self.fD_output, self.NGAmeta, \
	           self.CSmeta0, self.RupsMeta0, self.CSmeta, self.RupsMeta,  \
	           self.mapin_sites, self.mapin_rups,\
	           self.modelrup_plot_pth, self.modelrup_plot_pth0, self.modelrup_plot_pth00, \
	           self.flatfile_plot_pth, self.flatfile_plot_pth0, self.flatfile_plot_pth00, \
	           ]

	for ipth in xrange( len(newpth) ):
	    f = newpth[ipth]
	    if not os.path.exists( f ):
		os.mkdir( f )

	# NGA models (08 and SC08 available periods)
	self.NGAmodel = ['CB','BA','CY','AS',]
	self.Tlist = [2.0, 3.0, 5.0, 10.0]

	# CyberShake database
	hostname = 'focal.usc.edu'  # host                            
	username = 'cybershk_ro'    # usrname                         
	password = cybershk_database_Password   # password                        
	database = 'CyberShake'     # Database name 
	self.cybershk_database_psword = cybershk_database_Password
	
	try:
	    self.db = mdb.connect( host = hostname, user = username, \
	                           passwd = password, db = database )
	except:
	    self.db = None

	# CyberShake geographical domain (same for all application in the workflow!!!)
	zone = 11
	rot = -30  # degree
	origin = -119.337, 34.15
	inverse = False    # from lon/lat to xy
	self.kwds = {'zone':zone,'origin':origin,'rot':rot,'inverse':inverse}   # all the same projection (global)

	# CyberShake box (for interpolation)
	Lx, Ly = 199000, 113000    # in meters (to avoid the extrapolation error))
	delta = 5000   # in meters  grid spacing for the regular mesh
	shape = (
	    int( Lx/delta + 1.0 ),
	    int( Ly/delta + 1.0 ),
	)

	# Cartesian grid (for interpolation)
	nx = int(Lx/delta) + 1
	ny = int(Ly/delta) + 1
	x = np.arange( nx ) * delta
	y = np.arange( ny ) * delta

	xx,yy = np.meshgrid( x, y )

	kwds = {'zone':zone,'origin':origin,'rot':rot,'inverse':inverse}   # all the same projection (global)
	kwds['inverse'] = True
	slon2d, slat2d = projection(xx,yy,**kwds)
	Ny = xx.shape[0]
	Nx = xx.shape[1]

	# used in interpolation (ABF)
	self.Nloc = Ny*Nx

	# used in fD_compute
	self.dimS = Nx, Ny, delta/1000.

	self.slon1d = slon2d.reshape( self.Nloc )                                                 
	self.slat1d = slat2d.reshape( self.Nloc )

	# interpolation options (ABF)
	self.eps = 0.001
	smooth = 0.01
	self.method = {'name':'exp', 'smooth':smooth}

	# sitedata (used in cybershk_sites function)
	if self.vel_id in [2,4,]: 
	    # CVM-H Basin depth info for the same set of sites
	    self.cybershk_sitedata = wkd+'/data/Sites/cybershk_sites_info_cvmh.dat'
	if self.vel_id in [1,]: 
	    # CVM-S basin depth info (if for updated CVM-S, the basin depth value might change, keep this in mind when you show your results!!!)
	    self.cybershk_sitedata = wkd+'/data/Sites/cybershk_sites_info_cvms4.dat'

	if self.vel_id in [5,]: 
	    # CVM-S4i26 basin depth info (if for updated CVM-S, the basin depth value might change, keep this in mind when you show your results!!!)
	    self.cybershk_sitedata = wkd+'/data/Sites/cybershk_sites_info_cvms4i26.dat'

	if self.vel_id in [8,]: 
	    # BBP 1D velocity model, basin depth info 
	    self.cybershk_sitedata = wkd+'/data/Sites/cybershk_sites_info_bbp1d.dat'
	    print 'Site file:', self.cybershk_sitedata 

	self.sites = cybershk_sites( self.cybershk_sitedata )


    # Methods for the class
    def extract_database(self, sid, rid, RupVar=False, rup_model_ids=None):
	"""
	Source and Rupture based extraction (given one pair of source and rupture)
	Extract all verified CyberShake Sites
	"""
	if rup_model_ids == None: 
	    rup_model_ids = self.rup_model_ids 
	else: 
	    if len(rup_model_ids)<4: 
		print 'rup_model_ids could be None or (erf_id,sgt_id,rvid,vel_id)!'
		raise ValueError

	erf_id, sgt_id, rup_scenario_id, vel_id  = rup_model_ids 

	# original CyberShake IM 
	CSmeta = self.CSmeta + '%s/'%sid 
	if not os.path.exists( CSmeta ): 
	    os.mkdir(CSmeta)
	metafile = CSmeta + \
	    'meta_%s_%s_%s_%s_%s_%s.py'%(erf_id,sgt_id,rup_scenario_id,vel_id,sid,rid)

	# Original Sources 
	RupsMeta = self.RupsMeta + '%s/'%sid
	if not os.path.exists( RupsMeta ): 
	    os.mkdir(RupsMeta)
	metafile1 = RupsMeta + \
	    'meta_rup_%s_%s_%s.py'%(erf_id,rup_scenario_id,sid)

	# stdout file (for debug)
	fid_stdout = open( './stdout', 'w' )

	# Rupture Metafile (both extraction and saveing and loading metafiles are not time consuming)
	if not os.path.exists( metafile1 ):

	    # extract the rupture (sid,rid)
	    rups_info = {}
	    cursor = self.db.cursor()
	    rups_info = rup_gen( cursor, sid, rups_info, fid_stdout, erf_id=erf_id, rup_scenario_id=rup_scenario_id )
	    meta1 = dict( 
	        rups_info = rups_info, 
	    )
	    header = '# Rupture meta file for Source %s\n'%sid
	    save(metafile1,meta1,header=header)

	    # turn meta1 dictionary key to attributes to use below without load meta file
	    meta1 = namespace(meta1)

	else:
	    meta1 = load( metafile1 )

	# IM metafiles (Extraction takes time, and loading time is much longer than the extracting)
	if not os.path.exists( metafile ):

	    sites_info = {}
	    Ek = {}

	    fid_stdout.write('Meta.py generation = database extraction\n')
	    fid_stdout.write( 'Source_ID=%s\n'%(sid) )

	    # all sites
	    sites = self.sites     # cybershk_sites( self.cybershk_sitedata )

	    for stanam in sites.keys():
		# load one by one, for each one, it's fast, so the 
		fid_stdout.write('=====================\n')
		fid_stdout.write('Extraction for site %s\n'%stanam )

		cursor = self.db.cursor()  
		sites_info,Ek = \
		    im_gen( cursor,sid,rid,stanam,\
		            sites_info,Ek,fid_stdout,\
		            rup_model_ids = rup_model_ids,bb=self.bb)
		cursor.close()

		fid_stdout.write('Finish extraction for site %s\n'%stanam )

	    # save to *.py for further use
	    meta = dict( 
	        sites_info = sites_info,
	        im11 = Ek
	    )
	    header = '# Meta file for Source %s,%s \n'%(sid,rid)
	    save( metafile, meta, header = header )

	    # turn meta1 dictionary key to attributes to use below without load meta file
	    meta = namespace(meta)

	    # close what left
	    fid_stdout.close()

	else:
	    # if metafile exists, just load ( major time consuming and memory occupation )
	    # this will cause the sites_info key problem by load function
	    meta = load( metafile )


	if not RupVar:
	    return meta1, meta

	else:
	    # consider the rupture variance associated with one rupture and one hypocenter
	    # This is very important! 
	    # when deal with two CyberShake models, to keep using the same site set, you can select right here and make a unified site_info 
	    # dictionary for all later operations (to save slip variability for example) 
	    sites_run = meta.sites_info.keys()   # key to get site location and name corresponding for each sources (keep using this for the rest of workflow)

	    periods = meta.im11['periods']
	    Nsta = len(sites_run)
	    stanam = sites_run[0]
	    Nh, Nf, Nt = np.array( meta.im11[stanam] ).shape

	    ih_start = 0
	    ih_end = Nh 

	    E_xfs = {}
	    for ip in xrange( len(periods) ):
		Ti = periods[ip]
		Tkey = '%.2f'%Ti 
		tmp = np.zeros( (Nh,Nf,Nsta) )
		for ista in xrange( Nsta ):
		    stanam = meta.sites_info[sites_run[ista]]['name']
		    tmp[:,:,ista] = np.array(meta.im11[stanam])[:,:,ip]
		E_xfs[Tkey] = tmp

	    return E_xfs, meta1, meta


    def extract_ruptures(self,sid, meta_save=True):
	"""
	Extract rupture info and save as metadata
	"""
	rups_info = {}

	if meta_save == True:

	    metafile1 = self.RupsMeta + \
	        'meta_rup_%s_%s_%s.py'%(self.erf_id,self.rup_scenario_id,sid)

	    # get SQL cursor
	    cursor = self.db.cursor()
	    fid_stdout = open( './stdout', 'w' )
	    if not os.path.exists( metafile1 ):
		rups_info = rup_gen( cursor, sid, rups_info, fid_stdout, erf_id = self.erf_id, rup_scenario_id=self.rup_scenario_id )

		fid_stdout.close()
		cursor.close()

		# save rupture info into file
		meta1 = dict( 
		    rups_info = rups_info, 
		)
		header = '# Rupture meta file for rupture (%s, %s)\n'%(sid,rid)
		save(metafile1,meta1,header=header)
	else:
	    # get SQL cursor
	    cursor = self.db.cursor()
	    fid_stdout = open( './stdout', 'w' )
	    rups_info = rup_gen( cursor, sid, rups_info, fid_stdout, erf_id=self.erf_id, rup_scenario_id=self.rup_scenario_id )
	    fid_stdout.close()
	    cursor.close()

	return rups_info


    def cybershk_sites_rups(self, sid, meta_rups, meta):
	"""
	Generate GMT file for plot
	Runsites and ruptures
	"""

	# ==============
	# Rupture and Hypocenter Info
	# ==============
	rups_info = meta_rups.rups_info
	lon,lat,dep = rups_info['surface']
	fid_rup = open( self.mapin_rups+'/cybershk.%s.%s.%s'%(self.erf_id, self.rup_scenario_id,sid),'w' )  # fault trace file for GMT plot
	for ipoint in xrange( len(lon) ):
	    fid_rup.write('%s %s\n'%(lon[ipoint], lat[ipoint]))
	fid_rup.close()

	Nh, Nf, hypo_loc = rups_info['hypo_slip']

	ih_start = 0
	ih_end = Nh 

	fid=open(self.mapin_rups+'/cybershk.%s.%s.%s.hypo'%(self.erf_id, self.rup_scenario_id,sid),'w')
	for ih in range( ih_start, ih_end ):
	    hypo_key = 'h%4.4i'%ih
	    rlon, rlat, rdep = hypo_loc[hypo_key]
	    fid.write( '%s %s 10 0 0 LM %s\n'%(rlon,rlat,ih+1) )    # For GMT plot
	fid.close()

	# =====================
	# Run sites for GMT plots
	# =====================
	sites_run = meta.sites_info.keys()
	Nsta = len(sites_run)
	site_ids = []                                                   
	fid = open(self.mapin_sites+'/cybershk_sites_run.%s.%s'%(self.rup_scenario_id,sid),'w')   # run site locations
	fid1 = open(self.mapin_sites+'/cybershk_sites_run.%s.%s.pstext'%(self.rup_scenario_id,sid),'w')   # site names
	for i in xrange( Nsta ):
	    site_id = meta.sites_info[sites_run[i]]['id']
	    site_name = meta.sites_info[sites_run[i]]['name']
	    site_lon = meta.sites_info[sites_run[i]]['lon']
	    site_lat = meta.sites_info[sites_run[i]]['lat']

	    fid.write('%s %s\n'%(site_lon, site_lat) )
	    fid1.write('%s %s %s %s %s %s %s\n'%(site_lon+0.05,site_lat+0.05,10,0,0,'LM',site_name))
	fid.close()
	fid1.close()


    def rups_flatfile( self, meta1_rups, savefile=False ):
	"""
	Source info flatfile (all available ruptures)
	"""
	src_info = []
	if savefile == True:
	    flatfile_srcinfo = self.flatfile_pth + '/source_info'
	    fid_srcinfo = open( flatfile_srcinfo, 'w' )
	    header = 'SourceID RuptureID Mag FaultType HypoIndex'
	    fid_srcinfo.write( '%s\n'%header )

	for ik in xrange( len( self.sids ) ):
	    sid = self.sids[ik]
	    meta_rups = meta1_rups[ik]
	    rids = meta_rups['MR'][1]
	    Mws = meta_rups['MR'][0]

	    rake,dip,Ztor,Zbom = meta_rups['MR'][4:]
	    Nh,Nf,hypo_loc = meta_rups['hypo_slip']

	    ih_start = 0
	    ih_end = Nh 

	    Nhtmp = ''
	    for ih in range( ih_start, ih_end ):
		Nhtmp = Nhtmp + str(ih) + ' '

	    for irup in xrange( len(rids) ):
		src_info.append( [sid,rids[irup],Mws[irup],rake,dip,Ztor,Zbom,Nh,ih_start,ih_end] )

		if savefile == True:
		    fid_srcinfo.write( '%s %s %s %s %s %s %s %s\n'%(sid,rid,Mw,rake,dip,Ztor,Zbom,Nhtmp) )

	if savefile == True:
	    fid_srcinfo.close()

	return src_info


    def LocalModels(self,sites,meta):

	# sites comes form cybershk_sites function
	# ... use meta to specify the interpolation
	sites_run = meta.keys()
	z25_0 = []; z1_0 = []; vs30_0 = []
	slon0 = []; slat0 = []
	for i in xrange( len(sites_run) ):
	    stanam = meta[sites_run[i]]['name']
	    slon0.append( meta[stanam]['lon'] )
	    slat0.append( meta[stanam]['lat'] )
	    vs30_0.append( sites[stanam]['Vs30_Wills_2006'] )
	    z25_0.append( sites[stanam]['Z2.5'] )
	    z1_0.append( sites[stanam]['Z1.0'] )
	slon0 = np.array( slon0, 'f' )
	slat0 = np.array( slat0, 'f' )
	vs30_0 = np.array( vs30_0, 'f' ) 
	z25_0 = np.array( z25_0, 'f' )
	z1_0 = np.array( z1_0, 'f' )

	# Do the interpolation and save into file for GMT plot
	vs30 = interp( slon0, slat0, vs30_0, self.slon1d, self.slat1d, eps=self.eps, method=self.method )
	z25  = interp( slon0, slat0, z25_0, self.slon1d, self.slat1d, eps=self.eps, method=self.method )
	z1   = interp( slon0, slat0, z1_0, self.slon1d, self.slat1d, eps=self.eps, method=self.method )

	SiteGMTfile = self.mapin_sites + '/cybershk_site_local_models'
	fid = open( SiteGMTfile, 'w' )
	for i in xrange( len(vs30) ):
	    fid.write( '%s %s %s %s %s\n'%( self.slon1d[i], self.slat1d[i], vs30[i], z25[i], z1[i] ) )
	fid.close()



    def OpenSHA_nga_cpt( self, periods, rids, SiteName=None ):
	rids1 = []
	for ik in xrange( len( self.sids ) ):
	    rids1.append( rids[ik][-1] )   # compute the last rupture for each source (the largest one)
	cpt_OpenSHA_nga(self.scripts, self.NGAmeta, self.sids, rids1, periods, SiteName=SiteName, erf_id = self.erf_id )
	return 1


    # use OpenSHA NGA results (not use this one)
    def nga_OpenSHA(self,sid,rid,meta):
	"""
	Compute NGA from OpenSHA database for given rupture set
	"""

	OpenSHA_output = OpenSHA_nga_files(self.NGAmeta,sid,rid, self.Ti,erf_id=35)

	sites_run = meta.keys()
	Ns = len(sites_run)

	ngaOpenSHA = OpenSHA_nga( self.NGAmeta, self.NGAmodel, sid, rid, self.Ti, erf_id=35 )
	ngaO = {}
	for inga,nga in enumerate( self.NGAmodel ):
	    ngaO[nga] = []
	    for ista in xrange( Ns ):
		stanam = sites_run[ista]
		ngaO[nga].append(ngaOpenSHA[nga][stanam])    # rearange the site order in terms of metafile for CS
	    ngaO[nga] = np.exp( np.array( ngaO[nga] ) )

	return ngaO


    # extract site info from OpenSHA 
    def sites_flatfile( self, sid, meta, meta_rup, savefile=False):
	"""
	Write for site flatfile for given rupture set based on OpenSHA
	return dictionary for the sites 
	# This could be computed internally using the tools you developed
	"""
	rid = meta_rup[1][-1] 
	OpenSHA_output = OpenSHA_nga_files(self.NGAmeta, sid, rid, 3.0, erf_id=35)

	sites_flat_str = {}
	sites_flat = {}
	lines = open(OpenSHA_output, 'r').readlines()
	for il in range( 1, len(lines) ):                                                 
	    spl = lines[il].strip().split(',')
	    stanam = spl[1]

	    if stanam in self.sites.keys():
		# further selection from OpenSHA outputs to get the exact site set as inputed (the very first step)

		# attention to basin depth (depend on CVM models) 
		Z25 = float(self.sites[stanam]['Z2.5'])
		Z10 = float(self.sites[stanam]['Z1.0'])*1000.    # in meter
		Vs30 = float(self.sites[stanam]['Vs30_Wills_2006'])

		#  ID,  ShortName, Vs30 (m/s), Z2.5 (km),     Z1.0 (m),     Rjb,     Rrup,    Rx
		sites_flat_str[stanam] = ' '.join((spl[0],stanam,str(Vs30), str(Z25), str(Z10),spl[-12],spl[-11],spl[-9]))
		#sites_flat_str[stanam] = ' '.join((spl[0],stanam,spl[-18],spl[-17],spl[-15],spl[-12],spl[-11],spl[-9]))

		# ID Vs30 (m/s), Z2.5 (km),     Z1.0 (m),     Rjb,     Rrup,    Rx
		#sites_flat[stanam] = [  int(spl[0]),float(spl[-18]),float(spl[-17]),float(spl[-15]),\
		    #                  float(spl[-12]),float(spl[-11]),float(spl[-9]) ]
		sites_flat[stanam] = [  int(spl[0]), Vs30, Z25, Z10,\
		                        float(spl[-12]),float(spl[-11]),float(spl[-9]) ]

	# write flatfile (in the order of sites_run)
	sites_run = meta.keys()
	Ns = len(sites_run)

	sites_flat1 = []
	# get the order correct
	for ista in xrange(Ns):
	    stanam = sites_run[ista]
	    sites_flat1.append( sites_flat[stanam] )
	del( sites_flat )

	if savefile == True:

	    prefix = 'SourceID%s_RuptureID%s'%(sid,0)
	    flatfile = self.flatfile_pth+'/flatfile_%s/'%prefix
	    if not os.path.exists( flatfile ):
		os.mkdir( flatfile )
	    flatfile_siteinfo = flatfile + 'sites_info'

	    fid_sites = open( flatfile_siteinfo, 'w' )
	    #header = ' '.join(('ID','ShortName','Vs30 (m/s)','Z2.5 (km)','Z1.0 (m)','Rjb (km)','Rrup (km)','Rx (km)'))
	    header = ' '.join(('ID','ShortName','Vs30','Z2.5','Z1.0','Rjb','Rrup','Rx'))
	    fid_sites.write('%s\n'%header)
	    for ista in xrange( Ns ):
		stanam = sites_run[ista]
		fid_sites.write('%s\n'%sites_flat_str[stanam])
	    fid_sites.close()

	return np.array(sites_flat1)


    def nga_Py(self, sid, meta_rup, meta, sites_flat=None, ref=0, Ts=None, model_name='BA', NGAdict=None):
	"""
	Compute NGA model, or reference models for given rupture set 
	Using Python NGA classes to compute
	and return Ekxs which has length Ns
	"""
	# commonly used predictors in all NGA models
	Mws = meta_rup[0]
	rake,dip,Ztor,Zbom = meta_rup[4:]

	W = (Zbom-Ztor)/np.sin(dip*np.pi/180.)   #
	
	modelVersion = self.ngaModelVersion 

	if modelVersion == '2014': 
	    # check the application of pynga for NGA2
	    if model_name == 'BA': 
		model_name = 'BSSA'
	    if model_name == 'AS': 
		model_name = 'ASK'	    
	print 'NGA model: ', modelVersion 
	print model_name
	# sitesinfo
	if sites_flat == None: 
	    sites_flat = self.sites_flatfile(sid,meta,meta_rup)

	tmp1 = np.array( sites_flat )
	Vs30 = list( tmp1[:,1] )   # in m/s
	Z25 = list( tmp1[:,2] )    # in km
	Z10 = list( tmp1[:,3] )    # in m
	Rjb = list( tmp1[:,4] )
	Rrup = list( tmp1[:,5] )
	Rx = list( tmp1[:,6] )
	Fhw = 1*(Rx>0) + 0*(Rx<=0) 

	del( tmp1 )

	if ref == 0:
	    if NGAdict == None:
		dict1 = self.NGAs
	    else:
		dict1 = NGAdict # for updating purposes
	else:
	    dict1 = self.Reference

	if Ts != None: 
	    periods = Ts 
	else: 
	    periods = self.periods 

	# use the new NGA utils (combinations)
	ngaP = {}
	if modelVersion == '2008':
	    for ip in xrange( len(periods) ):
		Ti = periods[ip] 
		Tkey = '%.2f'%Ti    # attention to the 0.075
		ngaP[Tkey] = []
		for irup in xrange( len(Mws) ):
		    Mw = Mws[irup]
		    median, std, tau, sigma = NGA08(model_name, Mw, Rjb, Vs30, Ti, rake=rake,Mech=None,NGAs=dict1, \
			                            Rrup=Rrup, Rx=Rx, dip=dip,W=W,Ztor=Ztor,Z25=Z25,Z10=Z10,Fas=0,AB11=None,VsFlag=0)
		    
		    ngaP[Tkey].append( [list(median), list(np.log(std)), list(np.log(tau)), list(np.log(sigma))] )
	else: 
	    for ip in xrange( len(periods) ):
		Ti = periods[ip] 
		Tkey = '%.2f'%Ti    # attention to the 0.075
		ngaP[Tkey] = []
		for irup in xrange( len(Mws) ):
		    Mw = Mws[irup]		    
		    median, std, tau, sigma = NGA14(model_name, Mw, Rjb, Vs30, Ti, rake=rake,Mech=None,NGAs=dict1, \
			                            Rrup=Rrup, Rx=Rx, Fhw=Fhw, dip=dip,W=W,Ztor=Ztor,Z25=Z25,Z10=Z10,Fas=0,VsFlag=0)
		    
		    ngaP[Tkey].append( [list(median), list(np.log(std)), list(np.log(tau)), list(np.log(sigma))] )
	    
	return ngaP   # [Tkey][irup][ista]


    def SC08_cpt( self, rids ):

	rids1 = []
	for i in xrange( len(self.sids) ): 
	    rids1.append(rids[i][-1])
	cpt_SC08(self.sids, rids1, self.CSmeta, self.RupsMeta, self.cybershk_sitedata, \
	         self.scripts, self.fDcpt_pth, self.util, self.kwds,  self.dimS, \
	         Tlist=self.Tlist, NGAmodel=self.NGAmodel,\
	         rup_model_ids = self.rup_model_ids, mflag=self.mflag )
	return 1


    def D_flatfile( self, sid, meta, meta_rup, savefile=False ):
	"""
	Directivity Flatfile (similar to sites_flat)
	"""
	# use the computed SC (use original SC model) (just use rid=0)
	rid = meta_rup[1][-1]
	SC08_outputs = SC08_files( sid, rid, self.fD_output, 'BA', erf_id = self.erf_id, rup_scenario_id = self.rup_scenario_id )
	Nh = len(SC08_outputs)   # now Nh should include all hypocenters

	sites_run = meta.keys()
	Ns = len(sites_run)

	# ID and SiteName dict map
	#print sid,rid
	#site_ids = []; 
	sites_id1 = {}
	for i in xrange( len(sites_run) ):
	    site_id = meta[sites_run[i]]['id']
	#   site_ids.append( site_id )
	    if site_id == '297':
		print sites_run[i]
	    sites_id1['%s'%site_id] = sites_run[i]   # id and name

	D_flat0 = {}; D_flat_str = {}
	for ih in xrange( Nh ):
	    D_flat0['%s'%ih] = {}
	    D_flat_str['%s'%ih] = {}
	    lines = open(SC08_outputs['%s'%ih], 'r').readlines()
	    for il in range( 1,len(lines) ):
		spl = lines[il].strip().split(',')                                                        
		#site_id = int(spl[0])  

		site_id = spl[0]
		if sites_id1.has_key(site_id):
		    stanam = sites_id1[site_id]
		else: 
		    continue

		#if site_id in site_ids:
		#    stanam = sites_id1['%s'%site_id]
		    # keep the original c' here in the flatfile 
		    # but when you compute reference NGA model and target NGA models, you should change the value of c' for different purposes

		D_flat0['%s'%ih][stanam] = [float(spl[3]),float(spl[5]),float(spl[6]),float(spl[7]),float(spl[8]),float(spl[9])]   # return (Rrup,Rfn,Rfp,s,h,c')
		D_flat_str['%s'%ih][stanam] = ' '.join((stanam,spl[3],spl[5],spl[6],spl[7],spl[8],spl[9]))

	D_flat = {}
	for ih in xrange( Nh ):
	    # all hypocenter
	    D_flat['%s'%ih] = []
	    for ista in xrange( Ns ):
		stanam = sites_run[ista]
		D_flat['%s'%ih].append(D_flat0['%s'%ih][stanam])
	    D_flat['%s'%ih] = np.array( D_flat['%s'%ih] )

	del( D_flat0 )

	if savefile == True:
	    prefix = 'SourceID%s_RuptureID%s'%(sid,rid)
	    flatfile = self.flatfile_pth+'/flatfile_%s/'%prefix

	    if not os.path.exists( flatfile ):
		os.mkdir( flatfile )
	    for ih in xrange( Nh ):
		flatfile_Dinfo = flatfile + 'D_info_hypo%s'%ih
		fid_D = open( flatfile_Dinfo, 'w' )
		header = ' '.join(('ShortName','Rrup','Rfn','Rfp','s','h','ctildepr'))
		fid_D.write('%s\n'%header)
		for ista in xrange( Ns ):
		    stanam = sites_run[ista]
		    fid_D.write('%s\n'%D_flat_str['%s'%ih][stanam])
		fid_D.close()

	return D_flat


    def directivity_SC08(self, sid, Mws, Nh, meta, meta_rup, D_flat=None, ref=0, Ts=None, model_name='BA', NGAdict=None, IDPcpt=False):
	"""
	Extract from SC08 files for the directivity term (correction)
	or use updated model to compute directivity correction 
	"""

	sites_run = meta.keys()
	Ns = len(sites_run)

	if ref == 0:
	    # compute Non-reference model
	    if NGAdict == None:
		dict1 = self.NGAs
	    else:
		dict1 = NGAdict 
	else:
	    # compute reference model
	    dict1 = self.Reference

	if Ts != None: 
	    periods = Ts
	else: 
	    periods = self.periods 

	ih_start = 0
	ih_end = Nh 

	# D flatfile
	if D_flat == None: 
	    D_flat = self.D_flatfile( sid, meta, meta_rup )

	Rrup = {}; ctildepr = {}; s = {}; h = {}; Rfn = {}; Rfp = {}
	for ih in range( ih_start, ih_end ):
	    hkey = '%s'%ih
	    Rrup[hkey] = []; ctildepr[hkey] = []; s[hkey] = []; h[hkey] = []; Rfn[hkey] = []; Rfp[hkey] = []

	    tmp1 = D_flat[hkey]    
	    Rrup[hkey] = list( tmp1[:,0] )
	    Rfn[hkey] = list( tmp1[:,1] )
	    Rfp[hkey] = list( tmp1[:,2] )
	    s[hkey] = list( tmp1[:,3] )
	    h[hkey] = list( tmp1[:,4] )
	    ctildepr[hkey] = list( tmp1[:,5] )
	del( tmp1 )

	NewCoefs = dict1['SC']['NewCoefs'][model_name]
	cuteps = dict1['SC']['predictors']['cuteps']
	if NewCoefs == None and  cuteps== None:
	    #print 'Read from Matlab outputs'
	    # not available

	    # read fD directly from the fD_compute/outputs 
	    SC08_outputs = SC08_files( sid, meta_rup[1][-1], self.fD_output, model_name, erf_id = self.erf_id, rup_scenario_id = self.rup_scenario_id )    # has all hypocenters

	    sites_id1 = {}
	    for i in xrange( Ns ):
		key = '%s'%(meta[sites_run[i]]['id'])
		sites_id1[key] = sites_run[i]

	    for it in xrange( len(self.Tlist) ):
		if self.Ti == self.Tlist[it]:
		    index_t = it
		    break

	    SC08 = SC08_model( model_name+'08' )    # get nga directivity class instance
	    ind = (np.array(SC08.periods) == self.Ti).nonzero()[0]    # match the period
	    SC_a0 = SC08.a0s[ind]

	    ngaD = []
	    for ih in range( ih_start, ih_end ):
		ngaD.append([])
		lines = open(SC08_outputs['%s'%ih], 'r').readlines()
		tmp1 = {}
		for il in range( 1,len(lines) ):
		    spl = lines[il].strip().split(',')                                                        
		    site_id = spl[0]  
		    stanam = sites_id1[site_id]
		    tmp1[stanam] = np.exp(float(spl[index_t-len(self.Tlist)])+SC_a0)

		for ista in xrange( Ns ):
		    stanam = sites_run[ista]
		    ngaD[ih-ih_start].append( tmp1[stanam] )


	if NewCoefs != None or cuteps != None:
	    #print 'Use nga.SC08_model to compute SC'
	    ngaD = {} 
	    if IDPcpt: 
		IDP = []
	    for ip in xrange( len(periods) ):
		Ti = periods[ip]
		Tkey = '%.2f'%Ti 
		ngaD[Tkey] = []
		for irup in xrange( len(Mws) ):
		    Mw = Mws[irup]
		    ngaD[Tkey].append( [] )
		    for ih in range( ih_start, ih_end ):
			hkey = '%s'%ih
			SC08 = SC08_model(model_name+'08',cuteps=cuteps )
			kwds = {'NewCoefs':NewCoefs}
			try: 
			    mapfunc( SC08, Mw, Rrup[hkey], ctildepr[hkey], s[hkey], h[hkey], Rfn[hkey], Rfp[hkey], Ti,**kwds )
			except: 
			    print kwds 

			ngaD[Tkey][irup].append( SC08.fD )
			if ip == 0 and irup == 0 and IDPcpt:
			    IDP.append( SC08.IDP )
	if IDPcpt:
	    return ngaD, IDP  # ngaD dim: [Tkey][irup][ihypo][ista]; IDP dim: [ihypo][ista]
	else: 
	    return ngaD


    # ========================
    # Applications below:
    # like main function
    # ========================
    def IM_cpt(self, NGAcpt=False, MBF_S=False, MBF_D=False,MBF_C=False):

	CyberShakeRvar = []
	Sources = []

	if NGAcpt:
	    ngaP = {}; ngaD = {}
	    for inga, nga in enumerate( self.NGAmodel ):
		ngaP[nga] = []
		ngaD[nga] = []

	meta1_sites = []; meta1_rups=[]; Nhs = []
	sids0, rids = RupSelect( self.sids, self.cybershk_database_psword )
	for ik in xrange( len(self.sids) ):
	    sid = self.sids[ik]
	    Sources.append(sid)

	    BlockName = 'Working on Source %s...'%sid
	    start_time = HourMinSecToSec(BlockName=BlockName)
	    tmp = []
	    for rid in rids[ik]:
		print 'Rupture %s'%rid
		E_xfs, meta_rup, meta = self.extract_database( sid, rid, RupVar=True )
		tmp.append( E_xfs )
		if rid == rids[ik][-1]:
		    meta1_sites.append(meta.sites_info)
		    meta1_rups.append(meta_rup.rups_info['MR'])

	    CyberShakeRvar.append( tmp )

	    # use the last rupture
	    Nh = meta_rup.rups_info['hypo_slip'][0]

	    if NGAcpt:
		for inga, nga in enumerate( self.NGAmodel ):
		    ngaP[nga].append(self.nga_Py(sid, meta_rup.rups_info['MR'], meta.sites_info, model_name=nga))
		for inga, nga in enumerate( self.NGAmodel ):
		    ngaD[nga].append(self.directivity_SC08(sid, meta_rup.rups_info['MR'][0],Nh, meta.sites_info, meta_rup.rups_info['MR'], model_name=nga) )
	    # ngaP [nga][ik][Tkey][ir][ista] 
	    # ngaD [nga][ik][Tkey][ir][ih][ista] 

	    end_time = HourMinSecToSec()
	    hour,min,sec = SecToHourMinSec(end_time-start_time,BlockName='Source %s'%(sid))

	if MBF_S: 
	    nga = 'CB' 
	    T = 3.0
	    CBnga = CB08_nga() 

	    Fk = [] 
	    for ik in xrange(len(self.sids) ) : 
		sid = self.sids[ik] 

		# information for each source
		meta_rup = meta1_rups[ik] 
		meta = meta1_sites[ik] 

		Mws = meta_rup[0]
		rake,dip,Ztor,Zbom = meta_rup[4:]

		W = (Zbom-Ztor)/np.sin(dip*np.pi/180.)   #
		sites_flat = self.sites_flatfile(sid,meta,meta_rup)

		tmp1 = np.array( sites_flat )
		Vs30 = list( tmp1[:,1] )   # in m/s
		Z25 = list( tmp1[:,2] )    # in km
		Z10 = list( tmp1[:,3] )    # in m
		Rjb = list( tmp1[:,4] )
		Rrup = list( tmp1[:,5] )
		Rx = list( tmp1[:,6] )
		del( tmp1 )
		for ir in xrange( len(Mws) ):
		    Mw = Mws[ir]
		    kwds = {'Rrup':Rrup,'Ztor':Ztor,'dip':dip,'Z25':Z25,'W':W}
		    # site independent term
		    # ...


	if MBF_D: 
	    # organize informat
	    D = []; Nsta1 = []
	    for ik in xrange( len(self.sids) ): 
		D.append(ngaD['BA'][ik]['3.00']) # [ir][ih][ista]    r,k,x
		Nsta1.append( len(meta1_sites[ik].keys()) )
	    Nsta = np.min( Nsta1 )

	    # calculate D_x given ik and ir 
	    Dks = []
	    for ik in xrange( len(self.sids) ): 
		Nr = len( D[ik] ) 
		for ir in xrange( Nr ): 
		    tmp1 = np.log( np.array( D[ik][ir] )  )
		    Nh = tmp1.shape[0] 
		    tmp = np.average( tmp1, axis=0, weights = np.ones(Nh)*1.0/Nh ) 
		    Dks.append( tmp )   # [ik+ir][ista]

	    # calculate Ds
	    Ds = []
	    for ista in xrange( Nsta ): 
		ave = 0.0 
		ielem = 0
		for ik in xrange( len(self.sids) ): 
		    Nr = len( D[ik] ) 
		    for ir in xrange( Nr ): 
			ave += Dks[ielem][ista]  
			ielem += 1 
		ave = ave*1.0 / ielem 
		Ds.append( ave )

	    # calculate Dks-Ds 
	    dks = []; ielem = 0
	    for ik in xrange( len(self.sids) ): 
		Nr = len(D[ik]) 
		for ir in xrange( Nr ): 
		    tmp = [] 
		    for ista in xrange( Nsta ): 
			tmp.append( Dks[ielem][ista] - Ds[ista] ) 
		    ielem += 1
		    dks.append( tmp ) # ik,ir,ista

	    # save into file for GMT plot 
	    Gkxmfs0_pth= self.output_pth+'/Gkxmfs0'      # path to save the un-interpolated results
	    dks0_pth = Gkxmfs0_pth + '/Cks_MBF' 
	    if not os.path.exists( dks0_pth ): 
		os.mkdir( dks0_pth ) 

	    # save the first magnitude ir == 0
	    ielem = 0
	    for ik in xrange( len(self.sids) ): 
		meta = meta1_sites[ik]
		sites_run = meta.keys()
		sid = self.sids[ik]
		fid = open( dks0_pth + '/CyberShake.NGAs.3.00.Source%s.dks'%(sid), 'w' ) 
		for ir in xrange( Nr ): 
		    if ir == 0: 
			for ista in xrange( Nsta ): 
			    slon = meta[sites_run[ista]]['lon']
			    slat = meta[sites_run[ista]]['lat']
			    fid.write( '%s %s %s\n'%(slon,slat,dks[ielem][ista]) )
		    ielem += 1 
		fid.close()

	if MBF_C: 
	    nga = 'CB' 
	    T = 3.0
	    CBnga = CB08_nga() 

	    Fk = [] 
	    for ik in xrange(len(self.sids) ) : 
		sid = self.sids[ik] 

		# information for each source
		meta_rup = meta1_rups[ik] 
		meta = meta1_sites[ik] 

		Mws = meta_rup[0]
		rake,dip,Ztor,Zbom = meta_rup[4:]

		W = (Zbom-Ztor)/np.sin(dip*np.pi/180.)   #
		sites_flat = self.sites_flatfile(sid,meta,meta_rup)

		tmp1 = np.array( sites_flat )
		Vs30 = list( tmp1[:,1] )   # in m/s
		Z25 = list( tmp1[:,2] )    # in km
		Z10 = list( tmp1[:,3] )    # in m
		Rjb = list( tmp1[:,4] )
		Rrup = list( tmp1[:,5] )
		Rx = list( tmp1[:,6] )
		del( tmp1 )
		for ir in xrange( len(Mws) ):
		    Mw = Mws[ir]
		    kwds = {'Rrup':Rrup,'Ztor':Ztor,'dip':dip,'Z25':Z25,'W':W}
		    # path effect terms 
		    # ...




	if NGAcpt:
	    return CyberShakeRvar, ngaP, ngaD, meta1_sites, meta1_rups, Sources
	else: 
	    return CyberShakeRvar, meta1_sites, meta1_rups, Sources



    def GkxmfsAnalysis(self, CyberShakeRvar, sites_info, rups_info, Sources, srcPDF=None, Debug=False, \
                       ngaP0k=None, ngaD0k=None, hypoPDF={'CyberShake':{'Beta':(1.0,1.0),},'CB':1,'BA':1,'CY':1,'AS':1}, Ref='BA00',MBF=False):
	# MBF: Test the correspondance between ABF and MBF 

	# =============
	# PATH info
	# =============
	Gkxmfs0_pth= self.output_pth+'/Gkxmfs0'      # path to save the un-interpolated results

	# slip
	fkxmfs0_pth = Gkxmfs0_pth + '/Fkxmfs' 
	SFkxms0_pth = Gkxmfs0_pth + '/SFkxms'
	SFks0_pth = Gkxmfs0_pth + '/SFks'

	# magnitude
	ekxms0_pth = Gkxmfs0_pth + '/Ekxms' 
	SEkxs0_pth = Gkxmfs0_pth + '/SEkxs'
	SEks0_pth = Gkxmfs0_pth + '/SEks'

	# directivity
	dkxs0_pth = Gkxmfs0_pth + '/Dkxs' 
	SDks0_pth = Gkxmfs0_pth + '/SDks' 

	# Distance
	cks0_pth = Gkxmfs0_pth + '/Cks' 

	# site effects
	bs0_pth = Gkxmfs0_pth + '/Bs'

	# average level
	a0_pth = Gkxmfs0_pth + '/A'

	Gkxmfs_pth= self.output_pth+'/Gkxmfs'       # path to save the interpolated results

	# slip
	fkxmfs_pth = Gkxmfs_pth + '/Fkxmfs' 
	SFkxms_pth = Gkxmfs_pth + '/SFkxms'
	SFks_pth = Gkxmfs_pth + '/SFks'

	# magnitude
	ekxms_pth = Gkxmfs_pth + '/Ekxms' 
	SEkxs_pth = Gkxmfs_pth + '/SEkxs'
	SEks_pth = Gkxmfs_pth + '/SEks'

	# directivity
	dkxs_pth = Gkxmfs_pth + '/Dkxs' 
	SDks_pth = Gkxmfs_pth + '/SDks' 

	# Distance
	cks_pth = Gkxmfs_pth + '/Cks' 

	# basin
	bs_pth = Gkxmfs_pth + '/Bs'
	a_pth = Gkxmfs_pth + '/A'

	for f in [Gkxmfs0_pth, fkxmfs0_pth, SFks0_pth, SFkxms0_pth, 
	          ekxms0_pth, SEkxs0_pth, SEks0_pth, 
	          dkxs0_pth, SDks0_pth, 
	          cks0_pth, bs0_pth, a0_pth,
	          Gkxmfs_pth, fkxmfs_pth, SFks_pth, SFkxms_pth, 
	          ekxms_pth, SEkxs_pth, SEks_pth, 
	          dkxs_pth, SDks_pth, 
	          cks_pth, bs_pth, a_pth, ]:
	    if not os.path.exists(f):
		os.mkdir(f)

	if Debug: 
	    debug_pth = self.output_pth  + '/Debug' 
	    debug_hkxs_pth = debug_pth + '/Hkxs' 
	    for f in [debug_pth, debug_hkxs_pth,]: 
		if not os.path.exists( f ): 
		    os.mkdir( f ) 

	# ==================================
	# set up weighting functions (once)
	# ==================================
	Nk = len( CyberShakeRvar )
	print 'Source Number: ', Nk 

	Nstas = []

	Nhs = []
	pdf_f = []; pdf_m = []
	Tkey = '%.2f'%3.0   # just use one period
	for ik in xrange( Nk ):
	    Nm = len( CyberShakeRvar[ik] )
	    tmpCS = []
	    for irup in xrange( Nm ):
		tmpCS.append( np.log( CyberShakeRvar[ik][irup][Tkey] ) )   # [irup][ih,islip,ista]
	    tmp_CS = np.array( tmpCS )
	    Nm, Nh, Nf, Nsta = tmp_CS.shape

	    Nstas.append( Nsta )
	    print 'The source %s has %s sites'%(Sources[ik],Nsta)

	    del(tmpCS)
	    Nhs.append(Nh)  
	    pdf_f.append([])
	    pdf_m.append([])

	    # slip-distribution weighting function:
	    for ivar in xrange( Nf ):
		pdf_f[ik].append( 1./Nf ) 

	    # magnitude distribution
	    Mws, rids, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
	    AreaTmp = Nrow*Ncol  # in km
	    mu = 3.87 + np.log10(AreaTmp) * 1.05   # somerville 2006 M-A relationship
	    std = 0.2
	    prob0 = 1./(np.sqrt(2.*np.pi)*std) * np.exp(-0.5*((Mws-mu)/std)**2)
	    pdf_m[ik] = prob0 / sum(prob0)    # weighting should be summed up to 1
	    pdf_m[ik] = pdf_m[ik].tolist()

	Nsta0 = np.min( Nstas ) 

	# hypocenter location distribution (different half width)
	keytmp = hypoPDF['CyberShake'].keys()[0]
	pdf_x = []; Nsigma = 1

	# shape parameters (different for different distributions)
	sigmas = hypoPDF['CyberShake'][keytmp]  # for keytmp = Guassian: one sigma; for keytmp = Beta: (alpha, beta)
	try: 
	    Nsigma = len(sigmas) 
	except: 
	    Nsigma = 1 
	    sigmas = [sigmas,]

	if keytmp == 'Gaussian':
	    for isigma in xrange( Nsigma ):
		pdf_x.append([]) 
		for ik in xrange( Nk ):
		    Nh = Nhs[ik]   # different source has different hypocenter numbers
		    if sigmas[isigma] != 0.0: 
			xh = np.arange( Nh ) / float(Nh-1)
			mu = 0.5
			sigma = sigmas[isigma]   # input sigma should be in normalized unit (0,1), then it will be the same for each source!
			prob0 = np.exp(-0.5*((xh-mu)/sigma)**2)   # normalized Gaussian distribution
			pdf_x0 = (prob0/sum(prob0)).tolist()    # to make sure sum(prob) = 1
			pdf_x[isigma].append( pdf_x0 )
		    else: 
			# uniform distribution for hypocenters 
			pdf_x[isigma].append( np.repeat( 1./Nh, Nh).tolist() )

	elif keytmp == 'Beta': 
	    for isigma in xrange( Nsigma ): 
		pdf_x.append([])
		for ik in xrange( Nk ):
		    Nh = Nhs[ik]   # different source has different hypocenter numbers
		    xh = (np.arange( Nh )+1.0) / (Nh+1.0)
		    alpha, beta = sigmas[isigma]   # input sigma should be in normalized unit (0,1), then it will be the same for each source!
		    prob0 = stats.beta.pdf(xh,alpha,beta)     # use scipy.stats
		    pdf_x0 = (prob0/sum(prob0)).tolist()    # to make sure sum(prob) = 1
		    pdf_x[isigma].append( pdf_x0 )

	# source, site distribution
	if srcPDF == None: 
	    pdf_k = np.ones(Nk)/Nk
	else: 
	    print 'Use input srcPDF'
	    pdf_k = srcPDF

	pdf_s0 = np.ones(Nsta)/(1.0*Nsta) 
	pdf_s = np.ones(self.Nloc)/self.Nloc

	# =================
	# ABF analysis
	# =================
	for it, Ti in enumerate( self.periods ): 

	    Tkey = '%.2f'%Ti
	    print 'Analysis at periods = %s sec'%(Tkey)

	    print '='*30
	    start_time0 = HourMinSecToSec(BlockName='Prepare CyberShake, Ref, and NGAs...')
	    print '='*30

	    # 1. CyberShake G and decomposition using ABF
	    print 'CyberShake preparation'
	    Gkxmfs = []
	    for ik in xrange( Nk ):
		Nm = len( CyberShakeRvar[ik] )
		tmpCS = []
		for irup in xrange( Nm ):
		    tmpCS.append( np.log( CyberShakeRvar[ik][irup][Tkey] ) )   # [irup][ih,islip,ista]
		tmp_CS = np.array( tmpCS )
		Nm, Nh, Nf, Nsta = tmp_CS.shape
		Gkxmfs.append( tmp_CS )
		del(tmpCS, tmp_CS)

	    # 2. Use pre-defined NGA models as other reference model
	    # Note: you can just generate tmp[ik,ir,ih,ista] (and do the abf with pdf_f = None to reduce the memory and calculation time) 
	    # later you can just subtract between factors and save into file 
	    print 'Reference model NGA-type preparation'
	    RefModelName = Ref[:2]    # NGA model to select terms
	    DFlag = int( Ref[2] )     # 0 or 1 for directivity
	    BFlag = int( Ref[3] )     # 0 or 1 for directivity
	    RefModel = []; IDP0 = []
	    for ik in xrange( Nk ):
		Nm, Nh, Nf, Nsta = Gkxmfs[ik].shape
		tmp = np.zeros( (Nm,Nh,Nsta) )
		sid = Sources[ik] 
		print 'Source %s'%sid
		
		ngaP = self.nga_Py( sid, rups_info[ik], sites_info[ik], ref=1, Ts = [Ti,], model_name = RefModelName )
		Mws = rups_info[ik][0]
		ngaD, IDP00 = self.directivity_SC08(sid, Mws, Nh, sites_info[ik], rups_info[ik], ref=1, Ts = [Ti,], model_name = RefModelName, IDPcpt=True)

		# compute IDP (as the parameter for directivity effects)
		IDP0.append( np.array( IDP00 ) )

		for ir in xrange( Nm ):
		    for ista in xrange( Nsta ):
			for ih in xrange( Nh ): 
			    tmp[ir,ih,ista] = np.log( ngaP[Tkey][ir][0][ista] ) + np.log( np.array(ngaD[Tkey][ir])[ih,ista] ) * DFlag
		RefModel.append( tmp )
		del(tmp)
	    
	    print time.localtime() 

	    # 3. Prepare original NGA models (with directivity controlled)
	    print 'NGA model preparation'
	    TargetModel = {}
	    for inga in xrange( len(self.NGAmodel) ):
		nga = self.NGAmodel[inga]

		# for each nga model
		TargetModel[nga] = []
		for ik in xrange( Nk ):

		    Nm, Nh, Nf, Nsta = Gkxmfs[ik].shape
		    tmp = np.zeros( (Nm,Nh,Nsta) )
		    sid = Sources[ik] 

		    if ngaP0k == None:
			print 'compute NGA model'
			ngaP_tmp = self.nga_Py( sid, rups_info[ik], sites_info[ik], Ts = [Ti,], model_name = nga )
			ngaP = ngaP_tmp[Tkey]
		    else: 
			ngaP = ngaP0k[nga][ik][Tkey]
                    print time.localtime() 

		    if ngaD0k == None: 
			print 'compute NGA directivity model'
			Mws = rups_info[ik][0]
			ngaD_tmp = self.directivity_SC08(sid, Mws, Nh, sites_info[ik], rups_info[ik], Ts = [Ti,], model_name = nga )
			ngaD = ngaD_tmp[Tkey]
		    else: 
			ngaD = ngaD0k[nga][ik][Tkey]
	            print 'complete %s calculation for Source %s'%(nga, sid)
                    print time.localtime() 

		    for ir in xrange( Nm ):
			for ista in xrange( Nsta ):
			    for ih in xrange( Nh ): 
				tmp[ir,ih,ista] = np.log( ngaP[ir][0][ista] ) + np.log( np.array(ngaD[ir])[ih,ista] ) * hypoPDF[nga]
		    TargetModel[nga].append( tmp )
		    del(tmp,ngaP,ngaD)


	    print '='*30
	    end_time0 = HourMinSecToSec(BlockName='Preparation of CyberShake, Ref, and NGAs finished')
	    hour,min,sec = SecToHourMinSec(end_time0-start_time0,BlockName='IM Preparation')
	    print '='*30 + '\n'


	    for isigma in xrange( Nsigma ):
		try: 
		    # for distributions more than two parameters 
		    if len(sigmas[isigma]) == 2: 
			print sigmas[isigma][0], sigmas[isigma][1]
			sigma_str = 'Sigma%s_%s'%('%.2f'%sigmas[isigma][0],'%.2f'%sigmas[isigma][1])
		except:
		    sigma_str = 'Sigma%s'%('%.2f'%sigmas[isigma])

		print '='*30
		start_time0 = HourMinSecToSec(BlockName='ABF analysis at %s'%sigma_str)
		print '='*30

		pdf_x0 = pdf_x[isigma] 

		# 1. Decomposition of CyberShake (revise)
		Fkxmfs0, Ekxms0, Hkxs0 = im_decom1(Gkxmfs, pdf_m, pdf_x0, pdf_kf=pdf_f)   # first step (average out)
		Dkxs0, Cks0, Bs0, A0 = im_decom2( Hkxs0, Nsta0, pdf_x0, pdf_k, pdf_s0 )

		# interpolation
		Hkxs = [] 
		Rjb_0 = []; Rrup_0 = []; Rx_0 = []
		Rjb = []; Rrup = []; Rx = []
		for ik in xrange( Nk ):
		    Nm, Nh, Nf, Nsta = Gkxmfs[ik].shape
		    meta = sites_info[ik]
		    sites_run = meta.keys()
		    slon = []; slat = []
		    for i in xrange( Nsta ):
			slon.append( meta[sites_run[i]]['lon'] )
			slat.append( meta[sites_run[i]]['lat'] )
		    slon = np.array( slon )
		    slat = np.array( slat )

		    sites_info0 = self.sites_flatfile( Sources[ik], sites_info[ik], rups_info[ik] )    # Vs30 and Z1.0, Z2.5 (uninterpolated results)

		    if ik == 0:
			# compute once (all the same for all sources) 
			# interpolate the site info (Vs30, Z2.5, Z1.0)
			# add distance info for C-map test with NGA path effects

			# independent of velocity models
			Vs30_0 = sites_info0[:,1]
			Vs30 = interp( slon, slat, sites_info0[:,1], self.slon1d, self.slat1d, eps=self.eps, method=self.method )

			Z25_0 = sites_info0[:,2]
			Z10_0 = sites_info0[:,3]
			Z25 = interp( slon, slat, sites_info0[:,2], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			Z10 = interp( slon, slat, sites_info0[:,3], self.slon1d, self.slat1d, eps=self.eps, method=self.method )

		    Rjb_0.append(sites_info0[:,4])
		    Rrup_0.append(sites_info0[:,5])
		    Rx_0.append(sites_info0[:,6])

		    tmp = interp( slon, slat, sites_info0[:,4], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    Rjb.append( tmp ) 
		    tmp = interp( slon, slat, sites_info0[:,5], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    Rrup.append( tmp ) 
		    tmp = interp( slon, slat, sites_info0[:,6], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    Rx.append( tmp ) 

		    tmp = np.zeros((Nh,self.Nloc))
		    for ih in xrange( Nh ):
			try: 
			    tmp[ih,:] = interp( slon, slat, Hkxs0[ik][ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			except: 
			    print 'Source: %s, Hypocenter: %s'%(self.sids[ik], ih)
			    for ista in xrange( Nsta ): 
				print Hkxs0[ik][ih,ista]
			    raise ValueError
		    Hkxs.append( tmp )

		    if Debug and sigmas[isigma] == 0.0 and ik == 0 and Tkey == '3.00':
			sid = self.sids[ik]

			debug_hkxs_pth1 = debug_hkxs_pth + '/%s'%sigma_str
			debug_hkxs_pth0 = debug_hkxs_pth1 + '/%s'%sid
			for f in [debug_hkxs_pth1, debug_hkxs_pth0,]:
			    if not os.path.exists(f):
				os.mkdir(f)

			# save uninterpolated
			for ih in xrange(Nh):
			    hkxsFile = debug_hkxs_pth0 + '/CyberShake.%s.Source%s.Ih%s.Hkxs0'%(Tkey, sid, ih)
			    fid = open( hkxsFile, 'w' )
			    for iloc in xrange(Nsta): 
				fid.write( '%s %s %s\n'%(\
				    slon[iloc], slat[iloc], \
				    Hkxs0[ik][ih,iloc]))
			    fid.close()
			# save interpolated
			for ih in xrange(Nh):
			    hkxsFile = debug_hkxs_pth0 + '/CyberShake.%s.Source%s.Ih%s.Hkxs'%(Tkey, sid, ih)
			    fid = open( hkxsFile, 'w' )
			    for iloc in xrange(self.Nloc): 
				fid.write( '%s %s %s\n'%(\
				    self.slon1d[iloc], self.slat1d[iloc], \
				    tmp[ih,iloc]))
			    fid.close()
		    del( tmp )
		Dkxs, Cks, Bs, A = im_decom2( Hkxs, self.Nloc, pdf_x0, pdf_k, pdf_s )

		# 2.Decomposition of the reference model 
		Ekxms0R, Hkxs0R = im_decom1(RefModel, pdf_m, pdf_x0, pdf_kf=None)   # first step (average out)
		Dkxs0R, Cks0R, Bs0R, A0R = im_decom2( Hkxs0R, Nsta0, pdf_x0, pdf_k, pdf_s0 )

		# interpolation
		HkxsR = [] 
		fD0 = []; fD = []; IDP = []
		SC08 = SC08_model(RefModelName+'08') 
		coef_b = SC08.Coefs['%.3f'%Ti]['b']
		for ik in xrange( Nk ):
		    Nm, Nh, Nf, Nsta = Gkxmfs[ik].shape
		    meta = sites_info[ik]
		    sites_run = meta.keys()
		    slon = []; slat = []
		    for i in xrange( Nsta ):
			slon.append( meta[sites_run[i]]['lon'] )
			slat.append( meta[sites_run[i]]['lat'] )
		    slon = np.array( slon )
		    slat = np.array( slat )

		    tmpIDP_fD = np.zeros((Nh,Nsta))
		    for ih in xrange( Nh ):
			# remove the mean of IDP
			tmpIDP_fD[ih,:] = IDP0[ik][ih,:] - np.average( IDP0[ik], axis=0, weights=pdf_x0[ik] )
		    fD0.append( coef_b * tmpIDP_fD )   # normalized directivity effects (zero mean over hypocenters) to compare with d factor

		    tmp = np.zeros((Nh,self.Nloc))
		    tmpIDP = np.zeros((Nh,self.Nloc))
		    tmpfD = np.zeros((Nh,self.Nloc))
		    for ih in xrange( Nh ):
			tmp[ih,:] = interp( slon, slat, Hkxs0R[ik][ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			tmpIDP[ih,:] = interp( slon, slat, IDP0[ik][ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			tmpfD[ih,:] = interp( slon, slat, tmpIDP_fD[ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    IDP.append( tmpIDP )
		    fD.append( tmpfD ) 
		    HkxsR.append( tmp )
		    del( tmp, tmpIDP, tmpfD )

		DkxsR, CksR, BsR, AR = im_decom2( HkxsR, self.Nloc, pdf_x0, pdf_k, pdf_s )

		# Decomposition of NGA models
		Ekxms0N = []; EkxmsN = []
		Dkxs0N = []; Cks0N = []; Bs0N = []; A0N = []
		DkxsN = []; CksN = []; BsN = []; AN = []
		for inga in xrange( len(self.NGAmodel) ):
		    nga = self.NGAmodel[inga]

		    # interpolation of each NGA models
		    Ekxms0N1, Hkxs0N1 = im_decom1(TargetModel[nga], pdf_m, pdf_x0, pdf_kf=None)   # first step (average out)
		    Dkxs0N1, Cks0N1, Bs0N1, A0N1 = im_decom2( Hkxs0N1, Nsta0, pdf_x0, pdf_k, pdf_s0 )

		    Ekxms0N.append( Ekxms0N1) 
		    Dkxs0N.append( Dkxs0N1 ) 
		    Cks0N.append( Cks0N1 ) 
		    Bs0N.append( Bs0N1 ) 
		    A0N.append( A0N1 )

		    # interpolation:
		    HkxsN1 = []
		    for ik in xrange( Nk ):
			Nm, Nh, Nf, Nsta = Gkxmfs[ik].shape
			meta = sites_info[ik]
			sites_run = meta.keys()
			slon = []; slat = []
			for i in xrange( Nsta ):
			    slon.append( meta[sites_run[i]]['lon'] )
			    slat.append( meta[sites_run[i]]['lat'] )
			slon = np.array( slon )
			slat = np.array( slat )

			tmp = np.zeros((Nh,self.Nloc))
			for ih in xrange( Nh ):
			    tmp[ih,:] = interp( slon, slat, Hkxs0N1[ik][ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			HkxsN1.append( tmp )
			del( tmp )
		    DkxsN1, CksN1, BsN1, AN1 = im_decom2( HkxsN1, self.Nloc, pdf_x0, pdf_k, pdf_s )
		    DkxsN.append( DkxsN1 )
		    CksN.append( CksN1 )
		    BsN.append( BsN1  )
		    AN.append( AN1 )
		print '='*30
		end_time0 = HourMinSecToSec(BlockName='ABF analysis ends')
		hour,min,sec = SecToHourMinSec(end_time0-start_time0,BlockName='ABF analysis')
		print '='*30 + '\n'

		# Before save, think about what do you need to show and make them easy to compute and plot (more general)
		# Gkxmfs to get dimension and doing other operations
		print '='*30
		start_time0 = HourMinSecToSec(BlockName='ABF results saving at %s'%sigma_str)
		print '='*30

		################################################## 
		# Save un-interpolated metadata:
		################################################## 
		SigmaS = '/%s'%sigma_str

		# So Far
		# Compute SFks, SEks, SDks and save them 
		# save Dkxs, Cks, Bs, A
		SFks0_pth1 = SFks0_pth + SigmaS
		SEks0_pth1 = SEks0_pth + SigmaS 
		dkxs0_pth1 = dkxs0_pth + SigmaS
		SDks0_pth1 = SDks0_pth + SigmaS
		cks0_pth1  = cks0_pth  + SigmaS
		bs0_pth1   = bs0_pth   + SigmaS
		a0_pth1    = a0_pth    + SigmaS
		for f in [SFks0_pth1, SEks0_pth1, \
		          dkxs0_pth1, SDks0_pth1, cks0_pth1, bs0_pth1, a0_pth1, ]:
		    if not os.path.exists(f):
			os.mkdir(f)

		if isigma == 0: 
		    #    # for slip (m,f) related quantities, they don't really depends on hypocenter average except their averaged standard deviations
		    fkxmfs0_pth1 = fkxmfs0_pth + SigmaS
		#    SFkxms0_pth1 = SFkxms0_pth + SigmaS
		#    
		    ekxms0_pth1 = ekxms0_pth + SigmaS

		#    SEkxs0_pth1 = SEkxs0_pth + SigmaS 
		    for f in [fkxmfs0_pth1, ekxms0_pth1,]:
			if not os.path.exists(f):
			    os.mkdir(f)

		# common information
		Nk = len(Gkxmfs) 
		src_fid = open( Gkxmfs0_pth + '/SourceInfo', 'w' )
		srcsh_fid = open( Gkxmfs0_pth + '/SourceRuptureHypoInfo','w' )
		for ik in xrange( Nk ):
		    srcmw_fid = open( Gkxmfs_pth + '/SourceRuptureHypoInfo','w' )
		    sid = Sources[ik]
		    Mws, rids, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
		    Nm,Nh,Nf,Nsta = Gkxmfs[ik].shape
		    src_fid.write('%s %s %s\n'%(sid,Nrow*Ncol,pdf_k[ik]))
		    srcsh_fid.write('%s %s %s %s\n'%(sid,Nm,Nh,Nf))   # used later for standard deviation calculation
		src_fid.close() 
		srcsh_fid.close() 

		# a0 factor
		fid_a = open( a0_pth1 + '/CyberShake.NGAs.%s.a'%(Tkey), 'w' )
		#fid_a.write( '# %s %s %s %s %s %s\n'%( self.NGAmodel[0], self.NGAmodel[1], self.NGAmodel[2], self.NGAmodel[3], 'CyberShake', 'RefModel' ))
		fid_a.write( '%s %s %s %s %s %s\n'%( A0-A0N[0],A0-A0N[1],A0-A0N[2],A0-A0N[3], A0, A0-A0R) )   
		fid_a.close()

		SFks0 = []; SEks0 = []
		for ik in xrange( Nk ):
		    sid = Sources[ik]
		    dkxs0_pth0 = dkxs0_pth1 + '/%s'%sid
		    SFks0_pth0 = SFks0_pth1 + '/%s'%sid
		    SEks0_pth0 = SEks0_pth1 + '/%s'%sid
		    for f in [SFks0_pth0,SEks0_pth0,dkxs0_pth0,]:
			if not os.path.exists(f):
			    os.mkdir(f)

		    meta = sites_info[ik]
		    sites_run = meta.keys()
		    slon = []; slat = []
		    for i in xrange( Nsta0 ):
			slon.append( meta[sites_run[i]]['lon'] )
			slat.append( meta[sites_run[i]]['lat'] )
		    slon = np.array( slon )
		    slat = np.array( slat )

		    Mws, rids, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
		    Nm,Nh,Nf,Nsta = Gkxmfs[ik].shape
		    if isigma == 0:
			# F-map and E-map for only one sigma (isigma=0) 
			fkxmfs0_pth0 = fkxmfs0_pth1 + '/%s'%sid 
			ekxms0_pth0 = ekxms0_pth1 + '/%s'%sid 

			for f in [fkxmfs0_pth0, ekxms0_pth0,]:
			    if not os.path.exists(f):
				os.mkdir(f)

			srcmw_fid = open( ekxms0_pth0 + '/SourceRuptureMwInfo','w' )
			for ir in xrange(Nm): 
			    srcmw_fid.write( '%s\n'%(Mws[ir]) ) 
			    fkxmfs0_pth00 = fkxmfs0_pth0 + '/M%.2f'%Mws[ir] 
			    ekxms0_pth00 = ekxms0_pth0 + '/M%.2f'%Mws[ir] 
			    for f in [fkxmfs0_pth00, ekxms0_pth00,]:
				if not os.path.exists(f):
				    os.mkdir(f)

			    # E files
			    for ih in xrange( Nh ):
				ekxmsFile = ekxms0_pth00 + '/Period%s.Hypo%s.ekxms'%(Tkey, ih)
				fid = open( ekxmsFile, 'w' )
				tmp1 = []
				for inga in xrange( len(self.NGAmodel) ): 
				    tmp1.append(Ekxms0[ik][ir,ih,:] - Ekxms0N[inga][ik][ir, ih,:])
				tmp1.append( Ekxms0[ik][ir,ih,:] ) 
				tmp1.append( Ekxms0[ik][ir, ih,:] - Ekxms0R[ik][ir, ih,:] )
				for iloc in xrange(Nsta0): 
				    fid.write( '%s %s %s %s %s %s %s %s\n'%(\
				        slon[iloc], slat[iloc], \
				        tmp1[0][iloc], \
				        tmp1[1][iloc], \
				        tmp1[2][iloc], \
				        tmp1[3][iloc], \
				        tmp1[4][iloc], \
				        tmp1[5][iloc] )) 
				fid.close()

				# F-file
				fkxmfs0_pth01 = fkxmfs0_pth00 + '/Hypocenter%s'%(ih) 
				if not os.path.exists(fkxmfs0_pth01):
				    os.mkdir(fkxmfs0_pth01)
				for islip in xrange(Nf):
				    fkxmfsFile = fkxmfs0_pth01 + '/Period%s.Slip%s.fkxmfs'%(Tkey, islip)
				    fid = open( fkxmfsFile, 'w' )
				    for iloc in xrange(Nsta0): 
					fid.write( '%s %s %s\n'%(slon[iloc], slat[iloc], Fkxmfs0[ik][ir,ih,islip,iloc]) )
				    fid.close()
			srcmw_fid.close() 

		    # Calculate standard deviation of F-Map based Fkxmfs0 (attention: variation here, not standard deviation)
		    tmp_sd_xms = np.average( Fkxmfs0[ik]**2,axis=2, weights=pdf_f[ik] )
		    tmp_sd = np.average( np.average( tmp_sd_xms, axis=0, weights=pdf_m[ik]), axis = 0, weights = pdf_x0[ik] )
		    SFks0.append( tmp_sd )   # [ik][iloc] # only for CyberShake so far
		    SFksFile = SFks0_pth0 + '/CyberShake.NGAs.%s.Source%s.Sfks'%(Tkey, sid)
		    fid = open( SFksFile, 'w' )
		    for iloc in xrange(Nsta0): 
			fid.write( '%s %s %s\n'%(\
			    slon[iloc], slat[iloc],\
			    tmp_sd[iloc] ))
		    fid.close()

		    # Calculate standard deviation of E-Map based Ekxms0 
		    tmp_sd = [] 
		    for inga in xrange(len(self.NGAmodel) ): 
			#tmp_sd_xs = np.average( (Ekxms0[ik]-Ekxms0N[inga][ik])**2,axis=0, weights=pdf_m[ik] )
			tmp_sd_xs = np.average( (Ekxms0N[inga][ik])**2,axis=0, weights=pdf_m[ik] )
			tmp = np.average( tmp_sd_xs, axis=0, weights = pdf_x0[ik] )
			tmp_sd.append( tmp )

		    # CyberShake
		    tmp_sd_xs = np.average( Ekxms0[ik]**2,axis=0, weights=pdf_m[ik] )
		    tmp = np.average( tmp_sd_xs, axis=0, weights = pdf_x0[ik] )
		    tmp_sd.append( tmp )
		    SEks0.append( tmp_sd )    # [ik][inga][iloc]

		    # RefModel
		    #tmp_sd_xs = np.average( (Ekxms0[ik]-Ekxms0R[ik])**2,axis=0, weights=pdf_m[ik] )
		    tmp_sd_xs = np.average( (Ekxms0R[ik])**2,axis=0, weights=pdf_m[ik] )
		    tmp = np.average( tmp_sd_xs, axis=0, weights = pdf_x0[ik] )
		    tmp_sd.append( tmp )

		    SEksFile = SEks0_pth0 + '/CyberShake.NGAs.%s.Source%s.Seks'%(Tkey, sid)
		    fid = open( SEksFile, 'w' )
		    for iloc in xrange(Nsta0): 
			fid.write( '%s %s %s %s %s %s %s %s\n'%(\
			    slon[iloc], slat[iloc],\
			    tmp_sd[0][iloc] ,\
			    tmp_sd[1][iloc] ,\
			    tmp_sd[2][iloc] ,\
			    tmp_sd[3][iloc] ,\
			    tmp_sd[4][iloc] ,\
			    tmp_sd[5][iloc] ))
		    fid.close()

		    # D-Map
		    tmp2 = []
		    for ih in xrange(Nh):
			dkxsFile = dkxs0_pth0 + '/CyberShake.NGAs.%s.Source%s.Ih%s.dkxs'%(Tkey, sid, ih)
			fid = open( dkxsFile, 'w' )
			# compute residual 
			tmp1 = []
			for inga in xrange( len(self.NGAmodel) ): 
			    tmp1.append(Dkxs0[ik][ih,:] - Dkxs0N[inga][ik][ih,:])
			tmp1.append( Dkxs0[ik][ih,:] ) 
			tmp1.append( Dkxs0[ik][ih,:] - Dkxs0R[ik][ih,:] )
			tmp2.append( tmp1 ) 

			for iloc in xrange(Nsta0): 
			    fid.write( '%s %s %s %s %s %s %s %s %s %s\n'%(\
			        slon[iloc], slat[iloc],IDP0[ik][ih,iloc], fD0[ik][ih,iloc], \
			        tmp1[0][iloc], \
			        tmp1[1][iloc], \
			        tmp1[2][iloc], \
			        tmp1[3][iloc], \
			        tmp1[4][iloc], \
			        tmp1[5][iloc] )) 
			fid.close()
		    tmp2 = np.array( tmp2 ) 

		    # Standard deviation of D-Map based on tmp2 [ih][inga][iloc]
		    tmp1 = np.sqrt( np.average( tmp2**2, axis=0, weights=pdf_x0[ik] ) )
		    SDksFile = SDks0_pth1 + '/CyberShake.NGAs.%s.Source%s.Sdks'%(Tkey, sid)
		    fid = open( SDksFile, 'w' )
		    for iloc in xrange(Nsta0): 
			fid.write( '%s %s %s %s %s %s %s %s\n'%(\
			    slon[iloc], slat[iloc], \
			    tmp1[0][iloc], \
			    tmp1[1][iloc], \
			    tmp1[2][iloc], \
			    tmp1[3][iloc], \
			    tmp1[4][iloc], \
			    tmp1[5][iloc] )) 
		    fid.close()

		    # C-Map  (append the Distance information)
		    tmp1 = []
		    for inga in xrange( len(self.NGAmodel) ): 
			tmp1.append(Cks0[ik]-Cks0N[inga][ik])
		    tmp1.append( Cks0[ik] ) 
		    tmp1.append( Cks0[ik] - Cks0R[ik] )
		    cksFile = cks0_pth1 + '/CyberShake.NGAs.%s.Source%s.cks'%(Tkey, sid)
		    fid_cks = open( cksFile, 'w' )
		    for iloc in xrange(Nsta0): 
			fid_cks.write( '%s %s %s %s %s %s %s %s %s %s %s\n'%(\
			    slon[iloc], slat[iloc], Rjb_0[ik][iloc], Rrup_0[ik][iloc], Rx_0[ik][iloc], \
			    tmp1[0][iloc], \
			    tmp1[1][iloc], \
			    tmp1[2][iloc], \
			    tmp1[3][iloc], \
			    tmp1[4][iloc], \
			    tmp1[5][iloc] )) 
		    fid_cks.close()

		# write into Basin maps and overall averaging values
		tmp1 = []
		for inga in xrange( len(self.NGAmodel) ): 
		    tmp1.append(Bs0-Bs0N[inga])
		tmp1.append( Bs0 ) 
		tmp1.append( Bs0 - Bs0R )
		bsFile = bs0_pth1 + '/CyberShake.NGAs.%s.bs'%(Tkey)
		fid_bs = open( bsFile, 'w')
		for iloc in xrange( Nsta0 ):
		    fid_bs.write( '%s %s %s %s %s %s %s %s %s %s %s\n'%(\
		        slon[iloc], slat[iloc], Vs30_0[iloc], Z10_0[iloc], Z25_0[iloc], \
		        tmp1[0][iloc], \
		        tmp1[1][iloc], \
		        tmp1[2][iloc], \
		        tmp1[3][iloc], \
		        tmp1[4][iloc], \
		        tmp1[5][iloc] )) 
		fid_bs.close()


		################################################## 
		# Save interpolated metadata:
		################################################## 

		# So Far
		# Compute SFks, SEks, SDks and save them 
		# save Dkxs, Cks, Bs, A
		SFks_pth1 = SFks_pth + SigmaS
		SEks_pth1 = SEks_pth + SigmaS 
		dkxs_pth1 = dkxs_pth + SigmaS
		SDks_pth1 = SDks_pth + SigmaS
		cks_pth1  = cks_pth  + SigmaS
		bs_pth1   = bs_pth   + SigmaS
		a_pth1    = a_pth    + SigmaS
		for f in [SFks_pth1, SEks_pth1, \
		          dkxs_pth1, SDks_pth1, cks_pth1, bs_pth1, a_pth1, ]:
		    if not os.path.exists(f):
			os.mkdir(f)

		Nk = len(Gkxmfs) 
		src_fid = open( Gkxmfs_pth + '/SourceInfo', 'w' )
		srcsh_fid = open( Gkxmfs_pth + '/SourceRuptureHypoInfo','w' )
		for ik in xrange( Nk ):
		    sid = Sources[ik]
		    Mws, rids, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
		    Nm,Nh,Nf,Nsta = Gkxmfs[ik].shape
		    src_fid.write('%s %s\n'%(sid,Nrow*Ncol))
		    srcsh_fid.write('%s %s %s %s\n'%(sid,Nm,Nh,Nf))   # used later for standard deviation calculation
		src_fid.close() 
		srcsh_fid.close() 

		# a factor
		fid_a = open( a_pth1 + '/CyberShake.NGAs.%s.a'%(Tkey), 'w' )
		fid_a.write( '# %s %s %s %s %s %s\n'%( self.NGAmodel[0], self.NGAmodel[1], self.NGAmodel[2], self.NGAmodel[3], 'CyberShake', 'RefModel' ))
		fid_a.write( '%s %s %s %s %s %s\n'%( A-AN[0],A-AN[1],A-AN[2],A-AN[3], A, A-AR) )   
		fid_a.close()

		for ik in xrange( Nk ):
		    sid = Sources[ik]
		    SFks_pth0 = SFks_pth1 + '/%s'%sid
		    SEks_pth0 = SEks_pth1 
		    dkxs_pth0 = dkxs_pth1 + '/%s'%sid
		    for f in [SFks_pth0,SEks_pth0,dkxs_pth0,]:
			if not os.path.exists(f):
			    os.mkdir(f)

		    Mws, rids, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
		    Nm,Nh,Nf,Nsta = Gkxmfs[ik].shape

		    # Calculate standard deviation of F-Map based Fkxmfs0 and interpolate before save
		    tmp_sd = interp( slon, slat, SFks0[ik][:Nsta0], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    SFksFile = SFks_pth0 + '/CyberShake.NGAs.%s.Source%s.Sfks'%(Tkey, sid)
		    fid = open( SFksFile, 'w' )
		    for iloc in xrange(self.Nloc): 
			fid.write( '%s %s %s\n'%(\
			    self.slon1d[iloc], self.slat1d[iloc],\
			    tmp_sd[iloc] ))
		    fid.close()

		    # Calculate standard deviation of E-Map based Ekxms0 
		    tmp_sd = [] 
		    for inga in xrange(len(self.NGAmodel) ): 
			tmp = interp( slon, slat, SEks0[ik][inga][:Nsta0], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			tmp_sd.append( tmp )
		    tmp = interp( slon, slat, SEks0[ik][4][:Nsta0], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    tmp_sd.append( tmp )
		    tmp = interp( slon, slat, SEks0[ik][5][:Nsta0], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    tmp_sd.append( tmp )

		    SEksFile = SEks_pth0 + '/CyberShake.NGAs.%s.Source%s.Seks'%(Tkey, sid)
		    fid = open( SEksFile, 'w' )
		    for iloc in xrange(self.Nloc): 
			fid.write( '%s %s %s %s %s %s %s %s\n'%(\
			    self.slon1d[iloc], self.slat1d[iloc],\
			    tmp_sd[0][iloc] ,\
			    tmp_sd[1][iloc] ,\
			    tmp_sd[2][iloc] ,\
			    tmp_sd[3][iloc] ,\
			    tmp_sd[4][iloc] ,\
			    tmp_sd[5][iloc] ))
		    fid.close()

		    # D-Map
		    tmp2 = []
		    for ih in xrange(Nh):
			dkxsFile = dkxs_pth0 + '/CyberShake.NGAs.%s.Source%s.Ih%s.dkxs'%(Tkey, sid, ih)
			fid = open( dkxsFile, 'w' )
			tmp1 = []
			for inga in xrange( len(self.NGAmodel) ): 
			    tmp1.append(Dkxs[ik][ih,:] - DkxsN[inga][ik][ih,:])
			tmp1.append( Dkxs[ik][ih,:] ) 
			tmp1.append( Dkxs[ik][ih,:] - DkxsR[ik][ih,:] )
			tmp2.append( tmp1 ) 

			for iloc in xrange(self.Nloc): 
			    fid.write( '%s %s %s %s %s %s %s %s %s %s\n'%(\
			        self.slon1d[iloc], self.slat1d[iloc],IDP[ik][ih,iloc], fD[ik][ih,iloc], \
			        tmp1[0][iloc], \
			        tmp1[1][iloc], \
			        tmp1[2][iloc], \
			        tmp1[3][iloc], \
			        tmp1[4][iloc], \
			        tmp1[5][iloc] )) 
			fid.close()
		    tmp2 = np.array( tmp2 ) 

		    # Standard deviation of D-Map based on tmp2 [ih][inga][iloc]
		    tmp1 = np.sqrt( np.average( tmp2**2, axis=0, weights=pdf_x0[ik] ) )
		    SDksFile = SDks_pth1 + '/CyberShake.NGAs.%s.Source%s.Sdks'%(Tkey, sid)
		    fid = open( SDksFile, 'w' )
		    for iloc in xrange(self.Nloc): 
			fid.write( '%s %s %s %s %s %s %s %s\n'%(\
			    self.slon1d[iloc], self.slat1d[iloc], \
			    tmp1[0][iloc], \
			    tmp1[1][iloc], \
			    tmp1[2][iloc], \
			    tmp1[3][iloc], \
			    tmp1[4][iloc], \
			    tmp1[5][iloc] )) 
		    fid.close()

		    # C-Map  (append the Distance information)
		    tmp1 = []
		    for inga in xrange( len(self.NGAmodel) ): 
			tmp1.append(Cks[ik]-CksN[inga][ik])
		    tmp1.append( Cks[ik] ) 
		    tmp1.append( Cks[ik] - CksR[ik] )
		    cksFile = cks_pth1 + '/CyberShake.NGAs.%s.Source%s.cks'%(Tkey, sid)
		    fid_cks = open( cksFile, 'w' )
		    for iloc in xrange(self.Nloc): 
			fid_cks.write( '%s %s %s %s %s %s %s %s %s %s %s\n'%(\
			    self.slon1d[iloc], self.slat1d[iloc], \
			    Rjb[ik][iloc], Rrup[ik][iloc], Rx[ik][iloc], \
			    tmp1[0][iloc], \
			    tmp1[1][iloc], \
			    tmp1[2][iloc], \
			    tmp1[3][iloc], \
			    tmp1[4][iloc], \
			    tmp1[5][iloc] )) 
		    fid_cks.close()

		# write into Basin maps and overall averaging values
		tmp1 = []
		for inga in xrange( len(self.NGAmodel) ): 
		    tmp1.append(Bs-BsN[inga])
		tmp1.append( Bs ) 
		tmp1.append( Bs - BsR )
		bsFile = bs_pth1 + '/CyberShake.NGAs.%s.bs'%(Tkey)
		fid_bs = open( bsFile, 'w')
		for iloc in xrange( self.Nloc ):
		    fid_bs.write( '%s %s %s %s %s %s %s %s %s %s %s\n'%(\
		        self.slon1d[iloc], self.slat1d[iloc], Vs30[iloc], Z10[iloc], Z25[iloc], \
		        tmp1[0][iloc], \
		        tmp1[1][iloc], \
		        tmp1[2][iloc], \
		        tmp1[3][iloc], \
		        tmp1[4][iloc], \
		        tmp1[5][iloc] )) 
		fid_bs.close()

		print '='*30
		end_time0 = HourMinSecToSec(BlockName='ABF results saving ends')
		hour,min,sec = SecToHourMinSec(end_time0-start_time0,BlockName='ABF saving')
		print '='*30 + '\n'


    # Update model coefficients in NGA models (basin term, directivity)
    # need to be modified
    def HkxsAnalysis(self, CyberShakeRvar, sites_info, rups_info, Sources, ngaP0k=None, ngaD0k=None, \
                     hypoPDF= {'CyberShake':{'gaussian':0.6}, 'CB':1, 'BA':1, 'CY':1, 'AS':1}, \
                     NGAmodel='AS', opt='Bs', coefkey='a10'):

	# Regression
	# NGAmodel: 
	# opt: 
	# NewCoefs (just key)

	# PATH info
	Gkxmfs_pth= self.output_pth+'/GkxmfsDmap'      # path to save directivity analysis results (weighting functions)
	metapth = Gkxmfs_pth + '/Metadata' 

	# map file save pth
	dkxs_pth = Gkxmfs_pth + '/Dkxs' 
	SDks_pth = Gkxmfs_pth + '/SDks' 
	cks_pth = Gkxmfs_pth + '/Cks' 
	SCk_pth = Gkxmfs_pth + '/SCk' 
	bs_pth = Gkxmfs_pth + '/Bs'

	# plotpth (local plot, not map plots)
	plotpth0 = self.flatfile_plot_pth00 + '/Gkxmfs'
	plotpth = plotpth0 + '/' + opt
	pfmt = 'eps'

	for f in [Gkxmfs_pth, metapth, \
	          bs_pth, dkxs_pth, SDks_pth, \
	          plotpth0, plotpth, ]:
	    if not os.path.exists(f):
		os.mkdir(f)

	Nk = len( CyberShakeRvar )

	# prepare flatfile  (just once!)
	D_flat = []; sites_flat = []
	for ik in xrange( Nk ):
	    sid = Sources[ik] 
	    sites_flat.append( self.sites_flatfile( sid, sites_info[ik], rups_info[ik] ) )
	    D_flat.append( self.D_flatfile( sid, sites_info[ik] ) )

	# top loop
	for it, Ti in enumerate( self.periods ): 
	    # loop over periods (it could be one period)
	    Tkey = '%.2f'%Ti
	    print 'Analysis at periods = %s sec'%(Tkey)

	    # source dependent (list)
	    Gkxmfs = []
	    pdf_f = []; pdf_m = []; Nhs = []
	    for ik in xrange( Nk ):

		Nm = len( CyberShakeRvar[ik] )
		tmpCS = []
		for irup in xrange( Nm ):
		    tmpCS.append( np.log( CyberShakeRvar[ik][irup][Tkey] ) )   # [irup][ih,islip,ista]
		tmp_CS = np.array( tmpCS )
		Nm,Nh,Nf,Nsta = tmp_CS.shape
		Nhs.append( Nh )

		Gkxmfs.append( tmp_CS )
		del(tmpCS)

		# set up weighting functions for slip-distribution and stress drop
		pdf_f.append([])
		pdf_m.append([])

		# Generate the distribution of each variabiles 
		Mws, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
		AreaTmp = Nrow*Ncol  # in km
		mu = 3.87 + np.log10(AreaTmp) * 1.05   # somerville 2006 M-A relationship
		std = 0.2

		prob0 = 1./(np.sqrt(2.*np.pi)*std) * np.exp(-0.5*((Mws-mu)/std)**2)
		pdf_m[ik] = prob0 / sum(prob0)    # weighting should be summed up to 1
		pdf_m[ik] = pdf_m[ik].tolist()

		# slip distribution (just use uniform distribution)
		for ivar in xrange( Nf ):
		    pdf_f[ik].append( 1./Nf ) 

	    # weighting functions for source and sites
	    pdf_k = np.ones(Nk)/Nk
	    pdf_s = np.ones(self.Nloc)/self.Nloc

	    # hypocenter location distributions:
	    keytmp = hypoPDF['CyberShake'].keys()[0]
	    pdf_x = []; Nsigma = 1
	    sigmas = hypoPDF['CyberShake'][keytmp]
	    try: 
		Nsigma = len(sigmas) 
	    except: 
		Nsigma = 1 
		sigmas = [sigmas,]

	    for isigma in xrange( Nsigma ):
		pdf_x.append([]) 
		for ik in xrange( Nk ):
		    Nh = Nhs[ik]   # different source has different hypocenter numbers
		    xh = np.arange( Nh ) / float(Nh-1)
		    mu = 0.5
		    #sigma = sigmas[isigma] / float(Nh-1)    # from observations
		    sigma = sigmas[isigma]   # input sigma should be in normalized unit (0,1), then it will be the same for each source!
		    prob0 = np.exp(-0.5*((xh-mu)/sigma)**2)   # normalized Gaussian distribution
		    pdf_x0 = (prob0/sum(prob0)).tolist()    # to make sure sum(prob) = 1
		    pdf_x[isigma].append( pdf_x0 )

	    # use NGA-derived parameters (default) (to prepare for later usage)
	    # model dependent (dict) or just two models (you can use list)
	    IDPk = []
	    ngaPk = []; ngaDk = [] 
	    for ik in xrange( Nk ):
		Nm, Nh, Nf, Nsta = Gkxmfs[ik].shape
		tmpG = np.zeros( (Nm,Nh,Nf,Nsta) )
		sid = Sources[ik] 

		if ngaP0k == None:
		    ngaP = self.nga_Py( sid, rups_info[ik], sites_info[ik], sites_flat=sites_flat[ik], Ts = [Ti,], model_name = NGAmodel )
		    if opt == 'Dkxs': 
			ngaPk.append( ngaP )

		if ngaD0k == None:
		    Mws = rups_info[ik][0]
		    ngaD, IDP = self.directivity_SC08(sid, Mws, Nh, sites_info[ik], D_flat=D_flat[ik], Ts = [Ti,], model_name = NGAmodel, IDPcpt=True)
		    if opt == 'Bs': 
			ngaDk.append( ngaD )
		    if opt == 'Dkxs':
			IDPk.append( np.array(IDP) )

		# sites
		meta = sites_info[ik]
		sites_run = meta.keys()
		slon = []; slat = []
		for i in xrange( Nsta ):
		    slon.append( meta[sites_run[i]]['lon'] )
		    slat.append( meta[sites_run[i]]['lat'] )
		slon = np.array( slon )
		slat = np.array( slat )

		if ik == 0:
		    # interpolate the site info (Vs30, Z2.5, Z1.0)
		    sites_info0 = self.sites_flatfile( Sources[ik], sites_info[ik] )    # Vs30 and Z1.0, Z2.5 (uninterpolated results)
		    Vs30 = interp( slon, slat, sites_info0[:,1], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    Z25 = interp( slon, slat, sites_info0[:,2], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    Z10 = interp( slon, slat, sites_info0[:,3], self.slon1d, self.slat1d, eps=self.eps, method=self.method )

	    if ngaP0k != None: 
		ngaPk = ngaP0k[nga]
	    if ngaD0k != None: 
		ngaDk = ngaD0k[nga] 

	    # =========================================================
	    # the above will be all the same for different opts (Bs,Dkxs, or Cks)
	    # =========================================================
	    if opt == 'Bs': 
		NGAdict = self.NGAs
		NGAdict[NGAmodel]['NewCoefs'] = {}   # the default of NewCoefs == None

		# Generate new parameters
		if NGAmodel == 'AS': 
		    # generate grid search values around NGA values
		    ASnga = AS08_nga()
		    Nside = 5
		    dprm = 0.2
		    coef0 = ASnga.Coefs[Tkey]['a10']
		    coef00 = coef0
		    dtmp = abs(coef0 * dprm)
		    coefs = np.arange( -Nside*dtmp, Nside*dtmp+dtmp, dtmp ) + coef0

		# Loop over possible parameters (given by the opt dictionary)
		sigmaB = []; sigmaBmin = 10000 
		for item in xrange( len(coefs) ):

		    # use the updated coefficients:
		    NGAdict[NGAmodel]['NewCoefs'][coefkey] = coefs[item]

		    # recompute ngaPy and ngaD
		    Hkxs = []
		    for ik in xrange( Nk ):
			Nm, Nh, Nf, Nsta = Gkxmfs[ik].shape
			tmpG = np.zeros( (Nm,Nh,Nf,Nsta) )
			sid = Sources[ik] 
			ngaP = self.nga_Py( sid, rups_info[ik], sites_info[ik], Ts = [Ti,], model_name = NGAmodel,NGAdict=NGAdict )

			for ir in xrange( Nm ):
			    for ista in xrange( Nsta ):
				for ih in xrange( Nh ): 
				    tmpG[ir,ih,:,ista] = Gkxmfs[ik][ir,ih,:,ista] - np.log( ngaP[Tkey][ir][0][ista] ) - np.log( np.array(ngaDk[ik][Tkey][ir])[ih,ista] )
			tmp_xms = np.average( tmpG, axis=2, weights = pdf_f[ik] )
			Hkxs0 = np.average( tmp_xms, axis=0, weights=pdf_m[ik] )

			# sites
			meta = sites_info[ik]
			sites_run = meta.keys()
			slon = []; slat = []
			for i in xrange( Nsta ):
			    slon.append( meta[sites_run[i]]['lon'] )
			    slat.append( meta[sites_run[i]]['lat'] )
			slon = np.array( slon )
			slat = np.array( slat )

			tmp = np.zeros((Nh,self.Nloc))
			for ih in xrange( Nh ):
			    tmp[ih,:] = interp( slon, slat, Hkxs0[ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			Hkxs.append( tmp ) 

		    del( tmpG )

		    # we have TargetModel as (residual) that has shape: [Nk][Nm,Nh,Nf,Nsta]
		    # Decomposition of Hkxs
		    dkxs, Sdks, cks, cdks, bs, a = im_decom2( Hkxs, self.Nloc, pdf_x[0], pdf_k, pdf_s )
		    sigmaB0 = np.average( bs**2, axis=0, weights=pdf_s )   # minimize the variance (not the standard deviation)
		    if sigmaB0 <= sigmaBmin: 
			sigmaBmin = sigmaB0
			coef0 = coefs[item]
			bs0 = bs
		    sigmaB.append( sigmaB0 )

		# plot different parameters and sigmaB
		fig = plt.figure(1) 
		fig.clf()
		ax = fig.add_subplot(111) 
		ax.plot( coefs, sigmaB, 'ro' )
		ax.set_xlabel( '%s'%coefkey )
		ax.set_ylim([0,1])
		ax.set_ylabel( r'$\sigma_b$' )
		ax.set_title( 'NGA model: %s, period: %s sec, original a10 = %s '%(NGAmodel, '%.2f'%Ti, coef00) )
		ax.text( coef0, max(sigmaB)*0.5, '(%s,%s)'%('%.4f'%coef0,'%.4f'%sigmaBmin) )   # mark the minimum value 

		fig.savefig( plotpth + '/Model%s.SigmaB_%s_%s.pdf'%(NGAmodel, '%.2f'%Ti, coefkey), format='pdf')
		fig.savefig( plotpth + '/Model%s.SigmaB_%s_%s.%s'%(NGAmodel, '%.2f'%Ti, coefkey, pfmt), format=pfmt, dpi=300 )

		# find the smallest bs and write into Basin maps and overall averaging values
		bsFile = bs_pth + '/CyberShake.NGAs.%s.bs'%(Tkey)
		fid_bs = open( bsFile, 'w')
		for iloc in xrange( self.Nloc ):
		    fid_bs.write( '%s %s %s %s %s %s\n'%(\
		        self.slon1d[iloc], self.slat1d[iloc], \
		        Vs30[iloc], Z10[iloc], Z25[iloc], \
		        bs0[iloc]))    # write the minimun sigmaB for basin plots
		fid_bs.close()

		# write into metafile 
		metafile = metapth + '/CyberShake.%s.%s.%s.T%s.py'%(NGAmodel, opt, coefkey, '%.2f'%Ti) 
		meta = dict( 
		    coefs = coefs.tolist(), 
		    SigmaS = sigmaB, 
		)
		save( metafile, meta, header='#metafile for minimize %s by grid searching coeficient %s in %s\n'%(opt, coefkey, NGAmodel) )

		return coef0 

	    if opt == 'Dkxs': 
		# ==================================================
		# find corresponding parameters to miniminze sigmaD 
		# ==================================================
		dict1 = self.NGAs
		dict1['SC']['NewCoefs'][NGAmodel] = {}   # the default of NewCoefs == None
		cuteps = dict1['SC']['predictors']['cuteps']

		# generate grid search values around NGA values
		SC08 = SC08_model(NGAmodel+'08',cuteps = cuteps)
		Nside = 5
		coef0 = SC08.Coefs[Tkey][coefkey]
		dtmp = 0.2
		coefs = np.arange( -Nside*dtmp, Nside*dtmp+dtmp, dtmp ) + coef0

		coef00 = coef0   # original b

		# Loop over possible parameters (given by the opt dictionary)
		dkxsDmin = []; SdksDmin = []
		sigmaDmin = []
		sigmaS_D = {}
		fig = plt.figure(1) 
		sigmaStr = []; lines = []
		for isigma in xrange( Nsigma ):
		    #print 'deal with half width %s'%(sigmas[isigma])
		    sigmaDkey = '%s'%'%.2f'%sigmas[isigma]
		    sigmaD = []; sigmaD0min = 10000 
		    for item in xrange( len(coefs) ):


			#print 'deal with coef %s'%item 

			# use the updated coefficients:
			#dict1['SC']['NewCoefs'][NGAmodel][coefkey] = coefs[item]

			a = SC08.Coefs[Tkey]['a']    # keep other parameters (period dependent)
			b = coefs[item]    # b is changed
			a0 = SC08.Coefs[Tkey]['a0']    

			# recompute ngaD
			dkxs = []; Sdk = []; Sdks0 = []   # uninterpolated
			Hkxs = []     # interpolated
			for ik in xrange( Nk ):
			    Nm, Nh, Nf, Nsta = Gkxmfs[ik].shape
			    tmpG = np.zeros( (Nm,Nh,Nf,Nsta) )
			    sid = Sources[ik] 
			    Mws = rups_info[ik][0]
			    #print 'compute directivity for source %s'%sid

			    # computing time consuming:
			    #ngaD = self.directivity_SC08(sid, Mws, Nh, sites_info[ik], D_flat=D_flat[ik], Ts = [Ti,], model_name = NGAmodel, NGAdict = dict1)
			    for ir in xrange( Nm ):
				for ista in xrange( Nsta ):
				    for ih in xrange( Nh ): 
					#tmpG[ir,ih,:,ista] = Gkxmfs[ik][ir,ih,:,ista] - np.log( ngaPk[ik][Tkey][ir][0][ista] ) - np.log( np.array(ngaD[Tkey][ir])[ih,ista] )
					tmpG[ir,ih,:,ista] = Gkxmfs[ik][ir,ih,:,ista] - np.log( ngaPk[ik][Tkey][ir][0][ista] ) - (a+b*IDPk[ik][ih,ista]+a0)

			    tmp_xms = np.average( tmpG, axis=2, weights = pdf_f[ik] )
			    Hkxs0 = np.average( tmp_xms, axis=0, weights=pdf_m[ik] )

			    if 1:
				# sites for interpolation
				meta = sites_info[ik]
				sites_run = meta.keys()
				slon = []; slat = []
				for i in xrange( Nsta ):
				    slon.append( meta[sites_run[i]]['lon'] )
				    slat.append( meta[sites_run[i]]['lat'] )
				slon = np.array( slon )
				slat = np.array( slat )

				tmp = np.zeros((Nh,self.Nloc))
				for ih in xrange( Nh ):
				    tmp[ih,:] = interp( slon, slat, Hkxs0[ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
				Hkxs.append( tmp ) 
			    else:
				# compute dkxs (un-interpolated):
				# this may cause the problem that when you do the interpolation, the average of dkxs over x is not zeros.
				tmp = np.zeros( (Nh,Nsta) )
				for ih in xrange(Nh):
				    tmp[ih,:] = Hkxs0[ih,:] - np.average( Hkxs0, axis=0, weights = pdf_x[isigma][ik] )
				dkxs.append( tmp )
				Sdks00 = np.average( tmp**2, axis=0, weights=pdf_x[isigma][ik] ) 
				Sdks0.append( Sdks00 )

				Sdk0 = np.average( Sdks00, axis=0, weights=np.ones(Nsta)/Nsta )
				Sdk.append( Sdk0 )

			del( tmpG )
			if 1:
			    dkxs, Sdks, cks, cdks, bs, a000 = im_decom2( Hkxs, self.Nloc, pdf_x[isigma], pdf_k, pdf_s )
			    Sdks = np.array( Sdks, 'f' )   # standard deviation
			    sigmaD0 = np.average( np.average( Sdks**2, axis=0, weights = pdf_k ), axis=0, weights=pdf_s ) 
			else:
			    sigmaD0 = sum(Sdk) / Nk   # average over all Sources

			if sigmaD0 <= sigmaD0min: 
			    sigmaD0min = sigmaD0
			    coef0 = coefs[item]
			    dkxs0 = dkxs 
			    if 1:
				Sdks2 = Sdks
			sigmaD.append( sigmaD0 )

		    sigmaDmin.append( sigmaD0min )  
		    sigmaS_D[sigmaDkey] = sigmaD 

		    if 1:
			dkxsDmin.append( dkxs0 )
			SdksDmin.append( Sdks2 )
		    else: 
			# interpolate dkxs0 
			dkxs1 = []; Sdks1 = []
			for ik in xrange( Nk ):
			    meta = sites_info[ik]
			    sites_run = meta.keys()
			    slon = []; slat = []
			    for i in xrange( Nsta ):
				slon.append( meta[sites_run[i]]['lon'] )
				slat.append( meta[sites_run[i]]['lat'] )
			    slon = np.array( slon )
			    slat = np.array( slat )

			    tmp = np.zeros((Nh,self.Nloc))
			    for ih in xrange( Nh ):
				tmp[ih,:] = interp( slon, slat, dkxs0[ik][ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			    dkxs1.append( tmp ) 

			    # you can compute Sdks directly from the interpolated dkxs or you can interpolate Sdks0 !
			    tmp1 = np.average( tmp**2, axis=0, weights = pdf_x[isigma][ik] )
			    #tmp = interp( slon, slat, np.sqrt(Sdks0[ik]), self.slon1d, self.slat1d, eps=self.eps, method=self.method )
			    Sdks1.append( np.sqrt(tmp1) ) 

			dkxsDmin.append( dkxs1 ) 
			SdksDmin.append( Sdks1 )

		    # plot different parameters and sigmaD
		    ax = fig.add_subplot(111) 
		    line0 = ax.plot( coefs, sigmaD, '^-')
		    lines.append( line0 )
		    sigmaStr.append( r'$\gamma = %.3f$'%sigmas[isigma] )
		    ax.plot( coef0, sigmaD0min, 'k*' )

		lg = ax.legend( lines, sigmaStr, loc=0 )
		lg.draw_frame(False) 

		ax.set_xlabel( '%s'%coefkey )
		#ax.set_ylim([0,2])
		ax.set_ylabel( r'$\sigma_d$' )
		ax.set_title( 'NGA model: %s, period: %s sec, original b = %s'%(NGAmodel, '%.2f'%Ti, '%.3f'%coef00) )
		fig.savefig( plotpth + '/Model%s.SigmaDs_%s_%s.pdf'%(NGAmodel, '%.2f'%Ti, coefkey), format='pdf')
		fig.savefig( plotpth + '/Model%s.SigmaDs_%s_%s.%s'%(NGAmodel,   '%.2f'%Ti, coefkey, pfmt), format=pfmt, dpi=300 )

		# find the minimum sigmaD for all width 
		index = (np.array(sigmaDmin) == min( sigmaDmin )).nonzero()[0]   # find the index and corresponding half width
		sigmaGlobal = sigmas[index] 

		# D-map (write)
		srcsh_fid = open( Gkxmfs_pth + '/SourceRuptureHypoInfo','w' )
		for ik in xrange( Nk ):
		    sid = Sources[ik]
		    Nh = dkxs0[ik].shape[0]

		    Mws, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
		    Nm = len(Mws) 

		    # source info
		    srcsh_fid.write('%s %s %s\n'%(sid,Nm,Nh))

		    dkxs_pth0 = dkxs_pth + '/%s'%sid
		    for f in [dkxs_pth0,]:
			if not os.path.exists(f):
			    os.mkdir(f)
		    for ih in xrange(Nh):
			dkxsFile = dkxs_pth0 + '/CyberShake.%s.%s.Source%s.Ih%s.dkxs'%(NGAmodel,Tkey, sid, ih)
			fid = open( dkxsFile, 'w' )
			for iloc in xrange(self.Nloc): 
			    str0 = ''; str1 = ''
			    for isigma in xrange( Nsigma ):
				pdf_x0 = pdf_x[isigma][ik][ih]/max(pdf_x[isigma][ik])
				str0 = str0 + ' %s'%(dkxsDmin[isigma][ik][ih,iloc])
				str1 = str1 + ' %s'%(pdf_x0 * dkxsDmin[isigma][ik][ih,iloc])
			    fid.write( '%s %s%s%s\n'%(\
			        self.slon1d[iloc], self.slat1d[iloc],\
			        str0, str1 ))
			fid.close()

		    # Variance of D-Map
		    SDksFile = SDks_pth + '/CyberShake.%s.%s.Source%s.Sdks'%(NGAmodel, Tkey, sid)
		    fid = open( SDksFile, 'w' )
		    for iloc in xrange(self.Nloc): 
			str0 = ''
			for isigma in xrange( Nsigma ):
			    str0 = str0 + ' %s'%(SdksDmin[isigma][ik][iloc])
			fid.write( '%s %s%s\n'%(\
			    self.slon1d[iloc], self.slat1d[iloc],str0) )
		    fid.close()

		srcsh_fid.close() 

		# write into metafile 
		metafile = metapth + '/CyberShake.%s.%s.%s.T%s.py'%(NGAmodel, opt, coefkey, '%.2f'%Ti) 
		meta = dict( 
		    coefs = coefs.tolist(), 
		    SigmaS = sigmaS_D,
		)
		save( metafile, meta, header='# metafile for minimize %s by grid searching coeficient %s in %s\n'%(opt, coefkey, NGAmodel) )

		return coef0, sigmaGlobal



    #  After you do the Gkxmfs Analysis, you can use GMT to plot map-based factors, and then you can you this one to plot A and F_k
    def PlotAkE0(self, sigma, Ref = 'CyberShake'):
	# read Bs, Ak and E0 files
	# Bs like plots in R
	# refer to gmt way for different reference models
	sigmakey = 'Sigma%s'%sigma 

	# reference model consideration
	if Ref == 'CyberShake':
	    iref = 4
	elif Ref == 'NoRef':
	    iref = -1
	elif Ref == 'Ref':
	    iref = 5
	else:
	    for inga in xrange( len(self.NGAmodel) ):
		if Ref == self.NGAmodel[inga]:
		    iref = inga
		    break

	Np = len(self.periods) 
	if Np == 1:
	    Tstr='%.2f'%self.periods[0]
	else:
	    Tstr='periods'

	#Gkxmfs_pth=self.output_pth+'/Gkxmfs0'      # path to save the interpolated results
	Gkxmfs_pth=self.output_pth+'/Gkxmfs'      # path to save the interpolated results
	cks_pth = Gkxmfs_pth + '/Cks/'+sigmakey
	a_pth = Gkxmfs_pth + '/A/'+sigmakey

	src_info = Gkxmfs_pth + '/SourceInfo'
	src_info = np.loadtxt( src_info ) 

	# loop over all periods
	e0 = []; ak = []
	for ip in xrange( Np ):
	    period = self.periods[ip] 
	    Tkey = '%.2f'%period 
	    file_e0 = a_pth + '/CyberShake.NGAs.%s.a'%('%.2f'%period)
	    #tmp_e0 = np.loadtxt( file_e0, skiprows=1 )
	    tmp_e0 = np.loadtxt( file_e0 )
	    tmp_ak = []
	    for ik in xrange( len(self.sids) ):
		sid = self.sids[ik]
		cksFile = cks_pth + '/CyberShake.NGAs.%s.Source%s.cks'%(Tkey, sid)
		tmp = np.loadtxt( cksFile,usecols=(2,3,4,5,6,7) ) 
		tmp_ak.append( np.r_[src_info[ik], np.mean(tmp,0)] ) 
	    tmp_ak = np.array( tmp_ak )

	    try: 
		n1,n2 = tmp_ak.shape
	    except: 
		tmp_ak = tmp_ak.reshape( 1, len(tmp_ak) )

	    tmp_e0[:-2] = tmp_e0[-2] - tmp_e0[:-2]   # NGA models
	    tmp_e0[-1] = tmp_e0[-2] - tmp_e0[-1]   # Ref model 
	    for ik in xrange( len(tmp_ak) ):
		tmp_ak[ik,2:-2] = tmp_ak[ik,-2] - tmp_ak[ik,2:-2]  # NGA models
		tmp_ak[ik,-1] = tmp_ak[ik,-2] - tmp_ak[ik,-1]      # Ref model 
	    if iref not in [-1,4]:
		tmp_e0[:-1] = tmp_e0[:-1] - tmp_e0[iref]    # including CyberShake as the target model
		for ik in xrange( len( tmp_ak ) ):
		    tmp_ak[ik,2:-1] = tmp_ak[ik,2:-1]-tmp_ak[ik,iref+2]
	    if iref == 4: 
		tmp_e0[:-1] = tmp_e0[iref] - tmp_e0[:-1]    # including CyberShake as the target model
		for ik in xrange( len( tmp_ak ) ):
		    tmp_ak[ik,2:-1] = tmp_ak[ik,iref+2] - tmp_ak[ik,2:-1]

	    e0.append( tmp_e0 )
	    ak.append( tmp_ak )

	e0 = np.array( e0 )
	ak = np.array( ak )

	# plotpth
	plotpth0 = self.flatfile_plot_pth00 + '/Gkxmfs'
	plotpth = plotpth0 + '/Ak'
	plotpth1 = plotpth0 + '/A'
	for f in [plotpth0, plotpth, plotpth1]:
	    if not os.path.exists( f ):
		os.mkdir( f )

	plotpth0 = plotpth0 + '/' + sigmakey
	plotpth = plotpth + '/' + sigmakey
	plotpth1 = plotpth1 + '/' + sigmakey
	for f in [plotpth0, plotpth, plotpth1]:
	    if not os.path.exists( f ):
		os.mkdir( f )

	pfmt = 'eps'

	# e0
	fig1 = plt.figure(1,(10,8),dpi=600) 
	plt.rc('font',family='Arial')

	ax = fig1.add_subplot( 111 )
	ax1 = ax.twinx()

	clr = ['k','b','r','g','m']
	smb = ['o','^','s','v','*']
	lgnames = self.NGAmodel + ['CyberShake']
	lines = []
	for imodel in [3,1,0,2,4]:
	#for imodel in xrange( 5 ):
	    if imodel != 4:
		label0 = lgnames[imodel]+'08'
	    else: 
		label0 = lgnames[imodel]
	    print lgnames[imodel], e0[:,imodel]
	    if imodel == 4:
		ax.loglog( self.periods, np.exp(e0[:,imodel]), 'k-'+smb[imodel], ms=12, label=label0 )
	    else: 
		ax.loglog( self.periods, np.exp(e0[:,imodel]), 'k-'+smb[imodel], mfc='none', ms=8, label=label0 )
	    if imodel == 4:
		ax1.loglog( self.periods, np.exp(e0[:,imodel]), 'k-'+smb[imodel], ms=12, label=label0,basey=np.exp(1) )
	    else: 
		ax1.loglog( self.periods, np.exp(e0[:,imodel]), 'k-'+smb[imodel], mfc='none', ms=8, label=label0, basey=np.exp(1) )

	    if imodel==4 and 0:
		# linear fit 
		from scipy import polyfit, polyval
		npoly = 1
		w = polyfit(np.log10(self.periods),np.log10(np.exp(e0[:,imodel])),npoly)
		print w
		#line = w[0]*np.log10(self.periods)+w[1] 
		line = polyval( w, np.log10(self.periods) )
		ax.loglog( self.periods, 10**line, 'r--',lw=2 ) 

	lg = ax.legend( loc=3, numpoints=1 )
	#lg.draw_frame( False )
	ax.set_ylim([0.005,0.1])
	ax1.set_ylim([0.005,0.1])
	ax.set_xlim([1.8, 11.])
	ax.grid(True)
	#ax.set_ylabel('ln A')
	from matplotlib.ticker import LogLocator, FormatStrFormatter 
	minorLocator = LogLocator(subs = np.linspace(2,10,8,endpoint=False))
	ax.xaxis.set_minor_locator( minorLocator) 
	ax.grid(b=True,which='minor') 
	ax.yaxis.set_ticks_position('right')

	if 1:
	    subs = [-3, -4, -5]
	    majorLocator = LogLocator(subs = np.exp(subs) )
	    majorFormatter = FormatStrFormatter('$e^{%.2f}$')
	    ax1.yaxis.set_major_locator( majorLocator )
	    ax1.yaxis.set_major_formatter(majorFormatter) 
	    ax1.yaxis.set_ticks_position('left')

	ax.set_xlabel('Period (s)')
	#if Ref == 'NoRef':
	#    ax1.set_ylabel('A (g)' )
	#else: 
	#    ax1.set_ylabel('$a$')

	# save fig 
	fig1.savefig( plotpth1+ '/E0_ref%s.%s'%(Ref,pfmt), format=pfmt,dpi=600)

	# ak
	fig2 = plt.figure(2,(14,10))
	xt = np.arange( ak.shape[1] )
	Nk = len(xt)

	# sort sid by the Area
	Area0, Sid0, index = sorti( ak[0,:,1], list1=ak[0,:,0] )
	xtlab = np.array(Sid0,'i4') 

	fig2.clf()
	for ip in xrange( Np ):
	    ax = fig2.add_subplot( Np, 1, ip+1 )
	    for imodel in xrange( 5 ):
		Area0, ak0, index = sorti( ak[0,:,1], list1=ak[ip,:,imodel+2] )
		ax.plot( Area0, ak0, color=clr[imodel], marker=smb[imodel], label=lgnames[imodel] )

	    ax.set_xlim([0, 15000])
	    ax.set_ylim([-1,1])
	    ax.text( 10, 0.8, 'T=%s'%('%.2f'%self.periods[ip]), ha='left', va='top' )

	    if ip == 0:
		lg = ax.legend( loc = 0 )
		lg.draw_frame( False )
	    if ip == Np-1:
		tmpa = 0.15
		for i in xrange( len(Area0) ):
		    ax.text( Area0[i], -0.9+tmpa*np.mod(i,2), '%s'%'%4i'%Sid0[i], ha='left', va='bottom' )
		ax.set_xlabel( 'Rupture Area (km^2)' )
	    else:
		ax.set_xticks([])
	    ax.set_ylabel( 'fk (ln)' )

	# save fig 
	fig2.savefig( plotpth + '/Ak_%s.ref%s.%s'%(Tstr, Ref,pfmt), format=pfmt)


    def BasinDepthCorr(self, Zx0, sigma, PlotAll=False):

	sigmakey = 'Sigma%.2f'%sigma
	Gkxmfs_pth=self.output_pth+'/Gkxmfs'      # path to save the interpolated results
	bs_pth = Gkxmfs_pth + '/Bs/' + sigmakey 

	plotpth0 = self.flatfile_plot_pth00 + '/Gkxmfs'
	plotpth1 = plotpth0 + '/Bs'
	plotpth = plotpth1 + '/' + sigmakey
	for f in [plotpth0, plotpth1, plotpth,]:
	    if not os.path.exists( f ):
		os.mkdir( f )
	pfmt = 'eps'

	Np = len(self.periods) 
	if Np == 1:
	    Tstr='%.2f'%self.periods[0]
	else:
	    Tstr='periods'

	for ip in xrange( Np ):
	    period = self.periods[ip] 
	    Tkey = '%.2f'%period 
	    Binput = np.loadtxt( bs_pth + '/CyberShake.NGAs.%s.bs'%('%.2f'%period) )

	    # these could be read in from velocity models (...)
	    Vs30 = Binput[:,2]
	    Z10 = Binput[:,3]
	    Z25 = Binput[:,4] 
	    ZX = [Z10,Z25]
	    xtlab = ['Z1.0','Z2.5']
	    xt = np.arange( len(xtlab) )

	    Bs = Binput[:,-2]  # CyberShake Bs

	    # residual analysis
	    NGAbs = Binput[:,5:-2]

	    fig = plt.figure(1)  
	    for inga in xrange( 4 ):
		bs = NGAbs[:,inga] 
		ax = fig.add_subplot( 2,2,inga+1 )
		rho = []
		for iz in xrange( len(Zx0) ):
		    Zx = ZX[iz] 
		    Ztmp, bstmp, index = sorti( Zx, list1=bs )
		    index = (np.array(Ztmp)>=Zx0[iz]).nonzero()[0]
		    bstmp = np.array(bstmp)[index]
		    Ztmp = np.array(Ztmp)[index]
		    rho.append(np.corrcoef( bstmp, Ztmp )[0,1])
		print inga, rho 

		ax.plot( xt, rho, 'ro' )
		ax.set_xlim( [-1,2] )
		ax.set_ylim( [-1,1] )
		ax.set_xticks(xt)
		ax.set_title( '%s'%self.NGAmodel[inga] )
		ax.set_ylabel( r'$\rho$' )
		xticknames = plt.setp(ax,xticklabels=xtlab)
		plt.setp(xticknames, rotation=0, fontsize=10)

	    plotnam = '/BasinDepthCorrBs.NGAs.CyberShake.%s.%s'%(Tkey,pfmt)
	    fig.savefig( plotpth + plotnam, format = pfmt,dpi=300 )

	    if PlotAll:
		rs_Z10 = []; rs_Z25 = []
		fig.clf()
		for inga in xrange( 4 ):
		    bs = NGAbs[:,inga] 
		    rs_Z10.append([])
		    rs_Z25.append([])
		    for ista in xrange( len(Z10) ):
			Zx0 = [Z10[ista],Z25[ista]]  # loop over each criterion
			rho = []
			for iz in xrange( len(Zx0) ):
			    Zx = ZX[iz] 
			    Ztmp, bstmp, index = sorti( Zx, list1=bs )
			    index = (np.array(Ztmp)>=Zx0[iz]).nonzero()[0]
			    bstmp = np.array(bstmp)[index]
			    Ztmp = np.array(Ztmp)[index]
			    rho.append(np.corrcoef( bstmp, Ztmp )[0,1])

			rs_Z10[inga].append( rho[0] )
			rs_Z25[inga].append( rho[1] )

		    ax = fig.add_subplot( 2,2,inga+1 )
		    ax.plot( Z10, rs_Z10[inga], 'rx', Z25*1000, rs_Z25[inga], 'b+' )
		    ax.set_title(self.NGAmodel[inga])
		    ax.set_ylim([-1,1])
		    ax.set_ylabel( r'$\rho$' )
		    if inga in [2,3]:
			ax.set_xlabel( '>=BasinDepth (m)')
		    lg = ax.legend( ('Z1.0','Z2.5'), loc=0 )
		    lg.draw_frame(False)
		fig.savefig( plotpth + '/BasinDepthCorrBsAll.NGAs.CyberShake.%s.%s'%(Tkey,pfmt) )

		bsFile = bs_pth + '/CyberShake.NGAs.%s.rsZ10'%(Tkey)
		fid_bs = open( bsFile, 'w')
		bsFile1 = bs_pth + '/CyberShake.NGAs.%s.rsZ25'%(Tkey)
		fid_bs1 = open( bsFile1, 'w')
		for iloc in xrange( self.Nloc ):
		    fid_bs.write( '%s %s %s %s %s %s %s %s %s %s\n'%(\
		        self.slon1d[iloc], self.slat1d[iloc], \
		        Vs30[iloc], Z10[iloc], Z25[iloc], \
		        rs_Z10[0][iloc],rs_Z10[1][iloc],rs_Z10[2][iloc], rs_Z10[3][iloc], Bs[iloc] ) )
		    fid_bs1.write( '%s %s %s %s %s %s %s %s %s %s\n'%(\
		        self.slon1d[iloc], self.slat1d[iloc], \
		        Vs30[iloc], Z10[iloc], Z25[iloc], \
		        rs_Z25[0][iloc],rs_Z25[1][iloc],rs_Z25[2][iloc], rs_Z25[3][iloc], Bs[iloc] ) )
		fid_bs.close()
		fid_bs1.close()



    def RMSanalysis(self, sigmas, refmodel='BA', opt='Dkxs', wt=0): 
	""" 
	RMS analysis of each factors (d,c,b)
	"""

	# paths
	Gkxmfs_pth= self.output_pth+'/Gkxmfs'      # path to save directivity analysis results (weighting functions)

	# map file save pth
	dkxs_pth = Gkxmfs_pth + '/Dkxs' 
	cks_pth = Gkxmfs_pth + '/Cks' 
	bs_pth = Gkxmfs_pth + '/Bs'

	plotpth0 = self.flatfile_plot_pth00 + '/Gkxmfs'
	plotpth = plotpth0 + '/' + opt
	pfmt = 'eps'

	for f in [Gkxmfs_pth, \
	          cks_pth, dkxs_pth, bs_pth, \
	          plotpth0, plotpth, ]:
	    if not os.path.exists(f):
		os.mkdir(f)

	if wt: 
	    weighted='weighted'
	else: 
	    weighted='unweighted'

	# Reference model
	iref = 5   # no reference  
	for i in xrange( len(self.NGAmodel) ):
	    if refmodel == self.NGAmodel[i]: 
		iref = i
		break 
	if refmodel == 'CyberShake': 
	    iref = 4 

	try: 
	    Nsigma = len(sigmas) 
	except: 
	    sigams = [sigmas,]
	    Nsigma = len(sigmas) 

	TargetModels = self.NGAmodel + ['CyberShake']

	# weighting function for directivity
	# read in file
	srcsh_file = Gkxmfs_pth + '/SourceRuptureHypoInfo'
	srcNmh = np.loadtxt( srcsh_file, dtype={'names':('sids','Nms','Nhs'),'formats':('i4','i4','i4')} )  

	for it in xrange( len(self.periods) ):
	    Ti = self.periods[it]
	    Tkey = '%.2f'%Ti 

	    fig = plt.figure(1) 
	    lines = []; sigmaStr = []
	    for isigma in xrange( Nsigma ): 
		sigmakey = '%.2f'%sigmas[isigma] 
		SigmaS = '/Sigma%s'%(sigmakey)
		dkxs_pth1 = dkxs_pth + SigmaS
		cks_pth1  = cks_pth  + SigmaS 
		bs_pth1  = bs_pth  + SigmaS 

		Nk = len(srcNmh['sids'])
		Nhs = srcNmh['Nhs']

		if opt == 'Dkxs':
		    Ddata = []
		    for ik in xrange( Nk ):
			sid = srcNmh['sids'][ik]
			Nh = Nhs[ik]
			Ddata.append([])
			dkxs_pth0 = dkxs_pth1 + '/%s'%sid
			for ih in xrange(Nh):
			    dkxsFile = dkxs_pth0 + '/CyberShake.NGAs.%s.Source%s.Ih%s.dkxs'%(Tkey, sid, ih)
			    Alldata = np.loadtxt( dkxsFile )[:,2:]
			    Ddata0 = np.r_[ (Alldata[:,4+wt*8:5+wt*8] - Alldata[:, wt*8:4+wt*8]).T, Alldata[:,4+wt*8:5+8*wt].T ].T  # get capital D
			    if iref != 5: 
				Ddata0 = Ddata0 - Ddata0[:,iref:iref+1]    # residual
			    Ddata[ik].append( Ddata0 )

		    # Ddata has the dimension: [ik][ih][ista,imodel]
		    RMS = []   # five elements
		    for imodel in xrange( len(TargetModels) ):
			RMS_k = 0
			for ik in xrange( Nk ):
			    Nh = Nhs[ik]
			    RMS_h = 0
			    for ih in xrange( Nh ):
				RMS_h += sum( Ddata[ik][ih][:,imodel]**2 ) / self.Nloc 
			    RMS_k += RMS_h / Nh
			RMS.append( np.sqrt(RMS_k/Nk) )
		if opt == 'Cks':
		    Cdata = []
		    for ik in xrange( Nk ):
			sid = srcNmh['sids'][ik]
			cksFile = cks_pth1 + '/CyberShake.NGAs.%s.Source%s.cks'%(Tkey, sid)
			Alldata = np.loadtxt( cksFile )[:,2:]
			Cdata0 = np.r_[(Alldata[:,4:5]-Alldata[:,:4]).T, Alldata[:,4:5].T].T
			if iref != 5: 
			    Cdata0 = Cdata0 - Cdata0[:,iref:iref+1] 
			Cdata.append( Cdata0 )
		    RMS = []
		    for imodel in xrange( len(TargetModels) ):
			RMS_k = 0
			for ik in xrange( Nk ):
			    RMS_k += sum(Cdata[ik][:,imodel]**2) / self.Nloc 
			RMS.append( np.sqrt(RMS_k/Nk) )

		if opt == 'Bs':
		    bsFile = bs_pth1 + '/CyberShake.NGAs.%s.bs'%Tkey
		    Alldata = np.loadtxt( bsFile )[:,5:] 
		    Bdata = np.r_[(Alldata[:,4:5]-Alldata[:,:4]).T, Alldata[:,4:5].T].T
		    if iref != 5: 
			Bdata = Bdata - Bdata[:,iref:iref+1] 
		    RMS = []
		    for imodel in xrange( len(TargetModels) ):
			RMS.append( np.sqrt(sum(Bdata[:,imodel]**2)/self.Nloc) )

		# plot 
		#fig.clf()
		ax = fig.add_subplot( 111 )
		xt = np.arange( len( TargetModels ) )
		line0 = ax.plot( xt, RMS, 'o-' )
		lines.append( line0 ) 
		sigmaStr.append( r'$\gamma = %.3f$'%sigmas[isigma] )

	    lg = ax.legend( lines, sigmaStr, loc=0 )
	    lg.draw_frame(False) 
	    ax.set_title( 'Refmodel = %s, period = %s'%(refmodel, Tkey) )
	    ax.set_xlim([-1,max(xt)+1])
	    ax.set_ylim([-0.5,1.0])
	    ax.set_ylabel( '%s RMS'%opt )
	    ax.set_xticks(xt)
	    xticknames = plt.setp(ax,xticklabels=TargetModels)
	    plt.setp(xticknames, rotation=0, fontsize=10)
	    if opt == 'Dkxs':
		fig.savefig( plotpth + '/CyberShake.NGAs.Sigmas.%s.Ref%s.RMS.%s.%s.%s'%(Tkey, refmodel,opt,weighted, pfmt), format=pfmt )
		fig.savefig( plotpth + '/CyberShake.NGAs.Sigmas.%s.Ref%s.RMS.%s.%s.pdf'%(Tkey,refmodel,opt, weighted), format='pdf')
	    else: 
		fig.savefig( plotpth + '/CyberShake.NGAs.Sigmas.%s.Ref%s.RMS.%s.%s'%(Tkey, refmodel,opt,pfmt), format=pfmt )
		fig.savefig( plotpth + '/CyberShake.NGAs.Sigmas.%s.Ref%s.RMS.%s.pdf'%(Tkey,refmodel,opt), format='pdf')




    def BarPlot(self,sigma,refmodel='NoRef',wt=0): 

	# Plot A value, B,C,D map RMS values
	sigmakey = 'Sigma%.2f'%sigma
	Gkxmfs_pth=self.output_pth+'/Gkxmfs'      # path to save the interpolated results
	dkxs_pth1 = Gkxmfs_pth + '/Dkxs/'+sigmakey 
	cks_pth1 = Gkxmfs_pth + '/Cks/'+sigmakey
	bs_pth1 = Gkxmfs_pth + '/Bs/'+sigmakey
	a_pth1 = Gkxmfs_pth + '/A/'+sigmakey

	RMSpth = Gkxmfs_pth + '/RMS' 

	plotpth0 = self.flatfile_plot_pth00 + '/Gkxmfs'
	plotpth1 = plotpth0 + '/All'
	plotpth = plotpth1 + '/' + sigmakey
	pfmt = 'eps'

	for f in [ RMSpth, plotpth0, plotpth1, plotpth, ]:
	    if not os.path.exists(f):
		os.mkdir(f)

	if wt: 
	    weighted='weighted'
	else: 
	    weighted='unweighted'

	# Reference model
	iref = 5   # no reference  
	for i in xrange( len(self.NGAmodel) ):
	    if refmodel == self.NGAmodel[i]: 
		iref = i
		break 
	if refmodel == 'CyberShake': 
	    iref = 4 

	TargetModels = self.NGAmodel + ['CyberShake']
	Factors = ['d RMS','c RMS','b RMS','a']
	clr=['r','g','b','k']

	# weighting function for directivity
	# read in file
	srcsh_file = Gkxmfs_pth + '/SourceRuptureHypoInfo'
	srcNmh = np.loadtxt( srcsh_file, dtype={'names':('sids','Nms','Nhs'),'formats':('i4','i4','i4')} )  
	try: 
	    Nk = len(srcNmh['sids'])
	    Nhs = srcNmh['Nhs']
	except:
	    Nk = 1
	    srcNmh['sids'] = [srcNmh['sids'],]
	    Nhs = [srcNmh['Nhs'],]

	fig = init_fig(num=1,figsize=(12,8),dpi=100)
	if len(self.periods) == 1: 
	    axs = init_subaxes( fig, subs=(1,1),basic=(0.9,0.7,0.9,0.9) )
	    abcd = ['(a)']

	if len(self.periods) == 4:
	    axs = init_subaxes( fig, subs=(2,2),basic=(0.9,0.7,1.2,0.7) )
	    abcd = ['(a)','(b)','(c)','(d)']

	for it in xrange( len(self.periods) ):

	    Ti = self.periods[it]
	    Tkey = '%.2f'%Ti 

	    # Dkxs
	    Ddata = []
	    for ik in xrange( Nk ):
		sid = srcNmh['sids'][ik]
		Nh = Nhs[ik]
		Ddata.append([])
		dkxs_pth0 = dkxs_pth1 + '/%s'%sid
		for ih in xrange(Nh):
		    dkxsFile = dkxs_pth0 + '/CyberShake.NGAs.%s.Source%s.Ih%s.dkxs'%(Tkey, sid, ih)
		    Alldata = np.loadtxt( dkxsFile )[:,2:]
		    Ddata0 = np.r_[ (Alldata[:,4+wt*8:5+wt*8] - Alldata[:, wt*8:4+wt*8]).T, Alldata[:,4+wt*8:5+8*wt].T ].T  # get capital D
		    if iref != 5: 
			Ddata0 = Ddata0 - Ddata0[:,iref:iref+1]    # residual
		    Ddata[ik].append( Ddata0 )

	    # Ddata has the dimension: [ik][ih][ista,imodel]
	    RMSd = []   # five elements
	    for imodel in xrange( len(TargetModels) ):
		RMS_k = 0
		for ik in xrange( Nk ):
		    Nh = Nhs[ik]
		    RMS_h = 0
		    for ih in xrange( Nh ):
			RMS_h += sum( Ddata[ik][ih][:,imodel]**2 ) / self.Nloc 
		    RMS_k += RMS_h / Nh
		RMSd.append( np.sqrt(RMS_k/Nk) )

	    # Cks
	    Cdata = []
	    for ik in xrange( Nk ):
		sid = srcNmh['sids'][ik]
		cksFile = cks_pth1 + '/CyberShake.NGAs.%s.Source%s.cks'%(Tkey, sid)
		Alldata = np.loadtxt( cksFile )[:,2:]
		Cdata0 = np.r_[(Alldata[:,4:5]-Alldata[:,:4]).T, Alldata[:,4:5].T].T
		if iref != 5: 
		    Cdata0 = Cdata0 - Cdata0[:,iref:iref+1] 
		Cdata.append( Cdata0 )
	    RMSc = []
	    for imodel in xrange( len(TargetModels) ):
		RMS_k = 0
		for ik in xrange( Nk ):
		    RMS_k += sum(Cdata[ik][:,imodel]**2) / self.Nloc 
		RMSc.append( np.sqrt(RMS_k/Nk) )

	    # Bs
	    bsFile = bs_pth1 + '/CyberShake.NGAs.%s.bs'%Tkey
	    Alldata = np.loadtxt( bsFile )[:,5:] 
	    Bdata = np.r_[(Alldata[:,4:5]-Alldata[:,:4]).T, Alldata[:,4:5].T].T
	    if iref != 5: 
		Bdata = Bdata - Bdata[:,iref:iref+1] 
	    RMSb = []
	    for imodel in xrange( len(TargetModels) ):
		RMSb.append( np.sqrt(sum(Bdata[:,imodel]**2)/self.Nloc) )

	    # A value
	    file_e0 = a_pth1 + '/CyberShake.NGAs.%s.a'%Tkey
	    tmp_e0 = np.loadtxt( file_e0, skiprows=1 )
	    tmp_e0[:-2] = tmp_e0[-2] - tmp_e0[:-2]   # NGA models
	    if iref != 5: 
		tmp_e0[:-1] = tmp_e0[:-1] - tmp_e0[iref]    # including CyberShake as the target model

	    RMS = np.r_[[RMSd],[RMSc],[RMSb],[tmp_e0[:-1].tolist()]]
	    RMS = np.array(RMS)

	    # write RMS into file for Matlab plot 
	    outfile = RMSpth + '/CyberShake.NGAs.periods.%s.Ref%s.RMS.%s.txt'%(sigmakey,refmodel,weighted) 
	    np.savetxt( outfile, RMS, '%f' )

	    # plot 
	    xt = np.arange( len(TargetModels) )

	    # bar plot
	    if 0:
		yt = np.arange( 4 )   # A,B,C,D
		xxt, yyt = np.meshgrid( xt, yt ) 
		xpos = xxt.flatten()   # Convert positions to 1D array
		ypos = yyt.flatten()
		zpos = np.zeros(4*5)
		dx = 0.5 * np.ones_like(zpos)
		dy = dx.copy()
		dz = RMS.flatten()

		ax = Axes3D(fig)
		for item in xrange( 4 ):
		    index = 5*item
		    ax.bar3d( xpos[index:index+5], ypos[index:index+5], zpos[index:index+5], \
		              dx[index:index+5], dy[index:index+5], dz[index:index+5], color=clr[item] )
		ax.w_xaxis.set_ticklabels(TargetModels)
		ax.w_yaxis.set_ticklabels(Factors)
		ax.set_zlabel('RMS')
		ax.set_xlabel('Models')
		ax.set_ylabel('Factors')

	    if 1: 
		ax = fig.add_axes( axs[it] )
		lines = []
		for item in xrange( 4 ):
		    line0 = ax.plot( xt, RMS[item,:], clr[item]+'-o' )
		    lines.append( line0 )

		ax.set_xlim([-1,5])

		ax.set_xticks(xt)
		ax.set_xticklabels( TargetModels )
		xticknames = plt.setp(ax,xticklabels=TargetModels)
		plt.setp(xticknames, rotation=0, fontsize=10)
		ax.set_xlabel('Models')

		ax.set_ylabel('Factor values')
		ax.set_ylim([-1,1])

		ax.text( -0.5, 0.7, abcd[it], fontsize=14)

		if it == 0:
		    lg = ax.legend( lines, Factors, loc='lower center' )
		    lg.draw_frame(False)
		    ltext = plt.gca().get_legend().get_texts()
		    plt.setp(ltext, fontsize=8)

		#ax.set_title( r'Reference = %s, period = %s, $\gamma$ = %s'%(refmodel, Tkey, '%.2f'%sigma) )
		plt.grid(True)

		if len(self.periods) == 1:
		    fig.savefig( plotpth + '/CyberShake.NGAs.%s.%s.Ref%s.RMS.%s.pdf'%(Tkey,sigmakey,refmodel,weighted), format='pdf')
		    fig.savefig( plotpth + '/CyberShake.NGAs.%s.%s.Ref%s.RMS.%s.%s'%(Tkey,sigmakey,refmodel,weighted, pfmt), format=pfmt )
		else:
		    fig.savefig( plotpth + '/CyberShake.NGAs.periods.%s.Ref%s.RMS.%s.pdf'%(sigmakey,refmodel,weighted), format='pdf')
		    fig.savefig( plotpth + '/CyberShake.NGAs.periods.%s.Ref%s.RMS.%s.%s'%(sigmakey,refmodel,weighted, pfmt), format=pfmt )




    # =====================
    # seperate analysis: 
    # =====================
    # test different weighting functions weights for directivity effects
    # This step is to make CyberShake more like NGA data
    def DkxsAnalysis(self, CyberShakeRvar, ngaP, sites_info, rups_info, Sources, ngaD=None):

	# PATH info
	Gkxmfs_pth= self.output_pth+'/GkxmfsDmap'      # path to save directivity analysis results (weighting functions)

	# file save path
	dkxs_pth = Gkxmfs_pth + '/Dkxs' 
	SDks_pth = Gkxmfs_pth + '/SDks' 

	plotpth0 = self.flatfile_plot_pth00 + '/Gkxmfs'
	plotpth = plotpth0 + '/Dkxs'
	for f in [Gkxmfs_pth, dkxs_pth, SDks_pth, plotpth0, plotpth]:
	    if not os.path.exists(f):
		os.mkdir(f)
	pfmt = 'eps'

	Nk = len( CyberShakeRvar )


	# top loop
	for it, Ti in enumerate( self.periods ): 
	    # loop over periods (it could be one period)
	    Tkey = '%.2f'%Ti
	    print 'Analysis at periods = %s sec'%(Tkey)

	    # Use one NGA models to correct for directivity
	    TargetModel = {}
	    Hkxs =  {}; Dkxs = {}
	    Gkxmfs = []
	    pdf_f = []; pdf_m = []
	    for ik in xrange( Nk ):

		Nm = len( CyberShakeRvar[ik] )
		tmpCS = []
		for irup in xrange( Nm ):
		    tmpCS.append( np.log( CyberShakeRvar[ik][irup][Tkey] ) )   # [irup][ih,islip,ista]
		tmp_CS = np.array( tmpCS )
		Nm,Nh,Nf,Nsta = tmp_CS.shape

		Gkxmfs.append( tmp_CS )
		del(tmpCS)

		# set up weighting functions for slip-distribution and stress drop
		pdf_f.append([])
		pdf_m.append([])

		# Generate the distribution of each variabiles 
		Mws, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
		AreaTmp = Nrow*Ncol  # in km
		mu = 3.87 + np.log10(AreaTmp) * 1.05   # somerville 2006 M-A relationship
		std = 0.2

		prob0 = 1./(np.sqrt(2.*np.pi)*std) * np.exp(-0.5*((Mws-mu)/std)**2)
		pdf_m[ik] = prob0 / sum(prob0)    # weighting should be summed up to 1
		pdf_m[ik] = pdf_m[ik].tolist()

		# slip distribution (just use uniform distribution)
		for ivar in xrange( Nf ):
		    pdf_f[ik].append( 1./Nf ) 

	    # weighting functions for source and sites
	    pdf_k = np.ones(Nk)/Nk
	    pdf_s = np.ones(self.Nloc)/self.Nloc

	    TargetModel['CyberShake'] = Gkxmfs 
	    del( Gkxmfs )

	    Hkxs['CyberShake'] = []
	    Dkxs['CyberShake'] = []

	    # NGA Models
	    for inga in xrange( len(self.NGAmodel) ):
		nga = self.NGAmodel[inga]
		TargetModel[nga] = []
		Hkxs[nga] = []
		Dkxs[nga] = []

		for ik in xrange( Nk ):
		    Nm, Nh, Nf, Nsta = TargetModel['CyberShake'][ik].shape
		    tmp = np.zeros( (Nm,Nh,Nsta) )
		    for ir in xrange( Nm ):
			for ista in xrange( Nsta ):
			    if ngaD != None:
				for ih in xrange( Nh ): 
				    tmp[ir,ih,ista] = np.log( ngaP[nga][ik][Tkey][ir][0][ista] ) + np.log( np.array(ngaD[nga][ik][Tkey][ir])[ih,ista] )
			    else:
				tmp[ir,:,ista] =  np.log( ngaP[nga][ik][Tkey][ir][0][ista]  ) 
		    TargetModel[nga].append( tmp )
		    del(tmp)

	    # decomposition 
	    models = ['CyberShake'] + self.NGAmodel 
	    sumSDk = {}
	    sigmas = np.arange( 0.1, 1.1, 0.1 )  # normalized std
	    for imodel, model in enumerate( models ):
		for ik in xrange( Nk ):
		    if model == 'CyberShake':   
			Nm, Nh, Nf, Nsta = TargetModel[model][ik].shape

			tmp_xms = np.average( TargetModel[model][ik], axis=2, weights = pdf_f[ik] )
			Hkxs0 = np.average( tmp_xms, axis=0, weights=pdf_m[ik] )
		    else:
			Nm, Nh, Nsta = TargetModel[model][ik].shape
			tmp_xms = TargetModel[model][ik]
			Hkxs0 = np.average( tmp_xms, axis=0, weights=pdf_m[ik] )

		    # sites
		    meta = sites_info[ik]
		    sites_run = meta.keys()
		    slon = []; slat = []
		    for i in xrange( Nsta ):
			slon.append( meta[sites_run[i]]['lon'] )
			slat.append( meta[sites_run[i]]['lat'] )
		    slon = np.array( slon )
		    slat = np.array( slat )

		    tmp = np.zeros((Nh,self.Nloc))
		    for ih in xrange( Nh ):
			tmp[ih,:] = interp( slon, slat, Hkxs0[ih,:], self.slon1d, self.slat1d, eps=self.eps, method=self.method )
		    Hkxs[model].append( tmp ) 
		    del( tmp )

		# loop over half width of the hypocenter location distribution
		sumSDk0 = []; Dkxs0 = []
		for isigma in xrange(len(sigmas)): 
		    sigma = sigmas[isigma]
		    pdf_x = []
		    for ik in xrange( Nk ):
			Nh = np.array(Hkxs[model][ik]).shape[0]
			xh = np.arange( Nh ) / float(Nh-1)
			mu = 0.5
			# find proper Gaussian width which minimize the variances (Sdks0) of dkxs!
			prob0 = np.exp(-0.5*((xh-mu)/sigma)**2)   # normalized Gaussian distribution
			pdf_x.append( prob0/ sum(prob0) )   # to make sure sum(prob) = 1
		    Dkxs00, SDks0, cks0, bs0, a0 = im_decom2( Hkxs[model], self.Nloc, pdf_x, pdf_k, pdf_s )
		    tmp = []
		    for ik in xrange( Nk ):
			tmp.append( sum( SDks0[ik] )/self.Nloc)
		    sumSDk0.append(tmp)
		    Dkxs0.append(Dkxs00)

		sumSDk[model] = np.array(sumSDk0) 
		Dkxs[model] = Dkxs0 

	    # we now have sumSDk[model][sigma][ik]

	    # compute residual d-map[model] and sumSdk[model] and plot
	    # plot sigmas and sumSDks
	    sumSdk = {}
	    for inga,nga in enumerate( self.NGAmodel ):

		# compute standard deviation 
		sumSdk0 = []
		for isigma in xrange( len(sigmas) ): 
		    sigma = sigmas[isigma]
		    tmp = []
		    for ik in xrange( Nk ):
			Nh = np.array(Hkxs[model][ik]).shape[0]
			xh = np.arange( Nh ) / float(Nh-1)
			mu = 0.5
			prob0 = np.exp(-0.5*((xh-mu)/sigma)**2)   # normalized Gaussian distribution
			pdf_x = prob0/ sum(prob0)    # to make sure sum(prob) = 1
			DkxsTmp = Dkxs['CyberShake'][isigma][ik] - Dkxs[nga][isigma][ik]  # residual d
			tmp.append( np.sqrt( np.average( DkxsTmp**2, axis=0, weights=pdf_x ) ) ) 
		    tmp1 = []
		    for ik in xrange(Nk):
			tmp1.append( sum(tmp[ik])/self.Nloc )
		    sumSdk0.append(tmp1)
		sumSdk[nga] = sumSdk0 
	    # we now have sumSdk[nga][sigma][ik]

	    fig = plt.figure(1) 
	    for ik in xrange( Nk ):
		sid = Sources[ik]
		fig.clf()
		for inga, nga in enumerate( self.NGAmodel ):
		    ax = fig.add_subplot( 2,2, inga+1 )
		    ax.plot( sigmas, np.array(sumSDk['CyberShake'])[:,ik], 'bo' )
		    ax.plot( sigmas, np.array(sumSDk[nga])[:,ik], 'ro' )
		    ax.plot( sigmas, np.array(sumSdk[nga])[:,ik], 'gx' )
		    ax.set_xlim( [0,1.1] )
		    if inga == 0: 
			ax.set_title('%s, Source:%s'%(nga,sid))
		    else:
			ax.set_title('%s'%nga)

		    if inga in [2,3]:
			ax.set_xlabel( r'$\sigma$' )
		    if inga in [0,2]:
			ax.set_ylabel( r'$\sigma_d$' )
		fig.savefig( plotpth + '/SigmaDkxs.NGAs.Source%s.%s.%s'%(sid,Tkey,pfmt), format=pfmt, dpi=300 )

	    # ===============
	    # save to files
	    # ===============
	    if 0:
		src_fid = open( Gkxmfs_pth + '/SourceInfo', 'w' )
		srcsh_fid = open( Gkxmfs_pth + '/SourceRuptureHypoInfo','w' )

		for ik in xrange( Nk ):
		    sid = Sources[ik]
		    Mws, Nrow, Ncol, rake, dip, ztor, zbom = rups_info[ik]
		    Nm,Nh,Nf,Nsta = TargetModel[ik].shape

		    # source info
		    src_fid.write('%s %s\n'%(sid,Nrow*Ncol))
		    srcsh_fid.write('%s %s %s\n'%(sid,Nm,Nh))

		    Nh = dkxs0[ik].shape[0]

		    # D-Map
		    dkxs_pth0 = dkxs_pth + '/%s'%sid
		    for f in [dkxs_pth0,]:
			if not os.path.exists(f):
			    os.mkdir(f)

		    for ih in xrange(Nh):
			dkxsFile = dkxs_pth0 + '/CyberShake.%s.%s.Source%s.Ih%s.dkxs'%(nga, Tkey, sid, ih)
			fid = open( dkxsFile, 'w' )
			for iloc in xrange(self.Nloc): 
			    fid.write( '%s %s %s\n'%(\
			        self.slon1d[iloc], self.slat1d[iloc],\
			        dkxs0[ik][ih,iloc],\
			    ))
			fid.close()

		    # Standard deviation of D-Map
		    # sigma k_xs
		    SDksFile = SDks_pth + '/CyberShake.%s.%s.Source%s.Sdks'%(nga,Tkey, sid)
		    fid = open( SDksFile, 'w' )
		    for iloc in xrange(self.Nloc): 
			fid.write( '%s %s %s\n'%(\
			    self.slon1d[iloc], self.slat1d[iloc],\
			    Sdks0[ik][iloc],\
			) )
		    fid.close()

		src_fid.close()
		srcsh_fid.close()



