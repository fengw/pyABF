#!/usr/bin/env python 
""" 
Disaggregation Tools
""" 
from utils import * 


def ExtractDisaggInfo(DisaggFile, NsCut=100, IML=0):
    """ 
    DisaggFile: original disaggregation file (100 Sources)
    For each given site (one)
    """
    
    if NsCut > 100: 
	print 'Total Source Number exceeds 100'
	raise ValueError 

    lines = open(DisaggFile, 'r').readlines() 
    if IML == 0: 
	spl = lines[1].strip().split('=')
	IML = float( spl[-1][:-1] )
	Prob = float( spl[1].strip().split('(')[0] )
    else: 
	spl = lines[1].strip().split('=')
	Prob = float( spl[-1][:-1] )
	IML = float( spl[1].strip().split('(')[0] )
    HazardCurveInfo = [IML,Prob] 

    SourceDisaggInfo = []; SourceIDs = []
    il = 0
    while il < len(lines): 
	spl = lines[il].strip() 
	if spl == 'Disaggregation Source List Info:': 
	    il += 1   # to header
	    for ipoint in xrange( NsCut ): 
		spl0 = lines[il+ipoint+1].strip().split()
		sid0 = int(spl0[0])
		prob0 = float( spl0[1] )  # in percentage
		
		# for each site, each source has its own probability to contribution the hazard level given by IML or Prob
		SourceIDs.append( sid0 )
		SourceDisaggInfo.append( [sid0,prob0] )    
	    il = il + NsCut + 1
	else: 
	    il += 1
    return SourceDisaggInfo, SourceIDs, HazardCurveInfo


def GetSourceIDs(DisaggSourceFile, Ntop=10):

    # read source id  (functionalized in disaggregation)
    lines = open(DisaggSourceFile, 'r').readlines() 
    sids = []; pdf_k = []               
    if Ntop > len(lines)-1: 
	print 'Ntop larger than total number of Sources from disaggregation'
	raise ValueError 
    
    for il in range(1,Ntop+1):          
	spl = lines[il].strip().split(',') 
	sids.append( int(spl[0]) )      
	pdf_k.append( float(spl[-1]) )   
    print 'Contribution of %s sources: %s %%'%(Ntop, '%.2f'%sum(pdf_k))
    pdf_k = np.array( pdf_k )           
    pdf_k = pdf_k / sum(pdf_k)  # normalized
    return sids, pdf_k 


class SourceDisagg: 
    """ 
    Source Disaggregation Tool
    """
    def __init__(self, wrk, \
	         rup_model_ids = (35,5,3,1), period=3.0, Prob='1.0E-4', IML='0_G', ftype='txt' ):
	
	# Basin parameters
	wrk = wrk
	self.datapth = os.path.join( wrk, 'data' )
	
	self.metapth = os.path.join( wrk, 'metadata' ) 
	self.DisaggOutput = os.path.join( self.metapth , 'Disaggregation' )
	self.DisaggMeta = os.path.join( self.DisaggOutput , 'DisaggMeta' ) 
	self.DisaggSources = os.path.join( self.DisaggOutput, 'DisaggSources' )
	self.DisaggHazardMaps = os.path.join( self.DisaggOutput, 'DisaggHazardMap' ) 

	self.plotpth = os.path.join( wrk, 'plots' ) 
	self.DisaggPlotpth = os.path.join( self.plotpth, 'Disagg_plots' )
	for f in [self.metapth, \
		  self.DisaggOutput, self.DisaggMeta, self.DisaggSources, self.DisaggHazardMaps, \
		  self.plotpth,self.DisaggPlotpth ]: 
	    if not os.path.exists( f ): 
		os.mkdir(f)
	
        self.erf_id, self.sgt_id, self.rv_id, self.vel_id = rup_model_ids
        
	self.Ti = period 
	self.Prob = Prob     # string
	self.IML = float(IML.strip().split('_')[0])       # number

	self.ftype = ftype 
	if self.IML == 0.0: 
	    self.DisaggStr = 'DisaggPOE_%s'%Prob    
	else: 
	    self.DisaggStr = 'DisaggIML_%s'%IML
        
	self.PeriodStr = '%s'%('%.0f'%self.Ti)+'sec'

        # get all site information (disaggregation over all sites)
	self.sitedata = self.datapth + '/Sites/cybershk_sites_info.dat' 
	self.SiteMetafile = self.DisaggMeta + '/CyberShakeSitesInfo.py'
	
	if not os.path.exists( self.SiteMetafile ): 
	    sid = 10; rid = 0 
	    sites_info = {}; Ek = {} 
	    fid_stdout = open('./stdout','w')
	    site_only = True 
	    sites = cybershk_sites(self.sitedata)   # all CyberShake Site Info
	    for stanam in sites.keys(): 
		cursor = db.cursor() 
		sites_info,Ek = \
		im_gen( cursor,sid,rid,stanam,\
			sites_info,Ek,fid_stdout,\
			rup_model_ids=rup_model_ids,site_only = site_only )
		cursor.close()
	    fid_stdout.close() 
	    meta = dict( sites_info = sites_info,)
	    save(self.SiteMetafile, meta, header='# CyberShake Subsampling Sites Info\n') 
	else: 
	    sites_info = load( self.SiteMetafile ).sites_info 
        self.sites_info = sites_info 

	# Output Metadata (include all informations: all sites)
	DisaggMetaFile = 'DisaggMeta_ERF%s_%s_SA_%s.py'%(self.erf_id,self.DisaggStr,self.PeriodStr) 
	self.DisaggMetaFile = self.DisaggMeta + '/' + DisaggMetaFile
        
	DisaggSourceFile = 'DisaggSources_ERF%s_%s_SA_%s.txt'%(self.erf_id, self.DisaggStr, self.PeriodStr)
	self.DisaggSourceFile = self.DisaggSources + '/' + DisaggSourceFile 
	
	DisaggHMFile = 'DisaggHazardMaps_ERF%s_%s_SA_%s.txt'%(self.erf_id, self.DisaggStr, self.PeriodStr)
	self.DisaggHMFile = self.DisaggHazardMaps + '/' + DisaggHMFile 

        self.NsCut = 100   # default cutoff number for sources to be shown 




    def RunSourceDisagg(self, SiteDec = 1): 

	output_dir = self.DisaggOutput + '/Sites' 
	if not os.path.exists(output_dir): 
	    os.mkdir( output_dir )
	
	# loop over all sites 
	sites_info = self.sites_info
	SourceDisaggInfos = {}
	SiteHazardCurveInfos = {}
	SourceIDsForAllSites = []
	
	Ns = len(sites_info.keys())
	Nsite = int( Ns/SiteDec ) 
	print 'Number of Sites: ', Nsite
	for isite, stanam in enumerate( sites_info.keys() ): 
	    if np.mod(isite,SiteDec)==0:   # sub-sample of sites
		run_id = sites_info[stanam]['run_id'] 
		outdir  = output_dir + '/%s'%stanam 
		if not os.path.exists( outdir ): 
		    os.mkdir( outdir )
		SiteID = sites_info[stanam]['id']
		filetmp = '/%s_ERF%s_Run%s_%s_SA_%s_*.txt'%(stanam,self.erf_id,run_id,self.DisaggStr,self.PeriodStr) 
		DisaggFiles = glob.glob( outdir + filetmp )
		if len(DisaggFiles) == 0:
		    #print 'Run Java Disaggregation...'
		    # note: instead of using --probs, you can use --imls to select iml (directly)
		    # and this Plotter use default imr as CyberShake? (I think so)
		    if self.IML == 0:
			cmd = 'DisaggregationPlotter --erf-id %s --sgt-var-id %s --rv-id %s --vel-model-id %s\
			    --run-id %s --period %s --probs %s --output-dir %s --type %s > %s/Disagg.out'%(\
			    self.erf_id,self.sgt_id,self.rv_id,self.vel_id, run_id, '%.2f'%self.Ti, self.Prob, outdir, self.ftype,self.metapth )
		    else:
			cmd = 'DisaggregationPlotter --erf-id %s --sgt-var-id %s --rv-id %s --vel-model-id %s\
				--run-id %s --period %s --imls %s --output-dir %s --type %s > %s/Disagg.out'%(\
				self.erf_id,self.sgt_id,self.rv_id,self.vel_id, run_id, '%.2f'%self.Ti, self.IML, outdir, self.ftype,self.metapth )
		    os.system( cmd ) 
		    DisaggFiles = glob.glob( outdir + '/%s_ERF%s_Run%s_%s_SA_%s_*.txt'%(stanam,self.erf_id,run_id,self.DisaggStr,self.PeriodStr) )
		DisaggFile = DisaggFiles[0] 
		SourceDisaggInfo, SourceIDsForOneSite, HazardCurveInfo = ExtractDisaggInfo(DisaggFile,IML=self.IML,NsCut =self.NsCut)
                SourceIDsForAllSites.append(SourceIDsForOneSite)
		SourceDisaggInfos[stanam] = SourceDisaggInfo
		SiteHazardCurveInfos[stanam] = [SiteID,HazardCurveInfo[0],HazardCurveInfo[1]]    # for plot SiteID, IML, Prob (point in hazard curve)
		
		# progress stutas
		sys.stdout.write('\rFinished %s site of %s'%(isite+1,Nsite))
		sys.stdout.flush() 
		sleep(0.01)
	    else: 
		continue
	sys.stdout.write('\n') 

	SourceIDs = [tmp for tmp1 in SourceIDsForAllSites for tmp in tmp1]   # merge into one list (all sourceIDs)
	
	# save for other usage
	meta = dict( 
		SourceDisaggInfos=SourceDisaggInfos, 
		SourceIDs = SourceIDs, 
		SiteHazardCurveInfos = SiteHazardCurveInfos,
		)
	save(self.DisaggMetaFile, meta, header='#Source IDs from disaggregation\n')
        print len(SourceIDs) 

    def AverageSourceDisagg(self, Ntop=10):
	"""
	Compute Averaged source disaggregation over all used sites
	Return top 10 sources with their normalized contributions
	"""
	meta = load(self.DisaggMetaFile) 
	SourceDisaggInfos = meta.SourceDisaggInfos 
	SourceIDs = meta.SourceIDs
	SiteHazardCurveInfos = meta.SiteHazardCurveInfos
	
	# compute averaged disaggregation over all sites
	SourceIDgroup = np.array( group_list( SourceIDs ) )
	SepSourceIDs = SourceIDgroup[:,2]    # all considered source id
	AveProb = []
	AveProbDict = {}
	Nsite = len(SiteHazardCurveInfos.keys())
	for isource in xrange( len(SepSourceIDs) ): 
	    sid = SepSourceIDs[isource] 
	    Probs = []
	    for isite, SiteName  in enumerate(SiteHazardCurveInfos.keys()): 
		tmp = np.array( SourceDisaggInfos[SiteName] )
		SIDs = tmp[:,0]
		try: 
		    index = (SIDs == sid).nonzero()[0][0]
		    Probs.append( tmp[index,1] )
		except: 
		    Probs.append( 0.0 )

            Probs = np.array( Probs )    # p(r|k)
	    
	    if sum(Probs) == 0: 
		Pr = 0.0 
	    else: 
		Pr = Probs / sum(Probs)    # distribution depends on the contribution 
	    Pr = 1./Nsite    #uniform distribution (will give you some bias, for example, at one station, one source contributes more and other less)

	    AveProb0 = sum( Pr * Probs )   # average over all sites with proper weights to get independent p(k) (and then sort to choose the top Ntop sources for ABF analysis)
	    
	    AveProb.append( AveProb0 )
	    AveProbDict['%s'%sid] = AveProb0 

	# Sort:
	AveWeights, SourceIDk, index = sorti( AveProb, list1 = SepSourceIDs, descent = True )

        # write into file
	UCERF2 = load( self.datapth + '/Sources/meta_erf35_sources.py' )
	U = UCERF2.UCERF 
	fid = open( self.DisaggSourceFile, 'w' )
	
	# you could add Number of different site type that are affected by each of those sources
	fid.write( '# SourceID  SourceName  Contribution(%)\n' )   
	sids = []; probs = [] 
	for i in xrange( self.NsCut ): 
	    skey = '%s'%SourceIDk[i]
	    if i < Ntop and AveProbDict[skey]>0.0: 
		#select top 10 sources
		sids.append( int(skey) ) 
		probs.append( AveProbDict[skey] )
	    fid.write('%s,  %s,  %s\n'%(skey, U[skey]['name'], AveProbDict[skey]))
	fid.close()
	probs = (np.array(probs)/sum(probs)).tolist()   # Normalized weights for Ntop sources from disaggregation (will be used in ABF analysis)
	
	return sids, probs 


    def PlotMatrix(self, Ntop=20, pfmt='pdf'): 
	"""
	Plot Disaggregation Matrix and hazard map
	"""
	meta = load(self.DisaggMetaFile) 
	SourceDisaggInfos = meta.SourceDisaggInfos 
	SourceIDs = meta.SourceIDs
	SiteHazardCurveInfos = meta.SiteHazardCurveInfos    # used to extract the hazard map (two types, depending on the condition)
	
	if Ntop > self.NsCut: 
	    print 'Ntop = %s, larger than Source Cutoff number %s'%(Ntop,self.NsCut)
	    raise ValueError 

	if self.IML == 0: 
	    ylab2 = 'SA %s s (g)'%('%.3f'%self.Ti)
	    title2 = 'probability of exceedance = %1.0e'%(float(self.Prob))
	else: 
	    ylab2 = 'Probability of Exceedance' 
	    title2 = '$\it {y}$ = %.2f (SA %.3f sec, g)'%(self.IML,self.Ti)

	fig = plt.figure(100,(12,12),dpi=300) 
	fig.clf() 

	# Plot all matrix
	ax1dim = [0.12,0.15,0.67,0.6]
	ax1 = fig.add_axes(ax1dim)
	plt.rc('font',family='Arial') 
	
	texts = np.arange( 1,101,11 )
	y = 1
	yt = 1.15  # text location
	x = np.linspace(0,1,len(texts))
	ax5dim = [0.12,0.76,0.67,0.1]
	ax5 = fig.add_axes(ax5dim)
	ax5.scatter( x,np.repeat(y,len(x)),c='k',s=5*texts,facecolor='none' )
	for i in xrange(len(x)):
	    ax5.text( x[i],yt,str(texts[i]),ha='center' )
	ax5.set_ylim([0.9,1.3])
	ax5.axis('off')
        
	#axcbar = fig.add_axes([0.73,0.25,0.01,0.4])

	Nsite = len(SiteHazardCurveInfos.keys())
	SiteIDs = []; Yvalues = []; SumProbOverSource = []
	for isite, SiteName  in enumerate(SiteHazardCurveInfos.keys()):
	    tmp = np.array( SourceDisaggInfos[SiteName] )
	    SIDs = tmp[:,0]
	    Probs = tmp[:,1]
	    SiteID = SiteHazardCurveInfos[SiteName][0]
	    SiteID0 = np.repeat( SiteID, len(SIDs) )
	    #img = ax1.scatter( SiteID0, SIDs, c=Probs, s=5*Probs, edgecolor='none') 
	    img = ax1.scatter( SiteID0, SIDs, c='k',s=5*Probs, facecolor='none') 
	    img.set_alpha(0.75)
	    if self.IML==0:
		Yvalue = SiteHazardCurveInfos[SiteName][1]
	    else: 
		Yvalue = SiteHazardCurveInfos[SiteName][2] 
	    Yvalues.append( Yvalue ) 
	    SiteIDs.append( SiteID )
	    SumProbOverSource.append( sum(Probs)/100. )
	    #img.set_clim([0.0, 100.0])
        #cbar = plt.colorbar(img,cax=axcbar) 
	#cbar.set_label('$\it{p}(\it{k}|\it{r})$ (%)')

	
	ax1.text( -0.07,1.05, '(a)', fontsize=18, transform=ax1.transAxes )
	
	#ax1.set_xlim( [0,450] )   # for colorful version
	ax1.set_ylim([0,300])
	ax1.set_ylabel('SourceIDs ($\it{k}$)') 
	ax1.set_xlabel('SiteIDs ($\it{r}$)' )
	#ax1.set_xticklabels([])
	ax1.invert_yaxis()
	
	ylim1 = ax1.get_ylim() 
	xlim1 = ax1.get_xlim()
	#ax1.grid(True) 

        if 0:
	    # plot SiteID and HazardCurveInfo 
	    ax2dim = [0.12,0.77,0.7,0.2]
	    ax2 = fig.add_axes(ax2dim)
	    
	    ax2.plot( SiteIDs, Yvalues, 'ko' )
	    ax2.text( 15, 0.008, '(a)' )
	    ax2.set_ylabel(ylab2)
	    ax2.set_title(title2)
	    #ax2.set_ylim([-0.002,0.01])
	    ax2.set_xlim(xlim1)
	    ax2.set_xticklabels([]) 
	    ax2.grid(True) 
	
        if 0:
	    # Plot Sum of all source contributions for each site 
	    ax3dim = [0.12,0.05,0.7,0.08]
	    ax3 = fig.add_axes(ax3dim)
	    ax3.xaxis.set_ticks_position('bottom')
	    ax3.set_xlabel('SiteIDs ($\it{r}$)') 
	    ax3.plot( SiteIDs, SumProbOverSource, 'bo' )
	    ax3.text( 15, 1.05, '(c)' )
	    ax3.set_ylabel('Summed\nContribution') 
	    ax3.set_ylim([0.9,1.1])
	    ax3.set_xlim(xlim1)
	    ax3.grid(True) 
	
	# AveProbOverSitesPerSource
	SourceIDgroup = np.array( group_list( SourceIDs ) )
	SepSourceIDs = SourceIDgroup[:,2]    # all considered source id
	AveProb = []
	AveProbDict = {}
	for isource in xrange( len(SepSourceIDs) ): 
	    sid = SepSourceIDs[isource] 
	    Probs = 0.0
	    for isite, SiteName  in enumerate(SiteHazardCurveInfos.keys()): 
		tmp = np.array( SourceDisaggInfos[SiteName] )
		SIDs = tmp[:,0]
		try: 
		    index = (SIDs == sid).nonzero()[0][0]
		    Probs += tmp[index,1]
		except: 
		    Probs += 0.0 
	    AveProb0 = Probs / Nsite
	    AveProb.append( AveProb0 )
	    AveProbDict['%s'%sid] = AveProb0 

	# Sort:
	AveWeights, SourceIDk, index = sorti( AveProb, list1 = SepSourceIDs, descent = True )
	
	print 'Sum of averaged weights: %s'%sum(AveWeights) 

	# plot averaged over site weights for each source
	ax4dim = [0.85,0.15,0.1,0.6]
	ax4 = fig.add_axes(ax4dim) 
	AveProb = np.array(AveProb) 
	ax4.plot( AveProb, SepSourceIDs, 'ko', mec='k', mfc='none' ) 
	Top20Weights = 0.0
	for i in xrange( Ntop ): 
	    # mark top 10 Sources
	    #ax4.text( AveWeights[i], SourceIDk[i], '%s'%SourceIDk[i], color='b')
	    Top20Weights += AveWeights[i] 
	    ax4.plot(  AveWeights[i], SourceIDk[i], 'ko', mec='k', mfc='k' )
	print 'Sum of top %s averaged weights: %s'%( Ntop, Top20Weights)

	#ax4.text( -4.9, 10, '(b)' )
	ax4.text( -0.2,1.05, '(b)', fontsize=18, transform=ax4.transAxes )
	ax4.invert_yaxis()
	ax4.xaxis.set_ticks_position('top')
	ax4.set_yticklabels([])
	ax4.set_ylim(ylim1) 
	ax4.set_xlim([-5,max(AveProb)+5])
	ax4.set_xlabel( 'Averaged\nContribution (%)' )
	#ax4.grid(True) 

	plotnam = 'DisaggMatrix_ERF%s_%s_SA_%s.%s'%(self.erf_id,self.DisaggStr,self.PeriodStr,pfmt) 
	fig.savefig( self.DisaggPlotpth + '/' + plotnam, format=pfmt,dpi=300)
	#plt.show() 


    def GridGen(self): 
	""" 
	Generate regular grid for GMP grdimage Plotter 
	""" 
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
	self.eps = 0.0
	smooth = 0.01
	self.method = {'name':'exp', 'smooth':smooth}

	self.cybershk_sitedata = self.datapth+'/Sites/cybershk_sites_info.dat'
        self.sites = cybershk_sites(self.cybershk_sitedata) 


    def GetHazardMap(self): 
	""" 
	Prepare files for GMT plot
	""" 
	self.GridGen() 

	meta = load(self.DisaggMetaFile) 
	SiteHazardCurveInfos = meta.SiteHazardCurveInfos    # used to extract the hazard map (two types, depending on the condition)
	Yvalues0 = []; SiteNames = []; lons = []; lats = []
	for isite, SiteName  in enumerate(SiteHazardCurveInfos.keys()):
	    if self.IML==0:
		Yvalue = SiteHazardCurveInfos[SiteName][1]
	    else: 
		Yvalue = SiteHazardCurveInfos[SiteName][2] 
	    Yvalues0.append( Yvalue )
	    SiteNames.append( SiteName ) 
	    lons.append( float(self.sites[SiteName]['lon']) ) 
	    lats.append( float(self.sites[SiteName]['lat']) ) 

	# interpolation to generate file: 
	lons = np.array( lons ) 
	lats = np.array( lats )
	Yvalues0 = np.array( Yvalues0 ) 

	Yvalues = interp( lons, lats, Yvalues0, self.slon1d, self.slat1d, eps=self.eps, method=self.method )
        
	fid = open(self.DisaggHMFile, 'w') 
	for iloc in xrange(len(self.slon1d)):
	    fid.write( '%s %s %s\n'%(self.slon1d[iloc],self.slat1d[iloc],Yvalues[iloc]) ) 
	fid.close() 


if __name__ == '__main__':
    inpth = '/Users/fengw/work/Project/CyberShake_analysis/metadata/Disaggregation/DisaggSources/UniformSiteDistribution' 
    file0 = 'DisaggSources_ERF35_DisaggIML_0.3_G_SA_3sec.txt'
    file0 = 'DisaggSources_ERF35_DisaggPOE_0.01_SA_3sec.txt'
    print file0
    file0 = inpth + '/' + file0 
    for Ntop in [20, 25, 30]:
	s1,s2 = GetSourceIDs( file0, Ntop )


