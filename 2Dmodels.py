"""
2D IM models

"""

import dadi

def splitExpMig (params, ns, pts):
	"""
	Model from Filatov et al. 2016, Molecular ecology 25(11): 2467-2481.
	following the split migration grows exponentially from mb to me
	mb = migr begins = migration right after population split
	me = migr ends = migration at present
	"""
	nu1,nu2,T,mb,me = params
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	#define a function to describe the non-constant migration rate.
	m_func = lambda t: mb*(me/mb)**(t/T)
	phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m_func, m21=m_func)
	fs = dadi.Spectrum.from_phi (phi,ns, (xx ,xx))
	return fs

def uSplitExpMig (params, ns, pts):
	"""
	 Modified from from Filatov et al. 2016, Molecular ecology 25(11): 2467-2481.
	following the split migration grows exponentially from mb to me
	mb = migr begins = migration right after population split
	me = migr ends = migration at present
	s - size of pop1 after split
	"""
	nu1,nu2,T,mb,me, s = params
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: s * (nu1/s)**(t/T) 
	nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T) 
	#define a function to describe the non-constant migration rate.
	m_func = lambda t: mb*(me/mb)**(t/T)
	phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m_func, m21=m_func)
	fs = dadi.Spectrum.from_phi (phi,ns, (xx ,xx))
	return fs

def AM(params, (n1,n2), pts):
	nu1, nu2, m12, m21, Tam, Ts = params
	"""
	Ancient Migration model, from Tine et al. 2014, Nature communications 5: 5770.
	nu1: Size of population 1 after split.
	nu2: Size of population 2 after split.
	m12: Migration from pop 2 to pop 1 (2*Na*m12).
	m21: Migration from pop 1 to pop 2.
	Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
	Ts: The scaled time between the end of ancient migration and present.
	n1,n2: Size of fs to generate.
	pts: Number of points to use in grid for evaluation.
	"""
	# Define the grid we'll use
	xx = dadi.Numerics.default_grid(pts)
	
	# phi for the equilibrium ancestral population
	phi = dadi.PhiManip.phi_1D(xx)
	# Now do the divergence event
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21 
	phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=m12, m21=m21)
	# We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
	phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
	
	# Finally, calculate the spectrum.
	fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
	return fs

def SC(params, (n1,n2), pts):
	"""
	Secondary Contact model, from Tine et al. 2014, Nature communications 5: 5770.
	nu1: Size of population 1 after split.
	nu2: Size of population 2 after split.
	m12: Migration from pop 2 to pop 1 (2*Na*m12).
	m21: Migration from pop 1 to pop 2.
	Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
	Tsc: The scale time between secondary contact and present.
	n1,n2: Size of fs to generate.
	pts: Number of points to use in grid for evaluation.
	"""
	nu1, nu2, m12, m21, Ts, Tsc = params
	# Define the grid we'll use
	xx = dadi.Numerics.default_grid(pts)
	
	# phi for the equilibrium ancestral population
	phi = dadi.PhiManip.phi_1D(xx)
	# Now do the divergence event
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
	phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
	# We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
	phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
	
	# Finally, calculate the spectrum.
	fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
	return fs


def AM2M(params, (n1,n2), pts):
	"""
	Ancient Migration model with two categories of loci experiencing different migration rates, from Tine et al. 2014, Nature communications 5: 5770.
	nu1: Size of population 1 after split.
	nu2: Size of population 2 after split.
	m12: Migration from pop 2 to pop 1 (2*Na*m12).
	m21: Migration from pop 1 to pop 2.
	me12: Effective migration from pop 2 to pop 1 in genomic islands.
	me21: Effective migration from pop 1 to pop 2 in genomic islands.
	Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
	Ts: The scaled time between the end of ancient migration and present.
	P: The porportion of the genome evolving neutrally
	n1,n2: Size of fs to generate.
	pts: Number of points to use in grid for evaluation.
	"""
	nu1, nu2, m12, m21, me12, me21, Tam, Ts, P = params
	# Define the grid we'll use
	xx = dadi.Numerics.default_grid(pts)
	
	### Calculate the neutral spectrum
	# phi for the equilibrium ancestral population
	phiN = dadi.PhiManip.phi_1D(xx)
	# Now do the divergence event
	phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
	# We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
	phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
	# We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
	phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
	# calculate the spectrum.
	fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
	
	### Calculate the genomic island spectrum
	# phi for the equilibrium ancestral population
	phiI = dadi.PhiManip.phi_1D(xx)
	# Now do the divergence event
	phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
	# We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
	phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
	# We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
	phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
	# calculate the spectrum.
	fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
	
	### Sum the two spectra in proportion P
	fs = P*fsN+(1-P)*fsI
	return fs

def SC2M(params, (n1,n2), pts):
	"""
	Secondary Contact model with two categories of loci experiencing different migration rates, from Tine et al. 2014, Nature communications 5: 5770. 
	nu1: Size of population 1 after split.
	nu2: Size of population 2 after split.
	m12: Migration from pop 2 to pop 1 (2*Na*m12).
	m21: Migration from pop 1 to pop 2.
	me12: Effective migration from pop 2 to pop 1 in genomic islands.
	me21: Effective migration from pop 1 to pop 2 in genomic islands.
	Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
	Tsc: The scale time between the secondary contact and present.
	P: The porportion of the genome evolving neutrally
	n1,n2: Size of fs to generate.
	pts: Number of points to use in grid for evaluation.
	"""
	nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P = params
	# Define the grid we'll use
	xx = dadi.Numerics.default_grid(pts)
	
	### Calculate the neutral spectrum
	# phi for the equilibrium ancestral population
	phiN = dadi.PhiManip.phi_1D(xx)
	# Now do the divergence event
	phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
	# We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
	phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
	# We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
	phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
	# calculate the spectrum.
	fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
	
	### Calculate the genomic island spectrum
	# phi for the equilibrium ancestral population
	phiI = dadi.PhiManip.phi_1D(xx)
	# Now do the divergence event
	phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
	# We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
	phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
	# We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
	phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
	# calculate the spectrum.
	fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
	
	### Sum the two spectra in proportion P
	fs = P*fsN+(1-P)*fsI
	return fs


def IM_2M(params, ns, pts): 
	""" 
	Isolation-with-migration model with exponential pop growth and 2 classes of sites 
	by B. Nevado, adapted from dadi models & Tine et al. 2014
	s: Size of pop 1 after split. (Pop 2 has size 1-s) 
	nu1: Final size of pop 1. 
	nu2: Final size of pop 2. 
	T: Time in the past of split (in units of 2*Na generations)  
	m1: Migration in site class 1 (2*Na*m12) 
	m2: Migration in site class 2
	P: proportion of class 1 sites
	n1,n2: Sample sizes of resulting Spectrum 
	pts: Number of grid points to use in integration. 
	""" 
	s,nu1,nu2,T,m1,m2, P = params 
	
	xx = dadi.Numerics.default_grid(pts) 

	nu1_func = lambda t: s * (nu1/s)**(t/T) 
	nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T) 

	# site class 1
	phi = dadi.PhiManip.phi_1D(xx) 
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi) 
	phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m1, m21=m1) 
	fs1 = dadi.Spectrum.from_phi(phi, ns, (xx,xx)) 
	
	# site class 2
	phi2 = dadi.PhiManip.phi_1D(xx) 
	phi2 = dadi.PhiManip.phi_1D_to_2D(xx, phi2) 
	phi2 = dadi.Integration.two_pops(phi2, xx, T, nu1_func, nu2_func, m12=m2, m21=m2) 
	fs2 = dadi.Spectrum.from_phi(phi2, ns, (xx,xx)) 

	### Sum the two spectra in proportion P
	fs = P*fs1+(1-P)*fs2
	return fs 

def IM_2M_SC(params, ns, pts): 
	""" 
	Isolation-with-migration model with exponential pop growth, period of allopatry followed by Sec contact, 
	and 2 classes of sites. Migration is symetric.
	by B. Nevado, adapted from dadi models & Tine et al. 2014
	s: Size of pop 1 after split. (Pop 2 has size 1-s.) 
	nu1: Final size of pop 1. 
	nu2: Final size of pop 2. 
	T1: Duration of allopatric phase (in units of 2*Na generations)  
	T2: Duration of secondary contact phase (in units of 2*Na generations)  
	m1: Migration in site class 1 (2*Na*m12) 
	m2: Migration in site class 2
	P: proportion of class 1 sites
	n1,n2: Sample sizes of resulting Spectrum 
	pts: Number of grid points to use in integration. 
	""" 
	s,nu1,nu2,T1, T2,m1,m2, P = params 
	
	xx = dadi.Numerics.default_grid(pts) 

	nu1_func = lambda t: s * (nu1/s)**(t/(T1+T2)) 
	nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/(T1+T2)) 

	# site class 1
	phi = dadi.PhiManip.phi_1D(xx) 
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi) 
	phi = dadi.Integration.two_pops(phi, xx, T1, nu1_func, nu2_func, m12=0, m21=0) 
	phi = dadi.Integration.two_pops(phi, xx, T2, nu1_func, nu2_func, m12=m1, m21=m1) 
	fs1 = dadi.Spectrum.from_phi(phi, ns, (xx,xx)) 
	
	# site class 2
	phi2 = dadi.PhiManip.phi_1D(xx) 
	phi2 = dadi.PhiManip.phi_1D_to_2D(xx, phi2) 
	phi2 = dadi.Integration.two_pops(phi2, xx, T1, nu1_func, nu2_func, m12=0, m21=0) 
	phi2 = dadi.Integration.two_pops(phi2, xx, T2, nu1_func, nu2_func, m12=m2, m21=m2) 
	fs2 = dadi.Spectrum.from_phi(phi2, ns, (xx,xx)) 

	### Sum the two spectra in proportion P
	fs = P*fs1+(1-P)*fs2
	return fs 

def IM_2M_AM(params, ns, pts): 
	""" 
	ns = (n1,n2) 
	params = (s,nu1,nu2,T1,T2, m1, m2, P) 
	
	Isolation-with-migration model with exponential pop growth, period of gene flow followed by allopatry, 
	and 2 classes of sites. Migration is symetric. 
	by B. Nevado, adapted from dadi models & Tine et al. 2014
	s: Size of pop 1 after split. (Pop 2 has size 1-s.) 
	nu1: Final size of pop 1. 
	nu2: Final size of pop 2. 
	T1: Duration of gene flow after split (in units of 2*Na generations)  
	T2: Duration of allopatric phase (in units of 2*Na generations)  
	m1: Migration in site class 1 (2*Na*m12) 
	m2: Migration in site class 2
	P: proportion of class 1 sites
	n1,n2: Sample sizes of resulting Spectrum 
	pts: Number of grid points to use in integration. 
	""" 
	s,nu1,nu2,T1, T2,m1,m2, P = params 
	
	xx = dadi.Numerics.default_grid(pts) 

	nu1_func = lambda t: s * (nu1/s)**(t/(T1+T2)) 
	nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/(T1+T2)) 

	# site class 1
	phi = dadi.PhiManip.phi_1D(xx) 
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi) 
	phi = dadi.Integration.two_pops(phi, xx, T1, nu1_func, nu2_func, m12=m1, m21=m1) 
	phi = dadi.Integration.two_pops(phi, xx, T2, nu1_func, nu2_func, m12=0, m21=0) 
	fs1 = dadi.Spectrum.from_phi(phi, ns, (xx,xx)) 
	
	# site class 2
	phi2 = dadi.PhiManip.phi_1D(xx) 
	phi2 = dadi.PhiManip.phi_1D_to_2D(xx, phi2) 
	phi2 = dadi.Integration.two_pops(phi2, xx, T1, nu1_func, nu2_func, m12=m2, m21=m2) 
	phi2 = dadi.Integration.two_pops(phi2, xx, T2, nu1_func, nu2_func, m12=0, m21=0) 
	fs2 = dadi.Spectrum.from_phi(phi2, ns, (xx,xx)) 

	### Sum the two spectra in proportion P
	fs = P*fs1+(1-P)*fs2
	return fs 
