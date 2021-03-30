'''
FIELDS ADDED BY derive_slopes() and derive_brho()

Notes: 
- The fields are based on the radii r defined when formatting (format_data from format_functions), 
  which is a regularly log-spaced array of (default) nbins=150 bins between Rmin=s['eps'].min() and Rmax=10*Rvir
  p = pynbody.analysis.profile.Profile(s, min=Rmin, max=Rmax, type='log', nbins=nbins, ndim=3)
  
- rho_smooth, logrho_smooth, logsigmar_smooth, logsigmar2_smooth, sigmar_smooth, sigmar2_smooth, alpha, beta_smooth, 
  gamma, s and bs are only defined for log10(r/Rvir) within rlim ([-2.,0.] by default), they are nan outside
  
- The Rvir entering in the previous comment is the imprecise one [Rvir=max(h1['r'])]

- reduce_range_gl() reduces the ranges to [] if all radii are nan (as is the case when using Rvir from Fangzhou and this Rvir 
  is a nan)

'''

# PROFILE DATA WITH derive_slopes()

data.update({['all','d','g','s']:{'M':noinf(M),                                  # Enclosed mass (M=ss[c]['M'])
                                  'Mall':noinf(Mall),                            # Total enclosed mass (Mall=ss['all']['M'])
                                  'rho':10**noinf(logrho),                       # Density profile
                                  'logrho':noinf(logrho),                        # "
                                  'rho_smooth':10**noinf(logrho_smooth),         # "
                                  'logrho_smooth':noinf(logrho_smooth),          # "
                                  'sigmar':noinf(sigmar),                        # Velocity dispersion (sigmar=ss[c]['vr_disp'])
                                  'logsigmar':noinf(logsigmar),                  # "
                                  'logsigmar2':noinf(logsigmar2),                # "
                                  'logsigmar_smooth':noinf(logsigmar_smooth),    # "
                                  'logsigmar2_smooth':noinf(logsigmar2_smooth),  # "
                                  'sigmar_smooth':10**noinf(logsigmar_smooth),   # " 
                                  'sigmar2_smooth':10**noinf(logsigmar2_smooth), # "
                                  'alpha':noinf(alpha),                          # Log slope of the density profile
                                  'beta_smooth':noinf(beta_smooth),              # Smoothed anisotropy parameter
                                  'gamma':noinf(gamma)}})                        # Log slope of the sigmar^2 profile



# PROFILE DATA WITH derive_brho()

data.update({['all','d','g','s']:{'dM':dM,                                       # Uncertainty on M (dM = M/sqrt(cumsum(n)))
                                  'drho':drho,                                   # Uncertainty on rho (drho=rho/sqrt(n)) 
                                  'brho':brho,                                   # Mean density
                                  'dbrho':dbrho,                                 # Uncertainty on brho (brho/sqrt(cumsum(n)))
                                  's':s,                                         # Logarithmic slope of rho
                                  'bs':bs,                                       # Logarithmic slope of brho
                                  'Rmax':Rmax,                                   # Largest radius r[-1]
                                  'Mmax':Mmax,                                   # Largest mass M[-1]
                                  'Rthr':Rthr}})                                 # Rthr = nan
             

# newkeys=['M','dM','Mall','rho','drho','logrho','rho_smooth','logrho_smooth','brho','dbrho','sigmar','logsigmar','logsigmar2','logsigmar_smooth','logsigmar2_smooth','sigmar_smooth','sigmar2_smooth','alpha','beta_smooth','gamma','s','bs','Rmax','Mmax','Rthr']