'''
STRUCTURE OF THE DATA CREATED BY FORMAT_DATA FROM THE FORMAT_FUNCTIONS

Notes: 
- r is a regularly log-spaced array of (default) nbins=150 bins between Rmin=s['eps'].min() and Rmax=10*Rvir
  p = pynbody.analysis.profile.Profile(s, min=Rmin, max=Rmax, type='log', nbins=nbins, ndim=3)
  
- data['all']['r'], data['d']['r'], data['g']['r'], data['s']['r'] and data['flowdata']['r'] are similar

- quantities are defined for the following components: 'all','d' (dark matter),'g' (gas) and 's' (stars)

- in flowdata, such quantities are separated between 'in' (incoming) and 'out' (outflowing)
  
'''

# PROFILE DATA WITH getdata()
data = {    'sim':sim,                                                           # simulation name                         
            'a':a,                                                               # scale factor
            'z':z,                                                               # redshift
            't':t,                                                               # time
            'Rvir':Rvir,                                                         # virial radius (imprecise)
            'rho_crit':rho_crit,                                                 # critical density
            'fsv':fsv,                                                           # sum(h1.s['mass'])/(sum(h1['mass'])*fbar)
            'er':er,                                                             # h1.properties['Rmax']/max(h1['r']) 
            'Mstar':Mstar_tot,                                                   # stellar mass within 0.15 Rvir
            'SFR':SFR_tot,                                                       # SFR [40-80 Myr]
            'SFR_Ha':SFR_Ha_tot,                                                 # SFR [30-60 Myr]
            'SFR_UV':SFR_UV_tot,                                                 # SFR [80-120 Myr]
            ['all','d','g','s']:{'r':p['rbins'].in_units('kpc'),                 # radii
                                 'M':p['mass_enc'].in_units('Msol'),             # enclosed mass
                                 'n':cumsum(p['n']),                             # number of particles enclosed
                                 'Rvir':Rvir,                                    # Virial radius (imprecise), same as above
                                 'Mvir':float(sum(h1['mass'].in_units('Msol'))), # Virial mass (imprecise)
                                 'eps':min(s['eps']),                            # gravitational softening length
                                 'rho':p['density'].in_units('Msol kpc^-3'),     # density
                                 'vr_disp':p['vr_disp'].in_units('km s^-1'),     # radial velocity dispersion
                                 'vt_disp':p['vt_disp'].in_units('km s^-1'),     # tangential velocity dispersion
                                 'beta':p['beta']}}                              # anisotropy parameter

# STAR FORMATION BETWEEN TWO SNAPSHOTS WITH derive_SFRsnapshot()
data.update({'DM_snapshot':massform,
             'SFR_snapshot':SFR_snapshot,
             'SSFR_snapshot':SSFR_snapshot})

# FLOW DATA WITH derive_flowdata()

data.update({'flowdata':flowdata})

flowdata = {'r':r,                                                               # radii, same as data['all','d','g','s']['r']
            'z1':z1,                                                             # redshift of the previous timestep
            'a1':a1,                                                             # scale factor of the previous timestep
            't1':t1,                                                             # time of the previous timestep
            'z2':z2,                                                             # redshift of the current timestep
            'a2':a2,                                                             # scale factor of the current timestep
            't2':t2,                                                             # time of the current timestep
            'dt':dt,                                                             # time difference
            ['all','d','g','s']:{['in','out']:{'vr' = mean(vr),                  # flow velocity (at each r -- it is a list)
                                               'vr_ad' = mean(vr[ad]),           # adiabatic part's velocity (at each r)
                                               'vr_imp' = mean(vr[imp]),         # impulsive part's velocity (at each r)
                                               'vcirc' = vcirc,                  # circular velocity (at each r)
                                               'm' = sum(pclsc['mass']),         # flow mass (at each r)
                                               'm_ad' = sum(pclsc[ad]['mass']),  # adiabatic part's mass (at each r)
                                               'm_imp' = sum(pclsc[imp]['mass']) # impulsive part's mass (at each r)
                                               }}}
