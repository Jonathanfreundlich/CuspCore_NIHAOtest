'''
FIELDS ADDED BY do_fits()

Notes: 
- The fields are based on the quantities defined in gl, ie, with reduced range if reduce_range_gl() was carried out

- When the reduced range is [], the outputs are also []

'''

# PROFILE DATA WITH derive_slopes()

data.update({['all','d','g','s']:
            {['lsfit_brho_b2_g3_constrained','lsfit_brho_b2_g3_unconstrained']:                   
            {'p':p,                                                              # Fit parameters (c, a, b, g, Rvir, Mvir)
             'brho':brho,                                                        # Modeled brho (brho = prf.brho(r, p, model))
             'rho':rho,                                                          # Modeled rho (rho = prf.rho(r, p, model))
             'bs':bs,                                                            # Modeled bs (bs = prf.bs(r, p, model))
             's':s,                                                              # Modeled s (s = prf.s(r, p, model))
             'M':M,                                                              # Modeled M (M = prf.M(r, p, model))
             'V':V,                                                              # Modeled V (V = prf.V(r, p, model))
             's1':prf.s(xi1*Rvir, p, model),                                     # Slope at (default) 0.015 Rvir
             's2':prf.s(xi2*Rvir, p, model),                                     # Slope at (default) Rvir
             'bs1':prf.bs(xi1*Rvir, p, model),                                   # bs slope at (default) 0.015 Rvir
             'bs2':prf.bs(xi2*Rvir, p, model),                                   # bs slope at (default) Rvir
             'c2':c2,                                                            # Concentration (c2 = Rvir/r[bisect(bs,2)-1])
             'rms':rms,                                                          # RMS error between modeled brho and data
             'Merr':(prf.calc_total_mass(r, rho, Rvir)-Mvir)/Mvir}}}             # Relative Mvir error 
            
            