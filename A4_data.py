from scipy.interpolate import RegularGridInterpolator, interp1d, LinearNDInterpolator # psi_ytchen
import math
import CoolProp.CoolProp as CP
import numpy as np

#--------------------------------------------------------------------------------------------------
class Data:
    # Here is the testing change at 20240621
    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'data' object created in B
    def __init__(self, reactor):
        for matdic in reactor.control.input['mat']:
            if matdic['type'] == 'h2o':
                # psi_ytchen: creat interpolation functions for fluid properties
                print('Constructing IF97::Water database...')
                fluid = 'IF97::Water'
                hmax = CP.PropsSI('H','T',1073.0,'P',100.0E6, fluid)
                hmin = CP.PropsSI('H','T',300.00,'P', 100.0E3, fluid)
                pmax = 22.00E6
                pmin = 100.0E3
                pnp = np.linspace(pmin,pmax,num=100)
                hnp = np.linspace(hmin,hmax, num=800)
                xnp = np.zeros( len(pnp)*len(hnp) )
                rnp = np.zeros( len(pnp)*len(hnp) )
                vnp = np.zeros( len(pnp)*len(hnp) )
                mnp = np.zeros( len(pnp)*len(hnp) )
                knp = np.zeros( len(pnp)*len(hnp) )
                cnp = np.zeros( len(pnp)*len(hnp) )
                tnp = np.zeros( len(pnp)*len(hnp) )
                px =  np.zeros( len(pnp)*len(hnp) )
                hy =  np.zeros( len(pnp)*len(hnp) )
                value=np.zeros([len(pnp)*len(hnp),7] )
                
                rls = np.zeros( len(pnp) )
                rgs = np.zeros( len(pnp) )
                mls = np.zeros( len(pnp) )
                mgs = np.zeros( len(pnp) )
                cls = np.zeros( len(pnp) )
                cgs = np.zeros( len(pnp) )
                kls = np.zeros( len(pnp) )
                kgs = np.zeros( len(pnp) )
                ts  = np.zeros( len(pnp) ) 
                hgl = np.zeros( len(pnp) ) 
                sgm = np.zeros( len(pnp) ) 
                val1d = np.zeros( [11, len(pnp)] ) 
                
                for i in range(len(pnp)):
                    rls[i] = CP.PropsSI('D','Q',0.0,'P', pnp[i], fluid)
                    rgs[i] = CP.PropsSI('D','Q',1.0,'P', pnp[i], fluid)
                    mls[i] = CP.PropsSI('V','Q',0.0,'P', pnp[i], fluid)
                    mgs[i] = CP.PropsSI('V','Q',1.0,'P', pnp[i], fluid)
                    cls[i] = CP.PropsSI('C','Q',0.0,'P', pnp[i], fluid)
                    cgs[i] = CP.PropsSI('C','Q',1.0,'P', pnp[i], fluid)
                    kls[i] = CP.PropsSI('L','Q',0.0,'P', pnp[i], fluid)
                    kgs[i] = CP.PropsSI('L','Q',1.0,'P', pnp[i], fluid)
                    ts[i]  = CP.PropsSI('T','Q',0.5,'P', pnp[i], fluid)
                    hls    = CP.PropsSI('H','Q',0.0,'P',pnp[i], fluid)
                    hgs    = CP.PropsSI('H','Q',1.0,'P',pnp[i], fluid)
                    hgl[i] = hgs - hls
                    sgm[i] = CP.PropsSI('I','Q',0.0,'P',pnp[i], fluid)
                    for j in range(len(hnp)):
                        ij = i*len(hnp) + j
                        xnp[ij] = (hnp[j]-hls)/hgl[i]
                        tnp[ij] = CP.PropsSI('T','P',pnp[i],'H',hnp[j], fluid)
                        if xnp[ij] < 0. or xnp[ij] > 1. :
                            try:
                                rnp[ij] = CP.PropsSI('D','P',pnp[i],'H',hnp[j], fluid)
                                mnp[ij] = CP.PropsSI('V','P',pnp[i],'H',hnp[j], fluid)
                                knp[ij] = CP.PropsSI('L','P',pnp[i],'H',hnp[j], fluid)
                                cnp[ij] = CP.PropsSI('C','P',pnp[i],'H',hnp[j], fluid)
                            except ValueError:
                                if xnp[ij] < 0.:
                                    rnp[ij] = rls[i]
                                    mnp[ij] = mls[i]
                                    knp[ij] = kls[i]
                                    cnp[ij] = cls[i]
                                else:
                                    rnp[ij] = rgs[i]
                                    mnp[ij] = mgs[i]
                                    knp[ij] = kgs[i]
                                    cnp[ij] = cgs[i]
                        else:
                            rnp[ij] = rls[i]*(1. - xnp[ij]) + rgs[i]*xnp[ij]
                            rnp[ij] = 1./(1./rls[i] + xnp[ij]*(1./rgs[i] - 1./rls[i]))
                            mnp[ij] = 1.0/( 1.0/mls[i] + xnp[ij]*(1.0/mgs[i] - 1.0/mls[i]) )
                            knp[ij] = kls[i]*(1. - xnp[ij]) + kgs[i]*xnp[ij]
                            cnp[ij] = cls[i]*(1. - xnp[ij]) + cgs[i]*xnp[ij]
                            
                        vnp[ij] = mnp[ij]/rnp[ij]
                        px[ij]  = pnp[i]
                        hy[ij]  = hnp[j]
                
                value[:,0] = rnp.copy()
                value[:,1] = vnp.copy()
                value[:,2] = mnp.copy()
                value[:,3] = knp.copy()
                value[:,4] = cnp.copy()
                value[:,5] = tnp.copy()
                value[:,6] = xnp.copy()
                fph_h2o = LinearNDInterpolator(list( zip(px, hy) ), value)
                self.fph_h2o = fph_h2o
                hpt_h2o = LinearNDInterpolator(list( zip(px, tnp) ), hy)
                self.hpt_h2o = hpt_h2o
                
                val1d[0 , :] = hgl
                val1d[1 , :] = mls
                val1d[2 , :] = mgs
                val1d[3 , :] = rls
                val1d[4 , :] = rgs
                val1d[5 , :] = kls
                val1d[6 , :] = kgs
                val1d[7 , :] = ts 
                val1d[8 , :] = cls
                val1d[9 , :] = cgs
                val1d[10, :] = sgm
                
                fp_h2o = interp1d(pnp, val1d)
                
                self.fp_h2o = fp_h2o
                
                print('IF97::Water database done.')
                
            if matdic['type'] == 'lbe':
                print('Constructing Lead-Bismuth Eutectic database...')
                tm0 = 398.00
                tmin = 3.0E2
                tmax = 2.1E3 # temperatures [K]
                tnp = np.linspace(tmin,tmax,num=100)
                hnp = 164.8*(tnp-tm0)-1.97e-2*(tnp**2-tm0**2)+4.167e-6*(tnp**3-tm0**3)+4.56e5*(1.0/tnp-1.0/tm0)
                rnp = 11065-1.293*tnp
                vnp = (4.94e-4*np.exp(754.1/tnp))/rnp
                cnp = 164.8-3.94e-2*tnp+1.25e-5*tnp*tnp-4.56e5/tnp/tnp
                knp = 3.284 + 1.617e-2*tnp-2.305e-6*tnp*tnp
                values = np.zeros([5, 100])
                values[0,:] = tnp.copy()
                values[1,:] = rnp.copy()
                values[2,:] = vnp.copy()
                values[3,:] = cnp.copy()
                values[4,:] = knp.copy()
                fh_lbe = interp1d(hnp, values)
                self.fh_lbe = fh_lbe
                ht_lbe = interp1d(tnp, hnp)
                self.ht_lbe = ht_lbe
                print('Lead-Bismuth Eutectic database done.')
                
            if matdic['type'] == 'na':
                print('Constructing Sodium database...')
                t1  = 371.00
                t2  = 2000.0
                t3  = 2503.70
                tnp1 = np.linspace(t1, t2, num = 100, endpoint = False)
                hnp1 =-3.6577e5 + 1.6582e3*tnp1 - 4.2375e-1*tnp1**2 + 1.4847e-4*tnp1**3 + 2.9926e6/tnp1
                tnp2 = np.linspace(t2, t3, num = 100, endpoint = True)
                dhg  = ( 393.37*(1.0 - tnp2/t3) + 4398.6*(1.0 - tnp2/t3)**(0.29302) )*1000.0
                hnp2 = 2128.4e3 + 864.96*tnp2 - 0.5*dhg
                tnp = np.append(tnp1, tnp2)
                hnp = np.append(hnp1, hnp2)
                rnp = 219.0 + 275.32*(1.0 - tnp/t3) + 511.58*(1.0 - tnp/t3)**0.5
                vnp = np.exp(-6.4406 - 0.3958*np.log(tnp) + 556.835/tnp)/rnp
                cnp = 1658.2 - 0.84750*tnp + 4.4541e-04*tnp**2 - 2.9926e6/tnp**2
                knp = 124.67 - 0.11381*tnp + 5.5226e-5*tnp**2 - 1.1842e-8*tnp**3
                values = np.zeros([5, 200])
                values[0,:] = tnp.copy()
                values[1,:] = rnp.copy()
                values[2,:] = vnp.copy()
                values[3,:] = cnp.copy()
                values[4,:] = knp.copy()
                fh_na = interp1d(hnp, values)
                self.fh_na = fh_na
                ht_na = interp1d(tnp, hnp)
                self.ht_na = ht_na
                print('Sodium database done.')
                
    #----------------------------------------------------------------------------------------------
    # material properties: self is a 'data' object created in B, inp is a dictionary of input data dependent on the material
    def matpro(self, inp):

        # he: helium gas
        if inp['type'] == 'he':
            t = inp['t']
            k = 2.639e-3*t**0.7085
            return {'k':k}

        # mox: mixed uranium-plutonium oxide fuel
        if inp['type'] == 'mox':
            t,b,por,pu,x = inp['t'],inp['b'],inp['por'],inp['pu'],inp['x']
            # density (kg/m3)
            rho = (11460*pu + 10960*(1 - pu)) * (1 - por)
            # specific heat (J/kg-K), D.L. Hagrman, et al., "MATPRO-version 11", TREE-NUREG-1280, Rev 1, Idaho National Engineering Laboratory (1980).
            cp = 15.496*(19.53*539**2 * math.exp(539/t) / (t**2 * (math.exp(539/t) - 1)**2) + 2*9.25e-04*t + 6.02e06*4.01e4 / (1.987*t**2) * math.exp(-4.01e4/(1.987*t)))
            # thermal conductivity (W/m-K), Y. Philipponneau, J. Nuclear Matter., 188 (1992) 194-197
            k = (1/( 1.528*math.sqrt(x+0.00931) - 0.1055 + 0.44*b + 2.855e-4*t ) + 76.38e-12*t**3) * (1-por)/(1+por)/0.864
            return {'rho':rho, 'cp':cp, 'k':k}

        # na: liquid sodium
        # + psi_ytchen: enthalpy-based Na properties 07.02.2024
        elif inp['type'] == 'na':
            # Sodium two-phase properties should be added here
            ini = inp['ini']
            if ini ==0: # psi_ytchen: calculation based on h, vectorized in 20240312
                h = inp['h']
                hnp = np.array(h)
                value = self.fh_na(hnp)
                tf   = value[0,:]
                rhof = value[1,:]
                visf = value[2,:]
                cpf  = value[3,:]
                kf   = value[4,:]
                return {'tf':tf, 'rhof':rhof, 'visf':visf, 'cpf':cpf, 'kf':kf}
            elif ini > 0: # psi_ytchen: initialization run
            # psi_ytchen: calculation based on t
                t = inp['t']
                # psi_ytchen: International Atomic Energy Agency "Sodium Coolant Handbook: Physical and Chemical Properties", IAEA-TECDOC-XXXX, 2018
                hl = self.ht_na(t)
                return {'h':hl} 
        # - psi_ytchen: enthalpy-based Na properties 07.02.2024
        # lbe: liquid lead and bismuth (55%wt Bi, 45%wt Pb)
        # + psi_ytchen: enthalpy-based LBE properties 07.02.2024
        elif inp['type'] == 'lbe':
            ini = inp['ini']
            if ini ==0: # psi_ytchen: not initialization run
            # psi_ytchen: calculation based on h
                h = inp['h']
                hnp = np.array(h)
                value = self.fh_lbe(hnp)
                tf   = value[0,:]
                rhof = value[1,:]
                visf = value[2,:]
                cpf  = value[3,:]
                kf   = value[4,:]
                return {'tf':tf, 'rhof':rhof, 'visf':visf, 'cpf':cpf, 'kf':kf}
            elif ini > 0: # psi_ytchen: initialization run
            # psi_ytchen: calculation based on t
                t = inp['t']
                # OECD, Nuclear Energy Agency, Handbook on Lead-bismuth Eutectic Alloy and Lead Properties, Materials Compatibility, Thermalhydraulics and Technologies, OECD, 2015. https://doi.org/10.1787/42dcd531-en.
                # psi_ytchen: enthalpy (J/kg): @400-1100K equation from "Handbook on Lead-bismuth Eutectic Alloy and Lead Properties", p.98, same as the following ones
                hl = self.ht_lbe(t)
                return {'h':hl}
        # - psi_ytchen: enthalpy-based LBE properties 07.02.2024
        # + psi_ytchen: enthalpy-based h2o properties 06.02.2024
        elif inp['type'] == 'h2o':
            fluid = 'IF97::Water'
            ini   = inp['ini']
            if ini == 0: # psi_ytchen: not initialization run
            # psi_ytchen: calculation based on h
                h = inp['h']
                p = inp['p']
                hnp = np.array(h)
                pnp = np.array(p)
                value = self.fph_h2o( (pnp, hnp) )
                rhof =  value[ : ,0]
                visf =  value[ : ,1]
                miuf =  value[ : ,2]
                kf   =  value[ : ,3]
                cpf  =  value[ : ,4]
                tf   =  value[ : ,5]
                xe   =  value[ : ,6]
                
                
                val1d = self.fp_h2o( pnp )
                hgl = val1d[0 , :]
                mls = val1d[1 , :]
                mgs = val1d[2 , :]
                rls = val1d[3 , :]
                rgs = val1d[4 , :]
                kls = val1d[5 , :]
                kgs = val1d[6 , :]
                ts  = val1d[7 , :]
                cls = val1d[8 , :]
                cgs = val1d[9 , :]
                sgm = val1d[10, :]
                
                if np.max(xe) < 0. :
                    rhol = rhof.copy()
                    miul = miuf.copy()
                    kl   = kf.copy()
                    cpl  = cpf.copy()
                    rhog = rgs.copy()
                    miug = mgs.copy()
                    kg   = kgs.copy()
                    cpg  = cgs.copy()
                elif np.min(xe) > 1. :
                    rhol = rls.copy()
                    miul = mls.copy()
                    kl   = kls.copy()
                    cpl  = cls.copy()
                    rhog = rhof.copy()
                    miug = miuf.copy()
                    kg   = kf.copy()
                    cpg  = cpf.copy()
                else:
                    rhol = np.zeros( len(xe) )
                    rhog = np.zeros( len(xe) )
                    miul = np.zeros( len(xe) )
                    miug = np.zeros( len(xe) )
                    kl   = np.zeros( len(xe) )
                    kg   = np.zeros( len(xe) )
                    cpl  = np.zeros( len(xe) )
                    cpg  = np.zeros( len(xe) )
                    for i in range( len(xe) ):
                        if xe[i] < 0.: 
                            rhol[i] = rhof[i]
                            rhog[i] = rgs[i]
                            miul[i] = miuf[i] 
                            miug[i] = mgs[i]
                            kl[i]   = kf[i]
                            kg[i]   = kgs[i]
                            cpl[i]  = cpf[i]
                            cpg[i]  = cgs[i]
                        elif xe[i] > 1.: 
                            rhog[i] = rhof[i]
                            rhol[i] = rls[i]
                            miug[i] = miuf[i]
                            miul[i] = mls[i]
                            kg[i]   = kf[i]
                            kl[i]   = kls[i]
                            cpg[i]  = cpf[i]
                            cpl[i]  = cls[i]
                        else:
                            rhol[i] = rls[i]
                            rhog[i] = rgs[i]
                            miul[i] = mls[i]
                            miug[i] = mgs[i]
                            kl[i]   = kls[i]
                            kg[i]   = kgs[i]
                            cpl[i]  = cls[i]
                            cpg[i]  = cgs[i]
                return {'tf':tf,'rhol':rhol,'rhog':rhog,'rhof':rhof,'miul':miul,'miug':miug,'miuf':miuf,'visf':visf,\
                'cpl':cpl,'cpg':cpg,'cpf':cpf,'kl':kl,'kg':kg,'kf':kf,'xe': xe,'ts':ts,'hgl':hgl,'sgm':sgm}
            elif ini > 0: # psi_ytchen: initialization run
            # psi_ytchen: calculation based on t
                t = inp['t']
                p = inp['p']
                hh = self.hpt_h2o(p,t)
                return {'h':hh}
        # - psi_ytchen: enthalpy-based h2o properties 06.02.2024
        # ss316: stainless steel type of 316
        elif inp['type'] == 'ss316':
            t = inp['t']
            # density (kg/m3): @300K equation from Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.117
            rho = 7954.
            # specific heat (J/kg-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100. Note that 1 mol of SS316 = 10.165 kg (https://www.webqc.org/molecular-weight-of-SS316.html) and 1 cal = 4.184 J
            cp = (6.181 + 1.788e-3*t)*10.165*4.184
            # thermal conductivity (W/m-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100.
            k = 9.248 + 1.571e-2*t
            return {'rho':rho, 'cp':cp, 'k':k}
            
        # 9Cr1Mo: tube material of SGTF tube and shell material
        elif inp['type'] == '9Cr1Mo':
            t = inp['t']
            # density (kg/m3): from Vkiram SGTF article
            rho = 7600.
            # specific heat (J/kg-K): from Vkiram SGTF article
            cp = 600.
            # thermal conductivity (W/m-K): from Vkiram SGTF article
            A = 28.568
            B =-0.17859e-2
            C =-0.39325e-4
            D = 0.22380e-6
            E =-0.31979e-9
            F = 0.12585e-12
            tc = t - 273.15
            k = A + B*tc + C*tc**(2.0) + D*tc**(3.0) + E*tc**(4.0) + F*tc**(5.0)
            return {'rho':rho, 'cp':cp, 'k':k}
            
            # 2.25Cr1Mo: tube material of ETEC tube and shell material
        elif inp['type'] == '2.25Cr1Mo':
            t = inp['t']
            # density (kg/m3): from Vkiram SGTF article
            rho = 7750.
            # specific heat (J/kg-K): from Vkiram SGTF article
            cp = 600.
            # thermal conductivity (W/m-K): from Vkiram SGTF article
            A = 36.177
            B = 0.54737e-2
            C = 0.35874e-4
            D =-0.24726e-6
            E = 0.32447e-9
            F =-0.11300e-12
            tc = t - 273.15
            k = A + B*tc + C*tc**(2.0) + D*tc**(3.0) + E*tc**(4.0) + F*tc**(5.0)
            return {'rho':rho, 'cp':cp, 'k':k}
            
            
            
        # AISI316_powder: This material properties is only NACIE-UP benchmark
        elif inp['type'] == 'AISI316_powder':
            t = inp['t']
            rho = 7954.
            cp = (6.181 + 1.788e-3*t)*10.165*4.184
            # This material prporety is from NACIE-UP benchmark spec
            k = max(0.3 + 0.005*(t - 273.15 - 200), 0.30) 
            return {'rho':rho, 'cp':cp, 'k':k}

        # bn: boron nitide
        elif inp['type'] == 'bn':
            t = inp['t']
            tc = t -273.15
            # density (kg/m3): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            rho = 2000.
            # specific heat (J/kg-K): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            cp = 800.
            # thermal conductivity (W/m-K): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            k = 25.578 - 2.416*math.log(tc)
            return {'rho':rho, 'cp':cp, 'k':k}

        # cu: copper
        elif inp['type'] == 'cu':
            t = inp['t']
            # density (kg/m3): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            rho = 8933.
            # specific heat (J/kg-K): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            cp = 385.
            # thermal conductivity (W/m-K): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            k = 401
            return {'rho':rho, 'cp':cp, 'k':k}

    #----------------------------------------------------------------------------------------------
    # Nusselt number: self is a 'data' object created in B, inp is a dictionary of input data dependent on the case
    def qfluxcal(self, inp):
        HTCMOD = 0 # psi_ytchen: Heat Transfer Mode
        material_type = inp['type']
        Re = inp['re']
        Re = max(Re, 1.0) # psi_ytchen: to avoid math domain error
        Pr = inp['pr']
        Tf = inp['tmp'][0]
        Tw = inp['tmp'][1]
        dh = inp['dhyd']
        kf = inp['prop']['kf']
        PP = inp['p']
        if material_type == 'h2o':
            xe = inp['prop']['xe']
            Ts = inp['prop']['ts']
            Pmpa = PP/1.0e6
            if xe < 0.0:
              nuLam = 4.36
              nuTurb = 0.0233*Re**(0.80)*Pr**(0.40)
              nuSp = max(nuLam, nuTurb)
              hSp = nuSp*kf/dh # psi_ytchen: single-phase heat transfer coefficient
              qSp = hSp*(Tw - Tf)
            # psi_ytchen: the fluid is subcooled liquid
              if Tw < Ts: 
              # psi_ytchen: single-phase liquid heat transfer
                qflux = qSp
              else: # psi_ytchen: Tw >= Ts condition
                qThom= 2000.0*math.exp(Pmpa/4.34)*(Tw - Ts)**(2.0)
                if qThom < qSp: # Tw < Tonb condition
                  qflux = qSp
                else:           # Tw >= Tonb condition
                  hgl = inp['prop']['hgl']
                  sgm = inp['prop']['sgm']
                  rhol= inp['prop']['rhol']
                  rhog= inp['prop']['rhog']
                  qchf = 0.14*hgl*(sgm*9.8*rhog**(2.0)*(rhol-rhog))**(0.25)
                  Tchf = Ts + (qchf*math.exp(-Pmpa/4.34)/2000.0)**(0.50) # Obtained by making qchf .eq. qThom
                  if Tw < Tchf: # Tonb <= Tw < Tchf condition
                    qflux = qThom 
                    HTCMOD = 1 # subcooled nucleate boiling
                  else: # Tw > Tchf condition
                    kg  = inp['prop']['kg']
                    cpg = inp['prop']['cpg']
                    miug= inp['prop']['miug']
                    if Pmpa < 9.0:
                      Tmfb = 557.90 + 44.1*Pmpa - 3.72*Pmpa**(2.0)
                    else:
                      Tmfb = 647.09 + 0.71*Pmpa
                    # Tmfb is Leidenfrost temperature of water
                    if Tw < Tmfb: # Tchf < Tw < Tmfb condition, subcooled transition boiling
                      qmfb = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol - rhog)/miug)**(1.0/3.0)*(Tmfb - Ts)
                      ntmp = math.log( (Tmfb-Ts)/(Tw-Ts) )/math.log( (Tmfb-Ts)/(Tchf-Ts) )
                      qflux = qmfb*(qchf/qmfb)**(ntmp)
                      HTCMOD = 2 # subcooled transition boiling
                    else: # Tmfb < Tw condition, subcooled stable film boiling
                      qflux = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol - rhog)/miug)**(1.0/3.0)*(Tw - Ts)
                      HTCMOD = 3 # subcooled film boiling
            elif xe < 1.0: # 0.0 < xe < 1.0
              hgl = inp['prop']['hgl']
              sgm = inp['prop']['sgm']
              rhol= inp['prop']['rhol']
              rhog= inp['prop']['rhog']
              kg  = inp['prop']['kg']
              cpg = inp['prop']['cpg']
              miug= inp['prop']['miug']
              miul= inp['prop']['miul']
              G = inp['Gtot']
              
              alpha = rhol*xe/(rhol*xe + rhog*(1.0-xe)) # void fraction
              reg = G*dh/miug
              Prg = miug*cpg/kg
              
              #xcr = min( (0.44 - 0.006*Pmpa)*G**(0.114), 0.99 ) #for 0.012 < dh and dh < 0.013
              
              #omega = G*miul/(sgm*rhol)*(rhol/rhog)**(1.0/3.0)
              #xcr = min( 0.30 + 0.70*math.exp(-45.0*omega)*(0.008/dh)**(0.15), 0.99 )
              XX = Pmpa/9.8
              xcr = ( 0.39+1.57*XX-2.04*XX**(2.0)+0.68*XX**(3.0) )*(G/1000.)**(-0.5)*(dh/0.008)**(0.15)
              xcr = min(xcr, 0.99)
              #xcr = 0.99
              
              if xe < xcr: # 0.0 <= xe < xcr, not dryout condition
                qchf = 0.14*hgl*(sgm*9.8*rhog**(2.0)*(rhol-rhog))**(0.25)*(1.0 - alpha) # qchf considerting dryout effect
                Tchf = Ts + (qchf*math.exp(-Pmpa/4.34)/2000.0)**(0.50) # calculate Tchf
                # calculate qDrycr
                alpcr = rhol*xcr/(rhol*xcr + rhog*(1.0-xcr)) # void fraction
                Ycr = 1.0 - 0.1*((rhol/rhog - 1.0)*(1.0 - xcr))**(0.40)
                
                nuDrycr = 0.0230*(reg*(xcr/alpcr))**(0.80)*Prg**(0.40)*Ycr*(Tf/Tw)**(0.50) # Mitropolsky
                # nuDrycr = 0.052*(reg*(xcr/alpcr))**(0.688)*Prg**(1.26)*Ycr**(-1.06)*(Tf/Tw)**(0.50) # Groeneveld Annular
                # nuDrycr = 1.09e-03*(reg*(xcr/alpcr))**(0.989)*Prg**(1.41)*Ycr**(-1.15)*(Tf/Tw)**(0.50) # Groeneveld Circular
                # nuDrycr = 3.27e-03*(reg*(xcr/alpcr))**(0.901)*Prg**(1.32)*Ycr**(-1.50)*(Tf/Tw)**(0.50) # Groeneveld General
                hDrycr = nuDrycr*kg/dh
                qDrycr = hDrycr*(Tw - Tf) # Totally dryout heat flux
                phy = math.exp( xe/(xe - xcr) ) # psi_ytchen: interpolate factor for trial and error 20240324
                if Tw < Tchf:
                  qflux = 2000.0*math.exp(Pmpa/4.34)*(Tw - Ts)**(2.0) # Thom correlation
                  HTCMOD = 4 # saturated nucleate boiling
                else: # Tchf < Tw
                  if Pmpa < 9.0:
                    Tmfb = 557.90 + 44.1*Pmpa - 3.72*Pmpa**(2.0)
                  else:
                    Tmfb = 647.09 + 0.71*Pmpa
                  if Tw < Tmfb:
                    qmfb = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol - rhog)/miug)**(1.0/3.0)*(Tmfb - Ts)
                    ntmp = math.log( (Tmfb-Ts)/(Tw-Ts) )/math.log( (Tmfb-Ts)/(Tchf-Ts) )
                    qflux = qmfb*(qchf/qmfb)**(ntmp)
                    HTCMOD = 5 # saturated transion boiling
                  else: # Tmfb < Tw condition, subcooled stable film boiling
                    qflux = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol - rhog)/miug)**(1.0/3.0)*(Tw - Ts)
                    HTCMOD = 6 # saturated film boiling
                qflux = qflux*phy + (1.0-phy)*qDrycr # qflux is interpolated for trial and error 20240324
              else: # xcr <= xe < 1.0, dryout condition
                Y = 1.0 - 0.1*((rhol/rhog - 1.0)*(1.0 - xe))**(0.40)
                nuDry = 0.0230*(reg*(xe/alpha))**(0.80)*Prg**(0.40)*Y*(Tf/Tw)**(0.50) # Mitropolsky ！！！0.4 or 0.80
                # nuDry = 0.052*(reg*(xe/alpha))**(0.688)*Prg**(1.26)*Y**(-1.06)*(Tf/Tw)**(0.50) # Groeneveld Annular
                # nuDry = 1.09e-03*(reg*(xe/alpha))**(0.989)*Prg**(1.41)*Y**(-1.15)*(Tf/Tw)**(0.50) # Groeneveld Circular
                # nuDry = 3.27e-03*(reg*(xe/alpha))**(0.901)*Prg**(1.32)*Y**(-1.50)*(Tf/Tw)**(0.50) # Groeneveld General
                hDry = nuDry*kg/dh
                qDry = hDry*(Tw - Tf) # Totally dryout heat flux
                qflux = qDry
                HTCMOD = 7 # saturated dryout boiling
            else: # xe >= 1.0 condition
              nuLam = 4.36
              nuTurb = 0.0230*Re**(0.80)*Pr**(0.40)*(Tf/Tw)**(0.50)
              #nuTurb = 0.0270*Re**(0.80)*Pr**(1.0/3.0)*(Tf/Tw)**(0.50)
              nuSp = max(nuLam, nuTurb)
              hSp = nuSp*kf/dh # psi_ytchen: single-phase heat transfer coefficient
              qflux = hSp*(Tw - Tf)
              HTCMOD = 8 # single-phase steam

        elif material_type == 'na' or material_type == 'lbe' :
            pe = Re*Pr
            if 'p2d' in inp:
                # pin bundle
                p2d = inp['p2d']
                # forced convection in a pin bundle (Mikityuk, NED 2008)
                nu = 0.047*(1.0-math.exp(-3.8*(p2d-1.0))) * ((pe)**0.77 + 250.0)
            else:
                # round tube
                nu = 4.8 + 0.025 * (pe)**0.8
            hex = nu*kf/dh
            qflux = hex*(Tw - Tf)
        return [qflux,HTCMOD]

    #----------------------------------------------------------------------------------------------
    # Friction factor: self is a 'data' object created in B, inp is a dictionary of input data dependent on the case
    def fricfac(self, inp):



        if 'p2d' in inp: # psi_ytchen: wire-wrapped rod bundles or barerod
            
            #Chen, S. K., et al. (2018). "The upgraded Cheng and Todreas correlation for pressure drop 
            #in hexagonal wire-wrapped rod bundles." Nuclear Engineering and Design 335: 356-373.
            # psi_ytchen: for calculation of a single constant, 'math' is faster than 'numpy'
            p2d = inp['p2d']
            
            RebL = 320.0*(10**(p2d-1.0))
            RebT = 10000.0*(10**(0.7*(p2d-1.0)))
            
            re = np.array(inp['re'])
            re[re < 1.] = 1.
            phi = np.log10(re/RebL)/math.log10(RebT/RebL)

            if 'h2d' in inp:
            # psi_ytchen: wire-wrapped rod bundle
                h2d = inp['h2d']
                CfbL = (-974.6 + 1612.0*p2d - 598.5*p2d**2)*math.pow(h2d, 0.06-0.085*p2d)
                CfbT = (0.8063 - 0.9022*math.log10(h2d) + 0.3526*(math.log10(h2d)**2))*p2d**9.7*math.pow(h2d, 1.78-2.0*p2d)
            else: 
            # psi_ytchen: barerod bundle
                if 1.0 < p2d and p2d <= 1.1:
                    LL = [26.0000, 888.2, -3334.0]
                    TT = [0.09378, 1.398, -8.664]
                elif p2d <=1.5:
                    LL = [62.97,  216.9,   -190.2]
                    TT = [0.1458, 0.03632, -0.03333]
                else:
                    print('warning, barerod correlation input out of range, check input')
                    
                CfbL = LL[0] + LL[1]*(p2d-1.) + LL[2]*(p2d-1.)**(2.0)
                CfbT = TT[0] + TT[1]*(p2d-1.) + TT[2]*(p2d-1.)**(2.0)
            if np.max(phi) < 0.0:
                # psi_ytchen: pure laminar friction factor
                return CfbL/re
            elif np.min(phi) > 1.0:
                # psi_ytchen: pure turbulent friction factor
                return CfbT/re**0.18
            else:
                # psi_ytchen: transition friction factor
                phi[phi < 0.] = 0.0
                phi[phi > 1.] = 1.0
                return CfbL/RebL*(1 - phi)**(1/3)*(1 - phi**7) + (CfbT/RebT**0.18)*phi**(1/3)
        
        else: # psi_ytchen: changed to Churchill correlation to accelerate calculation 
            re = np.array(inp['re'])
            re[re < 1.] = 1.
            a = (2.457*np.log(1.0/((7.0/re)**(0.9)) ))**16.0
            b = (37530./re)**16.0
            fric = 8.0*( (8.0/re)**(12) + 1.0/(a + b)**1.50 )**(1.0/12.0)
            # psi_ytchen:
            return fric
