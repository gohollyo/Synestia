### LOAD MODULES
import numpy as np
import struct
#from scipy.interpolate import griddata #hydrostruct
"""Classes and functions for accessing and manipulating tabulated EOS data."""
### Module for accessing and manipulating tabulated EOS data
### STS 09/2019
###
# calculate the structure for one planet
# make a class to hold the PREM data
class isentrope_class:
    """Class to isentrope data extracted from EOS table."""  # this is a documentation string for this class
    def __init__(self): # self is the default name of the object for internal referencing of the variables in the class
        """A function to initialize the class object.""" # this is a documentation string for this function
        self.ND = 0 # number of radius points
        self.S = 0.0 #entropy of isentrope
        self.density     = []   
        self.pressure    = []
        self.temperature = []
        self.soundspeed  = []
        self.energy  = []
        self.partvel = []
        self.region = [] # Tillotson region flag
        # not going to use all the variables in the file
        self.units = '' # I like to keep a text note in a structure about the units
#
    def loadisen(self,Sdisk,extEOStable):
        """A function to load data into class object."""
        ### ST Stewart Oct 16, 2019 Disk Structure Test
        # loop across all densities and extract the values for the requested isentrope
        # Sdisk must be in units of MJ/K/kg
        self.S=Sdisk #MJ/K/kg
        ForsteriteEOS=extEOStable
        print('\n ...LOADING ISENTROPE')
        print('Check that your input entropy (MJ K^-1 kg^-1) is not outside the interpolation range.')
        print('Smin and Smax')
        print(np.amin(ForsteriteEOS.S[np.where(ForsteriteEOS.S > 0)]), np.amax(ForsteriteEOS.S[np.where(ForsteriteEOS.S > 0)]))
        for i in range(0,ForsteriteEOS.ND):
            ind = np.where(ForsteriteEOS.S[:,i] > 0)[0]
            interpfunctionP = interp1d(ForsteriteEOS.S[ind,i],ForsteriteEOS.P[ind,i]) # pressure in GPa
            self.pressure = np.append(self.pressure,interpfunctionP(Sdisk)) # pressure in GPa
            interpfunctionT = interp1d(ForsteriteEOS.S[ind,i],ForsteriteEOS.T[ind]) # temperature in K
            self.temperature = np.append(self.temperature,interpfunctionT(Sdisk)) # temperature in K
        self.density = ForsteriteEOS.rho #density in g/cm3
        self.units = 'Units: S (MJ K^-1 kg^-1); density(g m^-3); pressure (GPa); temperature (K)'
        print('ISENTROPE LOADED. \n')
        show_plot=0
        if show_plot == 1:
            plt.figure()
            plt.subplot(221)
            plt.loglog(self.temperature,self.pressure,markersize='1')
            plt.ylabel('pressure')
            plt.xlabel('temperature')
            plt.subplot(222)
            plt.loglog(self.density,self.pressure,markersize='1')
            plt.ylabel('pressure')
            plt.xlabel('density')
            plt.subplot(223)
            plt.loglog(self.temperature,self.density,markersize='1')
            plt.ylabel('density')
            plt.xlabel('temperature')
            plt.suptitle('S=%.1e MJ/K/kg Forsterite ANEOS' % Sdisk)
            plt.tight_layout()
            plt.subplots_adjust(top=0.9)
            plt.savefig('isentrope_fors_aneos_%.1e.pdf'%Sdisk,bbox_inches='tight')
            plt.show()
#        
    def find_rho_T(self,P):
        """A function to find the density and temperature on an isentrope
        given pressure and entropy."""
        ### ST Stewart Oct 16, 2019 Disk Structure Test
        #it is assumed the pressure input will be in Pa (EOS units)
        #density and temperature will be output in mks units
        rhodisk = np.interp(P/1.e9,self.pressure,self.density)*1e3 #kg/m3
        Tdisk = np.interp(P/1.e9,self.pressure,self.temperature) #K
        return rhodisk,Tdisk
#        
class hydrostat:
    """Class to store hydrostatic P calculations"""
    ###GOHollyday Oct 17, 2019 Hydrostatic Structure Calculations
    def __init__(self,flag,rxy,z):
        """Initialize hydrostatic calc class."""
        nrxy = len(rxy)
        nz = len(z)
        self.Nrxy = nrxy #integer
        self.Nz = nz #integer
        self.flag = flag #0 or 1
        if self.flag==1:
            #for when want to calc based on pre-determined rxy,z grid
            self.RXY,self.Z = np.meshgrid(rxy,z,indexing='ij') #m #grids
            self.P = np.zeros((nrxy,nz)) #grid
            self.S = np.zeros((nrxy,nz)) #grid
            self.T = np.zeros((nrxy,nz)) #grid
            self.g_z = np.zeros((nrxy,nz)) #grid
            self.rho = np.zeros((nrxy,nz)) #grid
        else:
            #for when want to match Gadget_snapshot points
            self.RXY = rxy #m #1-D array
            self.Z = z #m #1-D array
            self.P = np.zeros((nrxy,1)) #1-D array
            self.S = np.zeros((nrxy,1)) #1-D array
            self.T = np.zeros((nrxy,1)) #1-D array
            self.g_z = np.zeros((nrxy,1)) #1-D array
            self.rho = np.zeros((nrxy,1)) #1-D array
        self.J2Ma2 = 0 #value for structure in kg m^2
        self.Mp = 0 #value for structure in kg
        self.units = 'Units: rxy, z, RXY, Z (m); g_z (m s^-2); P (Pa); S (J K^-1 kg^-1); rho (kg m^-3); T (K); Nrxy, Nz are integers.'
#        
    def loadhyd(self,LSQ,Smid,J2Ma2,Mp):
        """Load grid and midplane data into class."""
        #inputs should be in units of Pa, MJ K^-1 kg^-1, kg m^2, kg accordingly
        if self.flag == 1:
            self.P[:,0] = 10**LSQ(self.RXY[:,0]/1e6) #Pa
        else:
            self.P = 10**LSQ(self.RXY/1e6) #Pa
        self.S = Smid #MJ/K/kg
        self.J2Ma2 = J2Ma2 #kg m^2
        self.Mp = Mp #kg
#    
    def calc_hydrostat_equil(self,isentrope_class):
        """Calculate hydrostatic equilibrium."""
        ### ST Stewart Oct 16, 2019 Disk Structure Test   
        G = 6.67408e-11 #gravitational constant in m^3 kg^-1 s^-2
        M = self.Mp
        J2Ma2 = self.J2Ma2
        if self.flag == 0:
            zmid = 1e6 #m
            self.rho,self.T = isentrope_class.find_rho_T(self.P) #kg/m3, K
            for i in range(self.Nrxy):
                rhomid = self.rho[i] #kg/m3
                pmid = self.P[i] #Pa
                zmax = np.abs(self.Z[i]) #m
                if (zmax < (zmid+500.)):
                    step = zmax/10. #m
                    #nz = int((zmax-step)//step)
                    self.Z[i] = step #m
                else:
                    step = 500. #m
                    #nz = int((ztemp-zmid)//step)
                    self.Z[i] = zmid #m
                while self.Z[i] < zmax:
                    zdiff = zmax - self.Z[i] #m
                    if zdiff >= step:
                        dz = step
                    else:
                        dz = zdiff
                    gz = (G*self.Mp*self.Z[i]/((np.sqrt(self.RXY[i]**2 + self.Z[i]**2))**3)) - (3.*G*self.J2Ma2*((self.RXY[i]**2 + self.Z[i]**2)**-2.5)*self.Z[i]*(2.5*((self.Z[i]**2)/(self.RXY[i]**2 + self.Z[i]**2)) - 1.5))
                    self.g_z = np.append(self.g_z, gz)
                    dp = self.rho[i]*self.g_z[i]*dz # Pa
                    self.P[i] = self.P[i] - dp # Pa
                    self.rho[i],self.T[i]=isentrope_class.find_rho_T(self.P[i]) #kg m^-3, K
                    self.Z[i] = self.Z[i] + dz #m

        elif self.flag == -1:
            #trying to match GADGET particles
            for i in range(self.Nrxy):
                rxy = self.RXY[i]
                P = self.P[i]
                z = 0.1
                dz = 500.
                while ((z - np.abs(self.Z[i])) < dz):
                    rho,T=isentrope_class.find_rho_T(P)
                    gz=(G*M*z/((np.sqrt(rxy**2 + z**2))**3)) - (3.*G*J2Ma2*((rxy**2 + z**2)**-2.5)*z*(2.5*((z**2)/(rxy**2 + z**2)) - 1.5))
                    P = P - (rho*gz*dz)
                    z = z + dz
                if ((z - np.abs(self.Z[i])) != 0.0):
                    dz = (z - np.abs(self.Z[i]))
                    rho,T=isentrope_class.find_rho_T(P)
                    gz=(G*M*z/((np.sqrt(rxy**2 + z**2))**3)) - (3.*G*J2Ma2*((rxy**2 + z**2)**-2.5)*z*(2.5*((z**2)/(rxy**2 + z**2)) - 1.5))
                    P = P - (rho*gz*dz)
                    z = z + dz
                rho,T=isentrope_class.find_rho_T(P)
                gz=(G*M*z/((np.sqrt(rxy**2 + z**2))**3)) - (3.*G*J2Ma2*((rxy**2 + z**2)**-2.5)*z*(2.5*((z**2)/(rxy**2 + z**2)) - 1.5))
                self.P[i]=P
                self.T[i]=T
                self.rho[i]=rho
                self.g_z[i]=gz

        else:
            self.rho[:,0],self.T[:,0] = isentrope_class.find_rho_T(self.P[:,0]) #kg/m3, K
            for i in range(1,self.Nz):
                self.g_z[:,i-1] = (G*self.Mp*self.Z[:,i-1]/((np.sqrt(self.RXY[:,i-1]**2 + self.Z[:,i-1]**2))**3)) - (3.*G*self.J2Ma2*((self.RXY[:,i-1]**2 + self.Z[:,i-1]**2)**-2.5)*self.Z[:,i-1]*(2.5*((self.Z[:,i-1]**2)/(self.RXY[:,i-1]**2 + self.Z[:,i-1]**2)) - 1.5))
                self.P[:,i] = self.P[:,i-1] - (self.rho[:,i-1]*self.g_z[:,i-1]*(self.Z[:,i]-self.Z[:,i-1])) # Pa
                self.rho[:,i],self.T[:,i]=isentrope_class.find_rho_T(self.P[:,i]) #kg m^-3, K
        Pmask=np.ma.masked_less(self.P,0.)
        Pinval=np.ma.fix_invalid(Pmask)
        self.P=Pinval.filled(0.)
#
    def calc_pdiff(self,SPHrxy,SPHz,SPHp):
        #returns pressure in bars!!!
        P2d=griddata((SPHrxy,SPHz),SPHp/1e5,(self.RXY,self.Z),method='linear')
        P2dmask=np.ma.masked_less(P2d,0.)
        Pinval=np.ma.fix_invalid(P2dmask)
        P2d=Pinval.filled(0.)
        return (self.P/1e5)*(self.P/1e5-P2d),P2d
#
class EOShugoniot:
    """Class for Hugoniot array from extEOStable."""	
    def __init__(self):
        self.NH = 0
        self.rho = np.zeros(self.NH)   
        self.T = np.zeros(self.NH)   
        self.P = np.zeros(self.NH)   
        self.U = np.zeros(self.NH)   
        self.S = np.zeros(self.NH)   
        self.up = np.zeros(self.NH)   
        self.us = np.zeros(self.NH)
        self.cs = np.zeros(self.NH)
        self.units = ''
#
class EOSvaporcurve:
    """Class for vapor curve from ANEOS."""	
    def __init__(self):
        self.NT = 0
        self.NV = 0
        self.T = np.zeros(self.NT)  
        self.rl = np.zeros(self.NT)  
        self.rv = np.zeros(self.NT)  
        self.Pl = np.zeros(self.NT)  
        self.Pv = np.zeros(self.NT)  
        self.Ul = np.zeros(self.NT)  
        self.Uv = np.zeros(self.NT)  
        self.Sl = np.zeros(self.NT)  
        self.Sv = np.zeros(self.NT)
        self.Gl = np.zeros(self.NT)  
        self.Gv = np.zeros(self.NT)
        self.units = ''
#
class EOSmeltcurve:
    """Class for melt curve from ANEOS."""	
    def __init__(self):
        self.NT = 0
        self.NV = 0
        self.T  = np.zeros(self.NT)  
        self.rl = np.zeros(self.NT)  
        self.rs = np.zeros(self.NT)  
        self.Pl = np.zeros(self.NT)  
        self.Ps = np.zeros(self.NT)  
        self.Ul = np.zeros(self.NT)  
        self.Us = np.zeros(self.NT)  
        self.Sl = np.zeros(self.NT)  
        self.Ss = np.zeros(self.NT)
        self.units = ''
#
class EOS1barcurve:
    """Class for 1bar curve from the EOS."""	
    def __init__(self):
        self.NT    = 0
        self.S     = np.zeros(self.NT)  
        self.T     = np.zeros(self.NT)  
        self.rho   = np.zeros(self.NT)
        self.Tvap  = 0.
        self.Tmelt = 0.
        self.Sim   = 0.
        self.Scm   = 0.
        self.Siv   = 0.
        self.Scv   = 0.
        self.rhoiv   = 0.
        self.rhocv   = 0.
        self.rhocm   = 0.
        self.rhoim   = 0.
        self.units = ''
#
class EOScriticalpoint:
    """Class for critical point state from the EOS."""	
    def __init__(self):
        self.P   = 0
        self.S   = 0  
        self.T   = 0 
        self.rho = 0
        self.U   = 0
        self.units = ''
#
class EOStriplepoint:
    """Class for triple point state from the EOS."""	
    def __init__(self):
        self.P   = 0
        self.T   = 0 
        self.Sim   = 0.
        self.Scm   = 0.
        self.Siv   = 0.
        self.Scv   = 0.
        self.rhol  = 0.
        self.units = ''
#
class EOSaneoshugoniot:
    """Class for Hugoniot calculated in ANEOS."""	
    def __init__(self):
        self.ND  = 0
        self.NV  = 0
        #self.all = np.zeros((self.ND,self.NV))
        self.rho = 0
        self.T   = 0
        self.P   = 0
        self.U   = 0
        self.S   = 0
        self.us  = 0
        self.up  = 0
        self.units = ''
#
class extEOStable:
    """Class for accessing EXTENDED SESAME-STYLE EOS tables output from ANEOS"""
    #     ANEOS KPA FLAG
    #                                TABLE          ANEOS
    #     KPAQQ=STATE INDICATOR      =1, 1p    =1, 1p    (eos without melt)
    #                                =2, 2p lv =2, 2p liquid/solid plus vapor
    #                                          =4, 1p solid  (eos with melt)
    #                                          =5, 2p melt   (eos with melt)
    #                                          =6, 1p liquid (eos with melt)
    #                                =-1 bad value of temperature
    #                                =-2 bad value of density
    #                                =-3 bad value of material number
    #
    def __init__(self):
        self.ND  = 0 # integer; number of density points in grid
        self.NT  = 0 # integer; number of temperature points in grid
        self.rho = np.zeros(self.ND)          # g/cm3, density values
        self.T   = np.zeros(self.NT)          # K, temperature values
        self.P   = np.zeros(self.ND*self.NT)  # GPA, pressure(T,rho)
        self.U   = np.zeros(self.ND*self.NT)  # MJ/kg, sp. internal energy(T,rho)
        self.A   = np.zeros(self.ND*self.NT)  # MJ/kg, Helmholtz free energy(T,rho)
        self.S   = np.zeros(self.ND*self.NT)  # MJ/K/kg, sp. entropy(T,rho)
        self.cs  = np.zeros(self.ND*self.NT)  # cm/s, sound speed(T,rho)
        self.cv  = np.zeros(self.ND*self.NT)  # MJ/K/kg, sp. heat capacity(T,rho)
        self.KPA = np.zeros(self.ND*self.NT)  # integer, ANEOS KPA flag(T,rho)
        self.MDQ = np.zeros(self.ND*self.NT)  # integer, Model Development Quality Flag(T,rho)
        self.units = ''
        self.hug = EOShugoniot()
        self.vc  = EOSvaporcurve()
        self.mc  = EOSmeltcurve()
        self.cp  = EOScriticalpoint()
        self.tp  = EOStriplepoint()
        self.onebar = EOS1barcurve()
        self.anhug = EOSaneoshugoniot()
        # these are variables needed for the sesame header
        self.MATID   = 0.
        self.DATE    = 0.
        self.VERSION = 0.
        self.FMN     = 0.
        self.FMW     = 0.
        self.R0REF   = 0.
        self.K0REF   = 0.
        self.T0REF   = 0.
        self.P0REF   = 0.
        # variables needed for the ANEOS gamma function
        self.gamma0 = 0.
        self.theta0 = 0.
        self.C24    = 0.
        self.C60    = 0.
        self.C61    = 0.
        self.beta   = 0.
        # model name/version string
        self.MODELNAME = ''

    def loadstdsesame(self, fname, unitstxt=None):
        """Function for loading STD SESAME-STYLE EOS table output from ANEOS"""
        data = ([])
        if unitstxt is None:
            self.units = 'Units: rho g/cm3, T K, P GPa, U MJ/kg, A MJ/kg, S MJ/K/kg, cs cm/s, cv MJ/K/kg, KPA flag. 2D arrays are (NT,ND).'
        else:
            self.units = unitstxt
        sesamefile = open(fname,"r")  
        sesamedata=sesamefile.readlines()
        sesamefile.close()
        nskip = 6 # skip standard header to get to the content of the 301 table
        # num.density, num. temps
        tmp = sesamedata[nskip][0:16]
        dlen = float(tmp)
        tmp = sesamedata[nskip][16:32]
        tlen = float(tmp)
        if (np.mod((dlen*tlen*3.0+dlen+tlen+2.0),5.0) == 0):
            neos = int((dlen*tlen*3.0+dlen+tlen+2.0)/5.0) 
        else:
            neos = int((dlen*tlen*3.0+dlen+tlen+2.0)/5.0) +1
        #print(dlen,tlen,neos,len(sesamedata))
        data = np.zeros((neos,5),dtype=float)
        for j in range(nskip,neos+nskip):
            tmp3 = sesamedata[j]
            tmp4 = list(tmp3.split())
            if len(tmp4) < 5:
                lentmp4 = len(tmp4)
                data[j-nskip,0:lentmp4] = np.asarray(tmp4[0:lentmp4])
            else:
                data[j-nskip,:] = np.asarray(tmp4)
            #print(j,eosarr[j,:])
        #print(data.shape)
        data=np.resize(data,(neos*5))
        #print(data.shape)
        self.ND  = data[0].astype(int)  # now fill the extEOStable class
        self.NT  = data[1].astype(int)
        self.rho = data[2:2+self.ND]
        self.T   = data[2+self.ND : 2+self.ND+self.NT]
        self.P   = data[2+self.ND+self.NT : 2+self.ND+self.NT+self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.U   = data[2+self.ND+self.NT+self.ND*self.NT
                            : 2+self.ND+self.NT+2*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.A   = data[2+self.ND+self.NT+2*self.ND*self.NT
                            : 2+self.ND+self.NT+3*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        #self.S   = data[2+self.ND+self.NT+3*self.ND*self.NT
        #                    : 2+self.ND+self.NT+4*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.cs  = data[2+self.ND+self.NT+4*self.ND*self.NT
        #                    : 2+self.ND+self.NT+5*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.cv  = data[2+self.ND+self.NT+5*self.ND*self.NT
        #                    : 2+self.ND+self.NT+6*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.KPA = data[2+self.ND+self.NT+6*self.ND*self.NT
        #                    : 2+self.ND+self.NT+7*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
#
    def loadextsesame(self, fname, unitstxt=None):
        """Function for loading EXTENDED SESAME-STYLE EOS table output from ANEOS"""
        data = ([])
        if unitstxt is None:
            self.units = 'Units: rho g/cm3, T K, P GPa, U MJ/kg, A MJ/kg, S MJ/K/kg, cs cm/s, cv MJ/K/kg, KPA flag. 2D arrays are (NT,ND).'
        else:
            self.units = unitstxt
        sesamefile = open(fname,"r")  
        sesamedata=sesamefile.readlines()
        sesamefile.close()
        nskip = 6 # skip standard header to get to the content of the 301 table
        # num.density, num. temps
        tmp = sesamedata[nskip][0:16]
        dlen = float(tmp)
        tmp = sesamedata[nskip][16:32]
        tlen = float(tmp)
        if (np.mod((dlen*tlen*4.0+dlen+tlen+2.0),5.0) == 0):
            neos = int((dlen*tlen*4.0+dlen+tlen+2.0)/5.0)
        else:
            neos = int((dlen*tlen*4.0+dlen+tlen+2.0)/5.0) +1
        #print(dlen,tlen,neos,len(sesamedata))
        data = np.zeros((neos,5),dtype=float)
        for j in range(nskip,neos+nskip):
            tmp3 = sesamedata[j]
            tmp4 = list(tmp3.split())
            if len(tmp4) < 5:
                lentmp4 = len(tmp4)
                data[j-nskip,0:lentmp4] = np.asarray(tmp4[0:lentmp4])
            else:
                data[j-nskip,:] = np.asarray(tmp4)        
            #print(j,eosarr[j,:])
        #print(data.shape)
        data=np.resize(data,(neos*5))
        #print(data.shape)
        self.ND  = data[0].astype(int)  # now fill the extEOStable class
        self.NT  = data[1].astype(int)
        self.rho = data[2:2+self.ND]
        self.T   = data[2+self.ND : 2+self.ND+self.NT]
        #self.P   = data[2+self.ND+self.NT : 2+self.ND+self.NT+self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.U   = data[2+self.ND+self.NT+self.ND*self.NT
        #                    : 2+self.ND+self.NT+2*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.A   = data[2+self.ND+self.NT+2*self.ND*self.NT
        #                    : 2+self.ND+self.NT+3*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        self.S   = data[2+self.ND+self.NT+0*self.ND*self.NT
                            : 2+self.ND+self.NT+1*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.cs  = data[2+self.ND+self.NT+1*self.ND*self.NT
                            : 2+self.ND+self.NT+2*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.cv  = data[2+self.ND+self.NT+2*self.ND*self.NT
                            : 2+self.ND+self.NT+3*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.KPA = data[2+self.ND+self.NT+3*self.ND*self.NT
                            : 2+self.ND+self.NT+4*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
#
    def calchugoniot(self, r0=None, t0=None, pmax=None, writefilename=None):
        """Function for calculating a Hugoniot from EXTENDED SESAME-STYLE EOS table."""
        if r0 is None:
            return 'Must provide r0 and t0.'
        if t0 is None:
            return 'Must provide r0 and t0.'
        if pmax is None:
            pmax=1.E4 # GPa
        self.hug.rho = []
        self.hug.P = []
        self.hug.T = []
        self.hug.U = []
        self.hug.S = []
        self.hug.up = []
        self.hug.us = []
        self.hug.cs = []
  
        it0 = int(np.round(np.interp(t0,self.T,np.arange(self.NT)))) # uses nearest value if t0 not in array
        ir0 = int(np.round(np.interp(r0,self.rho,np.arange(self.ND)))) # uses nearest value if r0 not in the array
        p0  = self.P[it0,ir0] # GPa
        #print(self.P[it0,ir0])
        e0  = self.U[it0,ir0]#np.interp(p0,self.P[it0,:],self.U[it0,:])
        s0  = self.S[it0,ir0]#np.interp(p0,self.P[it0,:],self.S[it0,:])
        up0 = 0. # no initial particle velocity
        us0 = self.cs[it0,ir0]/1.e5 # cm/s->km/s use sound velocity for initial
        cs0 = self.cs[it0,ir0]/1.e5 # cm/s->km/s use sound velocity for initial
        #print(ir0,it0,r0,t0,p0,e0,up0,us0)
        self.hug.rho = np.append(self.hug.rho, self.rho[ir0])
        self.hug.P = np.append(self.hug.P, p0)
        self.hug.T = np.append(self.hug.T, self.T[it0])
        self.hug.U = np.append(self.hug.U, e0)
        self.hug.S = np.append(self.hug.S, s0)
        self.hug.up = np.append(self.hug.up, up0)
        self.hug.us = np.append(self.hug.us, us0)
        self.hug.cs = np.append(self.hug.cs, cs0)
        
        #for iir in range(ir0+1,self.ND):
        iir=ir0+1
        pnew=p0
        while pnew<pmax:
            ediff =0.5*(self.P[it0::,iir]+p0)*(1./r0-1./self.rho[iir])+e0 -(self.U[it0::,iir])  # MJ/kg
            # np.interp wants x values increasing
            pnew = np.interp(0.,np.flipud(ediff),np.flipud(self.P[it0::,iir]))
            tnew = np.interp(0.,np.flipud(ediff),np.flipud(self.T[it0::]))
            enew = np.interp(0.,np.flipud(ediff),np.flipud(self.U[it0::,iir]))
            snew = np.interp(0.,np.flipud(ediff),np.flipud(self.S[it0::,iir]))
            upnew = np.sqrt((pnew-p0)*(1./r0-1./self.rho[iir]))
            usnew = (1./r0)*np.sqrt((pnew-p0)/(1./r0-1./self.rho[iir]))
            csnew = np.interp(0.,np.flipud(ediff),np.flipud(self.cs[it0::,iir]))/1.E5 # km/s
            #print(self.rho[iir],tnew,pnew,enew,upnew,usnew)
            self.hug.rho = np.append(self.hug.rho, self.rho[iir])
            self.hug.P = np.append(self.hug.P, pnew)
            self.hug.T = np.append(self.hug.T, tnew)
            self.hug.U = np.append(self.hug.U, enew)
            self.hug.S = np.append(self.hug.S, snew)
            self.hug.up = np.append(self.hug.up, upnew)
            self.hug.us = np.append(self.hug.us, usnew)
            self.hug.cs = np.append(self.hug.cs, csnew) # km/s
            iir += 1
        self.hug.NH=len(self.hug.P)
        self.hug.units='units: T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg, Up km/s, Us km/s, cs km/s'

        if writefilename:
            print('Writing Hugoniot to file: ',writefilename)
            hugoniotfile = open(writefilename,"w")  
            hugoniotfile.writelines('  Hugoniot \n') 
            hugoniotfile.writelines('  Temperature,    Density,        Pressure,       Int. Energy,    Sp. Entropy,    Part. Vel.,     Shock Vel. \n') 
            hugoniotfile.writelines('  K,              g/cm3,          GPa,            MJ/kg,          MJ/K/kg,        km/s,           km/s\n') 
            for iih in range(0,self.hug.NH):
                hugoniotfile.write("%14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e\n" % (
                    self.hug.T[iih],self.hug.rho[iih],self.hug.P[iih],self.hug.U[iih],self.hug.S[iih],self.hug.up[iih],self.hug.us[iih]))
            hugoniotfile.close() 
#
    def loadaneos(self, aneosinfname=None, aneosoutfname=None, silent=False):
        """Function for reading in ANEOS INPUT and OUTPUT FILE DATA into EOS structure."""
        if aneosinfname is None:
            return 'Must provide input file name.'
        if aneosoutfname is None:
            return 'Must provide output file name.'
        # function to gather data from ANEOS input and output files
        # SESAME FILE HEADER INFORMATION MUST BE LOADED INTO THE EOS STRUCTURE BEFORE CALLING THIS FUNCTION
        #
        # READ IN ANEOS INPUT FILE
        aneosinputfile = open(aneosinfname,"r")  
        testin=aneosinputfile.readlines()   # read in the whole ascii file at once because this is fatser
        aneosinputfile.close()
        # gather EOS information from the ANEOS.OUTPUT file
        aneosoutputfile = open(aneosoutfname,"r")  
        testout=aneosoutputfile.readlines() # read everything in at once because this is faster
        aneosoutputfile.close()
        if silent == False:
            print('Done loading ANEOS files.')

        # THIS CODE PARSES THE ANEOS.OUTPUT FILE INTO ARRAYS FOR USE IN PLOTTING/USING THE EOS
        if silent == False:
            print('ANEOS WAS CALLED WITH THE FOLLOWING INPUT, LOADED FROM FILE ',aneosinfname)
        # Gather parameters for the gamma function while printing the ANEOS INPUT FILE
        aneoscount=1
        for i in np.arange(len(testin)):
            if testin[i].find('ANEOS') == 0:
                if aneoscount<9:
                    if silent == False:
                        print(' '+testin[i-3],testin[i-2],testin[i-1],testin[i])
                    aneoscount=aneoscount+1
                else:
                    if silent == False:
                        print(' '+testin[i])
                if testin[i].find('ANEOS2') == 0:
                    tmp=testin[i]
                    nelem=int(tmp[10:20])
                    #print('nelem=',nelem)
                    rho0=float(tmp[30:40])
                    #print('rho0=',rho0)
                    gamma0=float(tmp[70:80])
                    #print('gamma0=',gamma0)
                    theta0=float(tmp[80:90])
                if testin[i].find('ANEOS3') == 0:
                    tmp=testin[i]
                    C24=float(tmp[20:30])/3.
                    #print('C24=',C24)
                if testin[i].find('ANEOS5') == 0:
                    tmp=testin[i]
                    C60=float(tmp[60:70])
                    C61=float(tmp[70:80])
                    #print('C60=',C60)
                if testin[i].find('ANEOS7') == 0:
                    tmp=testin[i]
                    betagamma=float(tmp[70:80])

        # some checks
        if rho0 != self.R0REF:
            print('WARNING: rho0 does not match. STOPPING THIS NOTEBOOK.')
            assert(False) # just a way to stop the notebook

        # GUESS A BIG ARRAY SIZE FOR THE PHASE BOUNDARIES AND HUGONIOT IN ANEOS.OUTPUT
        # the melt curve, vapor curve and Hugoniot curves are not fixed length outputs
        nleninit=300
        meltcurve = 0

        if silent == False:
            print('READING DATA FROM ANEOS OUTPUT FILE ',aneosoutfname)

        # Read in data from the ANEOS.OUTPUT FILE
        imc = -1 # flag for no melt curve in the model
        for i in np.arange(len(testout)):
            if testout[i].find('  Data for ANEOS number') == 0:
                tmp = testout[i+2][0:50]
                eosname = tmp.strip()
            if testout[i] == '  TWO-PHASE BOUNDARIES\n':
                nvc = nleninit
                ivc = i
                vcarrtmp = np.zeros((nvc,12),dtype=float)
                flag=0
                j=0
                while flag == 0:
                    if testout[j+i+4].find(' anphas') == 0:
                        print(testout[j+i+4])
                        vcarrtmp[j,:]=vcarrtmp[j-1,:]
                        j=j+1
                    else:
                        tmp=str.replace(testout[j+i+4],'D','E')
                        tmp3 = tmp[0:157]
                        tmp4 = list(tmp3.split())
                        if (len(tmp4) >0) and (float(tmp4[3]) > 0) and (float(tmp4[4]) > 0): # stop if the pressures become negative on the vapor curve
                            tmp5 = np.asarray(tmp4)
                            vcarrtmp[j,:] = tmp5[:]
                            j=j+1
                        else:
                            flag=1
                vcarr = np.zeros((j,12),dtype=float)
                vcarr[:,:] = vcarrtmp[0:j,:]
            if testout[i] == ' LIQUID/SOLID PHASE CURVE\n':
                nmc = nleninit
                imc = i
                meltcurve=1
                mcarrtmp = np.zeros((nmc,11),dtype=float)
                flag=0
                j=0
                while flag == 0:
                    tmp  = str.replace(testout[j+i+5],'D','E')
                    tmp3 = tmp[0:132]
                    tmp4 = list(tmp3.split())
                    if len(tmp4) > 0:
                        tmp5 = np.asarray(tmp4)
                        mcarrtmp[j,:] = tmp5[:]
                        j=j+1
                    else:
                        flag=1
                mcarr = np.zeros((j,11),dtype=float)
                mcarr[:,:] = mcarrtmp[0:j,:]
            if testout[i] == '   HUGONIOT\n':
                nhc = nleninit
                ihc = i
                hcarrtmp = np.zeros((nhc,9),dtype=float)
                flag=0
                j=0
                while flag == 0:
                    tmp=str.replace(testout[j+i+5],'D','E')
                    tmp3 = tmp[0:109]
                    tmp4 = list(tmp3.split())
                    if len(tmp4) > 0:
                        tmp4[3]='0.0' # this column often gives problems with exponential notation so don't read it
                        tmp5 = np.asarray(tmp4)
                        hcarrtmp[j,:] = tmp5[:]
                        j=j+1
                    else:
                        flag=1
                hcarr = np.zeros((j,9),dtype=float)
                hcarr[:,:] = hcarrtmp[0:j,:]

        # UPDATE THE MAIN EOS STRUCTURE WITH GATHERED INFORMATION
        # Add variables needed to calculate the ANEOS gamma function
        self.gamma0  = gamma0
        self.theta0  = theta0
        self.C24     = C24
        self.C60     = C60
        self.C61     = C61
        self.beta    = betagamma
        #
        # ANEOS.OUTPUT UNITS ARE NOT THE SAME AS THE SESAME TABLE!
        # add the vapor curve to this EOS object extracted from the ANEOS.OUTPUT
        #  TWO-PHASE BOUNDARIES
        #       T         RHOLIQ        RHOVAP        PLIQ         PVAP        ELIQ         EVAP         SLIQ         SVAP        GLIQ         GVAP         PSILIQ      PSIVAP         NTY
        #       K         kg/m**3       kg/m**3       GPa          GPa         J/kg         J/kg        J/kg-K       J/kg-K       J/kg         J/kg
        tmp = vcarr.shape
        #put vapor curve information in nicely named structure
        self.vc.NT  = tmp[0]
        self.vc.T   = vcarr[:,0] # K
        self.vc.rl  = vcarr[:,1]/1.E3 # g/cm3
        self.vc.rv  = vcarr[:,2]/1.E3 # g/cm3
        self.vc.Pl  = vcarr[:,3] # GPa
        self.vc.Pv  = vcarr[:,4] # GPa
        self.vc.Ul  = vcarr[:,5]/1.E6 # MJ/kg
        self.vc.Uv  = vcarr[:,6]/1.E6 # MJ/kg
        self.vc.Sl  = vcarr[:,7]/1.E6 # MJ/K/kg
        self.vc.Sv  = vcarr[:,8]/1.E6 # MJ/K/kg
        self.vc.Gl  = vcarr[:,9]/1.E6 # MJ/kg
        self.vc.Gv  = vcarr[:,10]/1.E6 # MJ/kg
        self.vc.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg, G MJ/kg'
        # np.interp wants increasing x values
        self.onebar.Tvap = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.T)) # extract boiling point temperature at 1 bar, K
        self.onebar.Siv  = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.Sl)) # extract liquid sp. entropy at 1 bar, MJ/K/kg
        self.onebar.Scv  = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.Sv)) # extract vapor sp. entropy at 1 bar, MJ/K/kg
        self.onebar.rhoiv  = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.rl)) # extract liquid density at 1 bar, g/cm3
        self.onebar.rhocv  = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.rv)) # extract vapor density at 1 bar, g/cm3
        #
        # add the ANEOS Hugoniot to this EOS object extracted from the ANEOS.OUTPUT
        #      RHO          T           P          PC           E           S           V           U       RHO/RHOO  #IT  STATE
        #    kg/m**3        K          GPa        GPa          J/kg      J/kg-K       km/sec      km/sec
        #self.anhug.all = hcarr  # 2D array of Hugoniot variables
        tmp = hcarr.shape
        self.anhug.ND = tmp[0] # number of density points on the Hugoniot
        self.anhug.rho = hcarr[:,0]/1.E3 # g/cm3
        self.anhug.T   = hcarr[:,1] # K
        self.anhug.P   = hcarr[:,2] # GPa
        self.anhug.U   = hcarr[:,4]/1.E6 # MJ/kg
        self.anhug.S   = hcarr[:,5]/1.E6 # MJ/K/kg
        self.anhug.us  = hcarr[:,6] # km/s
        self.anhug.up  = hcarr[:,7] # km/s
        self.anhug.units = 'vars: rho g/cm3, T K, P GPa, U MJ/kg, S MJ/K/kg, Us km/s, Up km/s'
        #
        # Add melt curve to EOS objects if available
        # LIQUID/SOLID PHASE CURVE
        #       T         RLIQ       RSOLID      PLIQ       PSOLID      ELIQ        ESOLID       SLIQ       SOLID        GLIQ       GSOLID        #ITER
        #       K        kg/m**3     kg/m**3      GPa         GPa       J/kg         J/kg        J/kg-K     J/kg-K       J/kg        J/kg
        if meltcurve == 1:
            # put the melt curve information in nicely named structure
            tmp=mcarr.shape
            self.mc.NT  = tmp[0]
            self.mc.T   = mcarr[:,0] # K
            self.mc.rl  = mcarr[:,1]/1.E3 # g/cm3
            self.mc.rs  = mcarr[:,2]/1.E3 # g/cm3
            self.mc.Pl  = mcarr[:,3] # GPa
            self.mc.Ps  = mcarr[:,4] # GPa
            self.mc.Ul  = mcarr[:,5]/1.E6 # MJ/kg
            self.mc.Us  = mcarr[:,6]/1.E6 # MJ/kg
            self.mc.Sl  = mcarr[:,7]/1.E6 # MJ/K/kg
            self.mc.Ss  = mcarr[:,8]/1.E6 # MJ/K/kg
            self.mc.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg'
            # NOTE THAT TRIPLE POINT AND VAPOR CURVE SOLUTIONS DO NOT ALWAYS MATCH PERFECTLY AT THE TRIPLE POINT
            tmp = np.where(mcarr[:,3] > 0.)[0] # find the triple point first entry with positive pressure
            self.tp.T = mcarr[tmp[0],0] # K
            self.tp.P = np.interp(self.tp.T,np.flipud(self.vc.T),np.flipud(self.vc.Pv))   # this has trouble for forsterite; use the vapor size of the VC mcarr[:,3] # GPa
            self.tp.Sim  = mcarr[tmp[0],8]/1.E6 # extract solid sp. entropy at tp, MJ/K/kg
            self.tp.Scm  = mcarr[tmp[0],7]/1.E6 # extract liquid sp. entropy at tp, MJ/K/kg
            # np.interp wants x values increasing
            self.tp.Siv  = np.interp(self.tp.T,np.flipud(self.vc.T),np.flipud(self.vc.Sl)) # extract liquid sp. entropy at tp, MJ/K/kg
            self.tp.Scv  = np.interp(self.tp.T,np.flipud(self.vc.T),np.flipud(self.vc.Sv)) # extract vapor sp. entropy at tp, MJ/K/kg
            self.tp.rhol = mcarr[tmp[0],1]/1.E3 # extract liquid density at tp, g/cm3
            self.tp.rhos = mcarr[tmp[0],2]/1.E3 # extract solid density at tp, g/cm3
            self.tp.rhov = np.interp(self.tp.T,np.flipud(self.vc.T),np.flipud(self.vc.rv)) # extract vapor density at tp, g/cm3
            self.tp.units = 'T K, P GPa, S MJ/K/kg, rho g/cm3'
            # Extract melting point
            self.onebar.Tmelt = np.interp(1.E-4,self.mc.Pl,self.mc.T) # extract melting point temperature at 1 bar, K
            self.onebar.Sim   = np.interp(1.E-4,self.mc.Pl,self.mc.Ss) # extract liquid sp. entropy at 1 bar BP, MJ/K/kg
            self.onebar.Scm   = np.interp(1.E-4,self.mc.Pl,self.mc.Sl) # extract vapor sp. entropy at 1 bar BP, MJ/K/kg
            self.onebar.rhoim   = np.interp(1.E-4,self.mc.Pl[3::],self.mc.rs[3::]) # extract solid density at 1 bar MP, MJ/K/kg
            self.onebar.rhocm   = np.interp(1.E-4,self.mc.Pl[3::],self.mc.rl[3::]) # extract liquid density at 1 bar MP, MJ/K/kg
        # put the data for the critical point in the EOS structure for easy access
        self.cp.T   = vcarr[0,0] # K
        self.cp.rho = vcarr[0,1]/1.E3 # g/cm3
        self.cp.P   = vcarr[0,3] # GPa
        self.cp.U   = vcarr[0,5]/1.E6 # MJ/kg 
        self.cp.S   = vcarr[0,7]/1.E6 # MJ/K/kg
        self.cp.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg'
        #------------------------------------------------------------------

###########################################################################################
########## GADGET STYLE TABLES
###########################################################################################
#Phil J Carter Sept 5, 2019 gadget.py
class GadgetHeader:
    """Class for Gadget snapshot header."""
    def __init__(self, t=0, nfiles=1, ent=1):
        self.npart = np.zeros(6)
        self.mass = np.zeros(6)
        self.time = t
        self.redshift = 0
        self.flag_sfr = 0
        self.flagfeedbacktp = 0
        self.npartTotal = np.zeros(6)
        self.flag_cooling = 0
        self.num_files = nfiles
        self.BoxSize = 0
        self.Omega0 = 0
        self.OmegaLambda = 0
        self.HubbleParam = 1
        self.flag_stellarage = 0
        self.flag_metals = 0
        self.nallhw = np.zeros(6)
        self.flag_entr_ics = ent
#
class Snapshot:
    """Gadget snapshot class
    Includes header and gas particle data, with functions for reading and writing snapshots.
    load() -- load Gadget snapshot data
    remove() -- remove particle from snapshot
    write() -- save snapshot
    identify() -- determine material types
    calc_vap_frac() -- calculate vapour fractions of particles
    #GOH 01/15/2020
    -- fit midplane density profile
    -- fit midplane entropy profile
    -- fit midplane pressure profile
    -- fit midplane temperature profile
    -- fit midplane velocity profile
    -- fit midplane sound speed profile
    -- fit scale height for density
    -- fit scale height for entropy
    -- fit scale height for pressure
    -- fit scale height for temperature
    -- fit scale height for velocity profile
    -- fit scale height for sound speed profile
    """
    def __init__(self):
        self.header = GadgetHeader()
        self.N = 0
        self.pos = np.zeros(3)
        self.vel = np.zeros(3)
        self.id = 0
        self.m = 0
        self.S = 0
        self.rho = 0
        self.hsml = 0
        self.pot = 0
        self.P = 0
        self.T = 0
        self.U = 0
        self.cs = 0
        #self.accel = 0
        #self.dt = 0
        #self.vapfrac = 0
        self.J2Ma2 = 0
        self.ind_outer_mid_spl = 0
        self.pmidfit = 0
        self.rhomidfit = 0,0,0
#
    def load(self, fname, thermo=False):
        f = open(fname, 'rb')
        struct.unpack('i', f.read(4))
        #HEADER
        self.header.npart = np.array(struct.unpack('iiiiii', f.read(24)))
        self.header.mass = np.array(struct.unpack('dddddd', f.read(48)))
        (self.header.time, self.header.redshift, self.header.flag_sfr,
         self.header.flag_feedbacktp) = struct.unpack('ddii', f.read(24))
        self.header.npartTotal = np.array(struct.unpack('iiiiii', f.read(24)))
        (self.header.flag_cooling, self.header.num_files, self.header.Boxsize,
         self.header.Omega0, self.header.OmegaLambda, self.header.HubbleParam,
         self.header.flag_stellarage,
         self.flag_metals) = struct.unpack('iiddddii', f.read(48))
        #print(self.header.Boxsize,self.header.flag_stellarage,self.flag_metals)
        self.header.nallhw = np.array(struct.unpack('iiiiii', f.read(24)))
        self.header.flag_entr_ics = struct.unpack('i', f.read(4))
        struct.unpack('60x', f.read(60))
        struct.unpack('i', f.read(4))
        if self.header.num_files != 1:
            print("WARNING! Number of files:", self.header.num_files,
                  ", not currently supported.\n")
        self.N = self.header.npart[0]
        count=str(self.N)
        count3=str(3*self.N)
        #PARTICLE DATA
        struct.unpack('i', f.read(4))
        self.pos = struct.unpack(count3 + 'f', f.read(3*self.N*4))
        struct.unpack('i', f.read(4))
        struct.unpack('i', f.read(4))
        self.vel = struct.unpack(count3 + 'f', f.read(3*self.N*4))
        struct.unpack('i', f.read(4))
        struct.unpack('i', f.read(4))
        self.id = np.array(struct.unpack(count + 'i', f.read(self.N*4)))
        struct.unpack('i', f.read(4))
        struct.unpack('i', f.read(4))
        self.m = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
        struct.unpack('i', f.read(4))
        struct.unpack('i', f.read(4))
        self.S = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
        struct.unpack('i', f.read(4))
        struct.unpack('i', f.read(4))
        self.rho = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
        struct.unpack('i', f.read(4))
        struct.unpack('i', f.read(4))
        self.hsml = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
        struct.unpack('i', f.read(4))
        struct.unpack('i', f.read(4))
        self.pot = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
        struct.unpack('i', f.read(4))
        if thermo:
            struct.unpack('i', f.read(4))
            self.P = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
            struct.unpack('i', f.read(4))
            struct.unpack('i', f.read(4))
            self.T = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
            struct.unpack('i', f.read(4))
            if len(f.read(4)) == 4:
                self.U = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
                struct.unpack('i', f.read(4))
            if len(f.read(4)) == 4:
                self.cs = np.array(struct.unpack(count + 'f', f.read(self.N*4)))
                struct.unpack('i', f.read(4))
        f.close()
        #REARRANGE
        self.pos = np.array(self.pos).reshape((self.N, 3))*(1e-2) #m
        self.x = self.pos.T[0]
        self.y = self.pos.T[1]
        self.z = self.pos.T[2]
        self.vel = np.array(self.vel).reshape((self.N, 3))*(1e-2) #m/s
        self.vx = self.vel.T[0]
        self.vy = self.vel.T[1]
        self.vz = self.vel.T[2]
        print("Read %d" % self.N, "particles from %s" % fname)
        #CALCULATE CENTER OF MASS
        N=25
        temp=np.argsort(self.pot)
        xcm=np.mean(self.x[temp[0:N]])
        ycm=np.mean(self.y[temp[0:N]])
        zcm=np.mean(self.z[temp[0:N]])
        vxcm=np.mean(self.vx[temp[0:N]])
        vycm=np.mean(self.vy[temp[0:N]])
        vzcm=np.mean(self.vz[temp[0:N]])
        #MOVE ONTO CENTRAL FRAME
        self.x=self.x-xcm
        self.y=self.y-ycm
        self.z=self.z-zcm
        self.vx=self.vx-vxcm
        self.vy=self.vy-vycm
        self.vz=self.vz-vzcm
        #CALCULATE BOUND MASS
        self.m = self.m*(1e-3) #kg
        #bndm=self.m[temp[0]]
        #G=6.67408E-11 #mks
        #oldm=bndm/10.
        #tol=1E-5
        #while (np.abs(oldm-bndm)>tol):
        #    oldm=bndm
        #    v2=np.add(np.add(np.power(self.vx,2.0),np.power(self.vy,2.0))np.power(self.vz,2.0))
        #    r=np.sqrt(np.add(np.add(np.power(self.x,2.0),np.power(self.y,2.0))np.power(self.z,2.0)))
        #    KE=0.5*np.multiply(self.m,v2)
        #    PE=-G*bndm*np.divide(self.m,r)
        #    bndm=np.sum(self.m[np.where((KE+PE)<0.)[0]])
        #CONVERT REST OF UNITS TO MKS
        self.rho = self.rho*(1e3) #kg/m3
        self.P = self.P*1e9 #Pa
        self.S = self.S*(1e-4) #J/K/kg
        self.pot = self.pot*(1e-4) #J/kg
        self.U = self.U*(1e-4) #J/kg
        self.cs = self.cs*(1e-2) #m/s
        self.rxy = np.add(np.power(self.x, 2), np.power(self.y, 2)) #m2
        radius2 = np.add(self.rxy,np.power(self.z,2)) #m2
        self.rxy = np.sqrt(self.rxy) #m
        self.J2Ma2 = -np.sum(0.5*np.multiply(self.m,radius2)*(3.0*np.divide(np.power(self.z,2),radius2) - 1.0)) #kg m2
        #print("Centered bound mass.\n")
#
    def indices(self,zmid,zmax,rxymin,rxymax,rxymida,rxymidb):
        #DETERMINE OUTER REGION PARTICLES (truncated at rxymin and rxymax)
        self.ind_outer=np.where((self.rxy >= rxymin) & (self.rxy <= rxymax) & (np.abs(self.z) <= zmax))
        self.ind_outer_1=np.where((self.rxy >= rxymin) & (self.rxy < rxymida) & (np.abs(self.z) <= zmax))
        self.ind_outer_2=np.where((self.rxy > rxymidb) & (self.rxy <= rxymax) & (np.abs(self.z) <= zmax))
        self.ind_outer_S=np.where(self.rxy >= rxymin)

        #DETERMINE MIDPLANE OUTER REGION PARTICLES
        self.ind_outer_mid=np.where((self.rxy >= rxymida) & (np.abs(self.z) <= zmid) & (self.rxy <= rxymidb))
        self.ind_outer_mid_spl = np.where((np.abs(self.z) <= zmid) & (self.rxy <= rxymax) & (self.rxy >= rxymin))
        self.ind_outer_mid_lsq=np.where((np.abs(self.z) <= zmid) & (self.rxy >= 9.4e6))
        self.ind_mid=np.where(np.abs(self.z) <= zmid)
        
    def fit_Pmid(self,knots,extra=None):
        #DETERMINE SPLINE FIT TO MIDPLANE PRESSURE CURVE
        ind_outer_mid_spl=np.where((np.abs(SNAP.z) <= zmid) & (SNAP.rxy <= rxymax) & (SNAP.rxy >= rxymin))
        indsort=np.argsort(SNAP.rxy[ind_outer_mid_spl])
        SPHrxyMm = SNAP.rxy[ind_outer_mid_spl][indsort]/1e6
        SPHplog = np.log10(SNAP.P[ind_outer_mid_spl][indsort])
        pknots=[*knots]
        pLSQUSPL = LSQUnivariateSpline(SPHrxyMm, SPHplog, t=pknots, k=3)
        if extra:
            print('knots for midplane pressure curve are rxy = {}'.format(pLSQUSPL.get_knots()))
            print('coefficients for midplane pressure curve are {}'.format(pLSQUSPL.get_coeffs()))
        return pLSQUSPL
    def fit_rhomid(self,extra=None):
        #DETERMINE LEAST-SQUARES FIT TO RESIDUAL OF MIDPLANE RHO S-CURVE 1
        params_guess=np.ones(4)
        res_lsq = least_squares(resfunc, params_guess, loss='soft_l1', f_scale=0.001, 
                                args=(np.log10(SNAP.rxy[ind_outer_2]/1e6), np.log10(SNAP.rho[ind_outer_2])))

        #DETERMINE LEAST-SQUARES FIT TO RESIDUAL OF MIDPLANE RHO S-CURVE 2
        params_guess_spl=np.array([150.,1.4,16.,-5.7])
        res_lsq_spl = least_squares(resfuncspl, params_guess_spl, loss='soft_l1', f_scale=0.001, 
                                    args=(np.log10(SNAP.rxy[ind_outer_mid]/1e6), np.log10(SNAP.rho[ind_outer_mid])))

        #DETERMINE LEAST-SQUARES FIT TO RESIDUAL OF MIDPLANE RHO LINE
        params_guess_lin=np.ones(2)
        res_lsq_lin = least_squares(resfunclin, params_guess_lin, loss='soft_l1', f_scale=0.001,
                                    args=(np.log10(SNAP.rxy[ind_outer_1]/1e6), np.log10(SNAP.rho[ind_outer_1])))

        if extra:
            print('Least Squares Fit to Midplane Density - S-curve \n')
            print(res_lsq)
            print('\n Least Squares Fit to Midplane Density - Spline \n')
            print(res_lsq_spl)
            print('\n Least Squares Fit to Midplane Density - Linear \n')
            print(res_lsq_lin)
            print('\n Params for midplane density:\n fit 0 {}\n fit 1 {}\n fit 2 {}\n Linear interpolation points are (x1_lim, y1_lim) = ({}, {}) and (x2_lim, y2_lim) = ({}, {})'.format(res_lsq_lin.x,res_lsq_spl.x,res_lsq.x,x1int,y1int,x2int,y2int))
            
        self.fitrhomid = res_lsq_lin.x,res_lsq_spl.x,res_lsq.x
#        
    def fit_Tmid(self,extra=None):
        params_guess_T=np.asarray([4.e12,-1.66,2.5])
        res_lsq_pow = least_squares(resfuncpow, params_guess_T, ftol=1e-10, xtol=1e-11, loss='soft_l1',
                                    args=(SNAP.rxy[ind_outer_mid_lsq], SNAP.T[ind_outer_mid_lsq]/1.e3))
        if extra:
            print('\n Least Squares Fit to Midplane Temperature - Power Law \n')
            print(res_lsq_pow)
#        
    def fit_smid(self,extra=None):
        params_guess_S = np.ones(5)
        res_lsq_lp = least_squares(resfunclinpiece, params_guess_S, ftol=1e-8, xtol=1e-8, loss='soft_l1', 
                                   args=(SNAP.rxy[ind_outer_mid_spl]/1e6, SNAP.S[ind_outer_mid_spl]))
        if extra:
            print('\n Least Squares Fit to Midplane Entropy - Linear Piecewise \n')
            print(res_lsq_lp)
#            
    def fit_zs_rho(self,extra=None):
        #SCALE HEIGHT FIT
        #bin by rxy, fit each bin's rho(z) profile and find z where rho/rho_mid=1/e
        ind_outer_offmid = np.where((SNAP.rxy >= 7.e6) & (np.abs(SNAP.z) > 1.e6))
        bins = np.arange(7.e6,np.amax(SNAP.rxy[ind_outer_offmid])+1.e6,1.e6)
        #bins_S = np.arange(7.e6,np.amax(SNAP.rxy[ind_outer_S])+1.e6,1.e6)
        ind_bin = np.digitize(SNAP.rxy[ind_outer_offmid],bins,right=False)
        #ind_bin_S = np.digitize(SNAP.rxy[ind_outer_S],bins_S,right=False)
        bins=(bins/1e6)-0.5 #convert to Mm
        #bins_S=(bins_S/1e6)-0.5
        params_guess_exp = 1.
        #def resfuncpieceS(params,x,y,z):
            #x is rxy, y is z, z is S
        #    f1 = params[0]
        #    f2 = params[1]
        #    f3 = lambda y: params[2]*y**2 + params[3]*y + params[4]
        #    return np.select([(y>params[5]),(y<=params[5])*(y>=params[6]),(y<params[6])*(x<10.)], [f1,f2,f3]) - z
        #params_guess_Spiece = np.asarray([4500.,8000.,1.,1.,4000.,15.,1.])
        res_lsq_exp = []
        #res_lsq_Spiece = []
        for i in range(len(bins)):
            ind_rxy = np.where(ind_bin == i)
            SNAP_rhodiv_offmid = np.log(SNAP.rho[ind_outer_offmid][ind_rxy]*(
                10**(-piece(np.log10(SNAP.rxy[ind_outer_offmid][ind_rxy]/1.e6),res_lsq_lin.x,res_lsq_spl.x,res_lsq.x))))
            reslsqexp = least_squares(resfuncexp, params_guess_exp, bounds=(1,np.inf), loss='soft_l1', f_scale=0.001, 
                                      args=(np.abs(SNAP.z[ind_outer_offmid][ind_rxy]/1e6),SNAP_rhodiv_offmid))
            if reslsqexp.active_mask[0] == -1:
                res_lsq_exp.append(np.nan)
            else:
                res_lsq_exp.append(reslsqexp.x[0])
        #for i in range(len(bins_S)):
        #    ind_rxy_S = np.where(ind_bin_S == i)
        #    print(ind_rxy_S.active_mask)
        #    if ind_rxy_S.active_mask == -1:
        #        res_lsq_Spiece.append(np.nan)
        #    else:
        #        reslsqSpiece = least_squares(resfuncpieceS, params_guess_Spiece, loss='soft_l1', f_scale=0.001, args=(SNAP.rxy[ind_outer_S][ind_rxy_S]/1.e6,np.abs(SNAP.z[ind_outer_S][ind_rxy_S])/1e6,SNAP.S[ind_outer_S][ind_rxy_S]))
        #        res_lsq_Spiece.append(reslsqSpiece.x)

        res_lsq_exp = np.asarray(res_lsq_exp) #Mm
        #res_lsq_Spiece = np.asarray(res_lsq_Spiece)
        print('\n Binned Rxy Scale Height Fits')
        print(res_lsq_exp)
        #print(res_lsq_Spiece)

        #MASKING SCALE HEIGHT FITS FOR NAN'S AND Z_S > 100 Mm
        res_lsq_exp_mask = np.ma.fix_invalid(res_lsq_exp)
        res_lsq_exp_compress = res_lsq_exp_mask.compressed()
        bins_mask = np.ma.masked_array(bins, mask=res_lsq_exp_mask.mask)
        bins_compress = bins_mask.compressed()
        res_lsq_exp_compress_mask = np.ma.masked_greater(res_lsq_exp_compress,100.)
        res_lsq_exp_compress2 = res_lsq_exp_compress_mask.compressed()
        bins_compress_mask = np.ma.masked_array(bins_compress, mask=res_lsq_exp_compress_mask.mask)
        bins_compress2 = bins_compress_mask.compressed()
        print('\n Masked Rxy Scale Heights')
        print(zip(bins_compress2,res_lsq_exp_compress2))
        knots=[25.5,30.5,32.5,39.5,45.5,57.5]
        LSQUSPL=LSQUnivariateSpline(bins_compress2, res_lsq_exp_compress2, t=knots, k=3)
        if extra:
            bknot = LSQUSPL.get_knots()
            bcoef = LSQUSPL.get_coeffs()
            print('\n LSQ Univariate Spline Fit to Scale Heights \n')
            print('knots are rxy = {} (Mm)'.format(bknot))
            print('coefficients are {}'.format(bcoef))
#
#FUNCTIONAL FORMS
def resfunc(params,x,y):
    return params[0]*x*(1. + np.exp(-x*x*params[1]))**(-params[2]) + params[3]*x - y
#
def resfuncspl(params,x,y):
    return params[0]/(x*np.exp(x*params[1]) + params[2]) + params[3] - y
#
def resfunclin(params,x,y):
    return params[0]*x + params[1] - y
#
def resfuncpow(params,x,y):
    return params[0]*(x**params[1]) + params[2] - y
#
def resfunclinpiece(params,x,y):
    return np.piecewise(x, [x < 10., ((x >= 10.)&(x <= 35.)), x > 35.],
                        [lambda x: params[0]*x + params[1], params[2],
                         lambda x: params[3]*x + params[4]]) - y
#
def resfuncexp(params,x,y):
    return x/params[0] + y
#
def piece(x,rxymidb,params0,params1,params2,extra=None):
    x2lim=np.log10(rxymidb/1e6 + 1.)
    x1lim=np.log10(rxymidb/1e6 - 1.)
    x0lim=1. #log10 of 10 Mm
    x2 = x[np.asarray(x >= x2lim)] #make sure none of these are empty; may throw error
    x1 = x[np.asarray((x <= x1lim) & (x >= x0lim))]
    x0 = x[np.asarray(x < x0lim)]
    y2=resfunc(params2, x2, np.zeros(len(x2)))
    y1=resfuncspl(params1, x1, np.zeros(len(x1)))
    y0=resfunclin(params0, x0, np.zeros(len(x0)))
    xmid = x[np.asarray((x < x2lim) & (x > x1lim))]
    y2lim = resfunc(params2, x2lim, 0.)
    y1lim = resfuncspl(params1, x1lim, 0.)
    ymid = np.interp(xmid,np.asarray([x1lim,x2lim]),np.asarray([y1lim,y2lim]))
    ytemp01=np.concatenate((y0,y1))
    ytemp1mid=np.concatenate((ytemp01,ymid))
    y=np.concatenate((ytemp1mid,y2))
    xtemp01=np.concatenate((x0,x1))
    xtemp1mid=np.concatenate((xtemp01,xmid))
    x=np.concatenate((xtemp1mid,x2))
    if extra:
        return y,10**x0,y0,10**x1,y1,10**xmid,ymid,10**x2,y2,x1lim,y1lim,x2lim,y2lim
    else:
        return y

##################################################################################
### BEGIN synfits.py RUN ###
##################################################################################

#LOAD IN GADGET-2 SNAPSHOTS
SNAP_Canup=Snapshot()
SNAP_Canup.load('TE_Example01_Cool05_snapshot_4096_long',thermo=True) #Canup 2012 style giant impact
SNAP_CukStewart=Snapshot()
SNAP_CukStewart.load('TE_Example03_Cool01_snapshot_10500_long',thermo=True) #Cuk & Stewart 2012 style giant impact
SNAP_Quintana=Snapshot()
SNAP_Quintana.load('TE_Example07_CoolB01_snapshot_7200_long',thermo=True) #Quintana style giant impact
zmax=100.e6 #m
zmid=1.e6 #m
rxymin=7.e6 #m
rxymax=100.e6 #m
rxymida=10.e6 #m
rxymidb=40.e6 #m
SNAP_Canup.indices(zmid,zmax,rxymin,rxymax,rxymida,rxymidb)
SNAP_CukStewart.indices(zmid,zmax,rxymin,rxymax,rxymida,rxymidb)
SNAP_Quintana.indices(zmid,zmax,rxymin,rxymax,rxymida,rxymidb)
#SNAP_Canup.fit_Pmid()
import sys
sys.exit()

#DETERMINE STEP FOR HYDROSTATIC EQUIL CALC
step=100.

#DETERMINE SPLINE FIT TO MIDPLANE PRESSURE CURVE
SNAP.fit_Pmid([10.,35.,55.]) #knots are in Mm
ptemplsq = pLSQUSPL(rxyarr)

#CALCULATE RESIDUAL LEAST SQUARES FIT TO MIDPLANE DENSITY CURVE
SNAP.fit_rhomid()
testy=piece(np.log10(rxyarr),rxymidb,*self.rhomidfit)

#DETERMINE LEAST-SQUARES FIT TO RESIDUAL OF MIDPLANE TEMPERATURE CURVE
SNAP.fit_Tmid()
rxyarr_T=np.linspace(9.4,100.,1000)
y_zero_T=np.zeros(len(rxyarr_T))
Tfit = resfuncpow(self.res_lsq_pow.x,rxyarr_T*1e6,y_zero_T)

#MIDPLANE VELOCITIES
SNAPvxy = np.sqrt(SNAP.vx[ind_outer_mid_spl]**2 + SNAP.vy[ind_outer_mid_spl]**2) #Mm/s?
vkep = np.sqrt((6.67408e-11)*np.sum(SNAP.m)/(rxyarr*1e6)) #m/s

#LEAST-SQUARES FIT TO RESIDUAL OF MIDPLANE ENTROPY CURVE


#CALCULATING SPLINE FIT TO SCALE HEIGHTS
temp4=LSQUSPL(rxyarr) #density

#LOAD IN EOS -- NEED TO PLOT AGAINST THERMAL PARAMETERS
NewEOS=extEOStable()
NewEOS.loadextsesame('NEW-SESAME-EXT.TXT')
NewEOS.loadstdsesame('NEW-SESAME-STD.TXT')
NewEOS.MODELNAME='GADGET2 Forsterite'
NewEOS.MDQ=np.zeros((NewEOS.NT,NewEOS.ND))
NewEOS.MATID=1.0
NewEOS.DATE=190802.
NewEOS.VERSION=1.0
NewEOS.FMN=70.
NewEOS.FMW=140.691
NewEOS.R0REF=3.32
NewEOS.K0REF=1.446E12
NewEOS.T0REF=298.
NewEOS.P0REF=1.E6
NewEOS.loadaneos(aneosinfname='ANEOS.INPUT',aneosoutfname='ANEOS.OUTPUT')
#isen=isentrope_class()
#isen.loadisen(s_isen,NewEOS)

want_hydro=None
if want_hydro:
    #INITIALIZE HYDROSTAT CLASS
    rxy=np.linspace(rxymin,rxymax,100)
    z=np.linspace(0.,zmax,int((zmax-.1)//step)+1)
    hyd=hydrostat(1,rxy,z) #first argument is grid flag: 1 for arbitrary reg. grid, 0 for GADGET positions

    #LOAD PRESSURE FIT COEFFICIENTS, ISO-ENTROPY, J2MA2, TOTAL MASS INTO HYDROSTAT CLASS
    s_isen=7.866e-3 #MJ/K/kg
    #put in entropy fit into hydrostatic calculations
    hyd.loadhyd(pLSQUSPL,s_isen,SNAP.J2Ma2,np.sum(SNAP.m))

    #CALCULATE HYDROSTATIC EQUIL FOR GRID
    hyd.calc_hydrostat_equil(isen)

    #CALC DIFF IN PRESSURE BETWEEN SNAPSHOT & SEMI-ANALYTIAL CALC
    pdiff,pgrid=hyd.calc_pdiff(SNAP.rxy[ind_outer],np.abs(SNAP.z[ind_outer]),SNAP.P[ind_outer]) #bar
        #second argument is for press grid (bar) computed by griddata
    A=hyd.RXY/1e6
    B=hyd.Z/1e6
    C=hyd.P/1e5
    D=pgrid-(hyd.P/1e5)

#PLOTS
#plot midplane fits
plot_mid=None
if plot_mid:
    rxyarr1 = np.linspace(7.,9.999999,50)
    y_zero_S1 = np.zeros(len(rxyarr1))
    Sfit1 = resfunclinpiece(res_lsq_lp.x,rxyarr1,y_zero_S1)
    rxyarr2 = np.linspace(10.,35.,100)
    y_zero_S2 = np.zeros(len(rxyarr2))
    Sfit2 = resfunclinpiece(res_lsq_lp.x,rxyarr2,y_zero_S2)
    rxyarr3 = np.linspace(35.01,100.,200)
    y_zero_S3 = np.zeros(len(rxyarr3))
    Sfit3 = resfunclinpiece(res_lsq_lp.x,rxyarr3,y_zero_S3)
    
    fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(6,9))
    ax[0,0].plot(SNAP.x/1e6,SNAP.z/1e6,'.',color='grey',markersize=1)
    ax[0,0].plot(SNAP.x[ind_outer_mid_spl]/1e6,SNAP.z[ind_outer_mid_spl]/1e6,'k.',markersize=1)
    ax[0,0].set_xlabel('x (Mm)')
    ax[0,0].set_ylabel('z (Mm)')
    ax[0,0].set_aspect(1)
    ax[0,1].plot(SNAP.rxy[ind_outer_mid_spl]/1e6,SNAP.S[ind_outer_mid_spl],'k.',markersize=1)
    ax[0,1].plot(rxyarr1,Sfit1,'r',markersize=1,label='piece-wise 0')
    ax[0,1].plot(rxyarr2,Sfit2,'b',markersize=1,label='piece-wise 1')
    ax[0,1].plot(rxyarr3,Sfit3,'g',markersize=1,label='piece-wise 2')
    ax[0,1].set_xlabel('r$_{xy}$ (Mm)')
    ax[0,1].set_ylabel('outer region midplane\n specific entropy (J K$^{-1}$ kg$^{-1}$)',fontsize=8)
    ax[0,1].legend(loc=4,fontsize='small')
    ax[1,0].plot(SNAP.rxy[ind_outer_mid_spl]/1e6,np.log10(SNAP.P[ind_outer_mid_spl]),'k.',markersize=1)
    ax[1,0].plot(rxyarr,ptemplsq,'.',markersize=1)
    ax[1,0].set_xlabel('r$_{xy}$ (Mm)')
    ax[1,0].set_ylabel('outer region midplane\n pressure (log$_{10}$ Pa)',fontsize=8)
    ax[1,1].plot(SNAP.rxy[ind_outer_mid_spl]/1e6,np.log10(SNAP.rho[ind_outer_mid_spl]),'k.',markersize=1)
    ax[1,1].plot(testx0,testy0,'r',markersize=1,label='piece-wise 0')
    ax[1,1].plot(testx1,testy1,'b',markersize=1,label='piece-wise 1')
    ax[1,1].plot(testxmid,testymid,markersize=1,color='orange',label='piece-wise mid')
    ax[1,1].plot(testx2,testy2,'g',markersize=1,label='piece-wise 2')
    ax[1,1].set_xlabel('r$_{xy}$ (Mm)')
    ax[1,1].set_ylabel('outer region midplane\n density (log$_{10}$ kg/m$^3$)',fontsize=8)
    ax[1,1].legend(fontsize='small')
    ax[2,0].plot(SNAP.rxy[ind_outer_mid_spl]/1e6,SNAP.T[ind_outer_mid_spl],'k.',markersize=1)
    ax[2,0].plot(rxyarr_T,Tfit*1e3,'.',markersize=1)
    ax[2,0].set_xlabel('r$_{xy}$ (Mm)')
    ax[2,0].set_ylabel('outer region midplane\n temperature (K)',fontsize=8)
    ax[2,1].plot(rxyarr,vkep,'--',color='purple',markersize=1)
    ax[2,1].plot(SNAP.rxy[ind_outer_mid_spl]/1e6,SNAPvxy,'k.',markersize=1)
    ax[2,1].set_xlabel('r$_{xy}$ (Mm)')
    ax[2,1].set_ylabel('outer region midplane\n in-plane velocity v$_{xy}$ (m/s)',fontsize=8)
    fig.tight_layout()
    plt.savefig('Ex01Cool05Snap4096midfits.pdf',dpi=300,bbox_inches='tight')
    plt.show()
    plt.close()

#plot scale height of density
plot_zs=None
if plot_zs:
    #GRIDDATA OFF-THE-MIDPLANE DENSITY
    xi=SNAP.rxy[ind_outer_S]/1e6
    yi=np.abs(SNAP.z[ind_outer_S])/1e6
    test0,test1=np.meshgrid(np.linspace(7.,125.,100),np.linspace(0.,100.,100))
    test2=griddata((xi,yi),np.log10(SNAP.rho[ind_outer_S]),(test0,test1),method='linear')
    
    plt.figure()
    plt.pcolormesh(test0,test1,test2,cmap='viridis_r',vmin=-4,vmax=3)
    plt.colorbar(label='log density (kg/m$^3$)')
    plt.plot(bins,res_lsq_exp,'ko',markersize=1,label='z$_s$ bin')
    plt.plot(temp0,temp4,'m',markersize=1,label='cubic spline fit')
    plt.ylim([0,100])
    plt.xlim([7,125])
    plt.xlabel('r$_{xy}$ (Mm)')
    plt.ylabel('z (Mm)')
    plt.legend(loc=2,fontsize='xx-small')
    plt.savefig('Ex01Cool05Snap4096_rhofit_zscalefit.pdf',dpi=300,bbox_inches='tight')
    plt.show()
    plt.close()
    
#plot error analytic vs SPH pressures
plot_diff=None
if plot_diff:
    from matplotlib.colors import Normalize,SymLogNorm
    cmap1=plt.get_cmap('viridis_r')
    norm1=SymLogNorm(linthresh=1.,vmin=0.,vmax=1.e5)
    extent=[np.amin(A),np.amax(A),np.amin(B),np.amax(B)]
    """plt.figure()
    plt.imshow(pgrid.T,cmap=cmap1,norm=norm1,aspect='auto',origin='lower',interpolation='none',extent=extent)
    plt.colorbar(label='pressure (log bar)')
    plt.xlabel('r$_{xy}$ (Mm)')
    plt.ylabel('z (Mm)')
    plt.xlim([7,80])
    plt.ylim([0,40])
    plt.savefig('Ex01Cool05snap4096pressure.pdf',dpi=300,bbox_inches='tight')
    plt.show()
    plt.close()

    plt.figure()
    plt.imshow(C.T,cmap=cmap1,norm=norm1,aspect='auto',origin='lower',interpolation='none',extent=extent)
    plt.colorbar(label='pressure (log bar)')
    plt.xlabel('r$_{xy}$ (Mm)')
    plt.ylabel('z (Mm)')
    plt.xlim([7,80])
    plt.ylim([0,40])
    plt.savefig('Ex01Cool05snap4096hydrostatp.pdf',dpi=300,bbox_inches='tight')
    plt.show()
    plt.close()"""
    
    plt.figure()
    cmap2=plt.get_cmap('seismic')
    pdiffmax=np.amax(D)
    pdiffmin=np.amin(D)
    if pdiffmax < -pdiffmin:
        pdifflim = -pdiffmin
    elif pdiffmax >= -pdiffmin:
        pdifflim = pdiffmax
    norm2=SymLogNorm(linthresh=10.,vmin=-pdifflim,vmax=pdifflim)
    plt.imshow(D.T,cmap=cmap2,norm=norm2,aspect='auto',origin='lower',interpolation='none',extent=extent)
    plt.colorbar(label='diff. in press. between SPH and hyd. equil. (bar)')
    plt.xlabel('r$_{xy}$ (Mm)')
    plt.ylabel('z (Mm)')
    plt.xlim([7,80])
    plt.ylim([0,40])
    plt.savefig('Ex01Cool05snap4096hydrostatdiff.pdf',dpi=300,bbox_inches='tight')
    plt.show()
    plt.close()
    
####################################################################################
### END synfits.py FILE ###
####################################################################################
