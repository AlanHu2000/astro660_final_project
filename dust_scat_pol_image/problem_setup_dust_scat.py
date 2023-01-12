#
# Import NumPy for array handling
#
import numpy as np
import math
#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ss  = 5.6703e-5      # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
kk  = 1.3807e-16     # Bolzmann's constant     [erg/K]
mp  = 1.6726e-24     # Mass of proton          [g]
mH  = 1.66053e-24    # Mass of atomic hydrogen [g] CF
GG  = 6.67408e-08    # Gravitational constant  [cm^3/g/s^2]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
#
# Monte Carlo parameters
#
nphot    = 10000000
nphot_scat = 5000000
#
# Grid parameters
#
nx       = 100
ny       = 60
nz       = 10
#
# Model parameters
#
rin      = 0.25*au  ## change to fit HL tau
rout     = 150*au   ## change to fit HL tau
rho0     = 1e-16 * 5 *10/7.0905721352689118 * 10/7.4419820333705298
#prho     = -2.e0
H100min  = 1*au    ## change to fit HL tau
beta     = 1.15     ## change to fit HL tau
rc       = 78.9*au  ## change to fit HL tau
gamma    = -0.22    ## change to fit HL tau
amax     = 100      ## maximum grain size (mum)
#hpr      = 0.1e0
#
# Star parameters
#
mstar    = 1.7*ms
rstar    = 7*rs
tstar    = 4000     ## change to fit HL tau
pstar    = [0.,0.,0.]

#
# Vertical grid parameters (theta-grid in spherical coordinates)
#
zrmax    = 0.5
thetaup  = np.pi*0.5 - zrmax


#
# Make the coordinates
#
xi       = rin * (rout/rin)**(np.linspace(0.,nx,nx+1)/(nx-1.0))
yi       = np.linspace(thetaup,np.pi-thetaup,2*ny+1)
zi       = np.linspace(0.,math.pi*2,nz+1)
xc       = 0.5e0 * ( xi[:-1] + xi[1:] )
yc       = 0.5e0 * ( yi[:-1] + yi[1:] )
zc       = 0.5e0 * ( zi[:-1] + zi[1:] )


#
# Make the dust density model
#
rr,tt    = np.meshgrid(xc,yc,indexing='ij')
zzr      = math.pi/2 - tt

hpr      = H100min*((amax/0.01)**-0.1)*((rr/(100*au))**beta)
 
rhod     = rho0 * ((rr/rc)**(-1*gamma))*np.exp(-1*(rr/rc)**(2-gamma))
rhod     = rhod * np.exp(-0.50*(zzr/hpr)**2)

for i in range(rr.shape[0]):
    for j in range(rr.shape[1]):
        rr_value = rr[i][j]
        if rr_value >= 7.2*au and rr_value <= 20.2*au:
            rhod[i][j] = rhod[i][j]/18
        elif rr_value >= 26.8*au and rr_value <= 37.8*au:
            rhod[i][j] = rhod[i][j]/16
        elif rr_value >= 38.7*au and rr_value <= 45.3*au:
            rhod[i][j] = rhod[i][j]/6.9
        elif rr_value >= 47.75*au and rr_value <= 52.25*au:
            rhod[i][j] = rhod[i][j]/3.8
        elif rr_value >= 58.2*au and rr_value <= 70.2*au:
            rhod[i][j] = rhod[i][j]/8
        elif rr_value >= 69.65*au and rr_value <= 77.75*au:
            rhod[i][j] = rhod[i][j]/12
        elif rr_value >= 86.05*au and rr_value <= 95.95*au:
            rhod[i][j] = rhod[i][j]/11
        else:
            pass

rhod_extend = np.repeat(rhod[:, :, np.newaxis], nz, axis=2)        

Tc       = 70 * (rr/(10*au))**(-0.5)
Tc_extend = np.repeat(Tc[:, :, np.newaxis], nz, axis=2)


#
# Write the wavelength_micron.inp file
#
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
nlam     = lam.size


#
# Write the wavelength file
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    for value in lam:
        f.write('%13.6e\n'%(value))
#
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    for value in lam:
        f.write('%13.6e\n'%(value))
    f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('1 1 1\n')                   # Include x,y,z coordinate
    f.write('%d %d %d\n'%(nx,2*ny,nz))     # Size of grid
    for value in xi:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in yi:
        f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
    for value in zi:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*2*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rhod_extend.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')

with open('dust_temperature.dat','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*2*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = Tc_extend.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')


#
# Dust opacity control file
#
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('10               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('new_100        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')



#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('nphot_scat = %d\n'%(nphot_scat))
    f.write('scattering_mode_max = 5\n')   # Put this to 1 for isotropic scattering
    f.write('setthreads = 8\n')
    f.write('istar_sphere = 1\n')

    
with open('disk_properties.dat','w+') as f:
    for cR in xc:
        omega = np.sqrt(GG*mstar/cR**3)                                             # angular velocity
        Tmp = 70 * (cR/(10*au))**(-0.5)
        cs = np.sqrt(kk*Tmp/(2.3333*mH))
        rhomp = rho0 * ((cR/rc)**(-1*gamma))*np.exp(-1*(cR/rc)**(2-gamma))          # dust mass in the midplane
        if cR >= 7.2*au and cR <= 20.2*au:
            rhomp = rhomp/18
        elif cR >= 26.8*au and cR <= 37.8*au:
            rhomp = rhomp/16
        elif cR >= 38.7*au and cR <= 45.3*au:
            rhomp = rhomp/6.9
        elif cR >= 47.75*au and cR <= 52.25*au:
            rhomp = rhomp/3.8
        elif cR >= 58.2*au and cR <= 70.2*au:
            rhomp = rhomp/8
        elif cR >= 69.65*au and cR <= 77.75*au:
            rhomp = rhomp/12
        elif cR >= 86.05*au and cR <= 95.95*au:
            rhomp = rhomp/11
        else:
            pass
        hp = H100min*((amax/0.01)**-0.1)*((cR/(100*au))**beta)                       # scale height in cm
        hm = hp
        integx = math.erf(hm/hp)
        sigma = np.sqrt(2.0*np.pi)*integx*rhomp*hp*100.0                               # convert to gas mass
        QToomre = cs * omega / (np.pi*GG*sigma)
        f.write ("% 13.6e %13.6e %13.6e %13.6e %13.6e\n"%(cR/au, omega,cs,sigma,QToomre))
