from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter


import argparse


parser = argparse.ArgumentParser(description='Configure the parameters of the execution.')
parser.add_argument("--RA", help="RA of the target",default=181.4184,type=float)
parser.add_argument("--Dec", help="Dec of the target",default=49.1748,type=float)
parser.add_argument("--zl", help="Redshift of the deflector",default=0.405, type=float)
parser.add_argument("--zs", help="Redshift of the source",default=0.405, type=float)
parser.add_argument("--lmin", help="Minimum wavelength (Angstroms)",default=3800, type=float)
parser.add_argument("--lmax", help="Maximum wavelength (Angstroms)",default=9500, type=float)
parser.add_argument("--plotlinessource", help="True to plot lines on the residuals",default=False, type=bool)
parser.add_argument("--plotlineslens", help="True to plot lines on the data",default=False, type=bool)
args = parser.parse_args()
linenames = ['Lyb', 'OVI', 'Lya', 'NV', 'OI_1', 'CII_1', 'SiIV', 'SiIV\ \n+OIV', 'CIV', 'HeII', 'OIII_1', 'AlIII', 'CIII', 'CII_2', 'FeII_1', 'NeIV', 'FeII_2', 'MgII', 'NeV', 'NeVI', 'OII_1', 'OII_2', 'NeIII', 'HeI', 'SII_1', 'Hd', 'Hg', 'OIII_2', 'Hb', 'OIII_3', 'OIII_4', 'OIII_5', 'HeI_1', 'OI_2', 'OI_3', 'NI', 'NII_1', 'Ha', 'NII_2', 'SII_2', 'SII_3', 'Ca H', 'K', 'G', 'Mg', 'Na', 'NIII']
lines = [1025.05577, 1033.160518, 1215.075915, 1240.220377, 1304.948362, 1334.730576, 1397.032745, 1399.222767, 1548.897336, 1639.808838, 1665.255924, 1856.778311, 1908.103559, 2325.285565, 2382.037919, 2438.760078, 2599.395781, 2798.292088, 3345.828232, 3425.867762, 3726.032278, 3728.814555, 3868.760, 3887.898111, 4071.150115, 4101.732081, 4340.45917, 4363.209158, 4861.32092, 4931.225284, 4958.909898, 5006.842106, 5875.624, 6300.300802, 6363.773682, 6527.223571, 6548.047948, 6562.793967, 6583.448389, 6716.432467, 6730.808582, 3933.663149, 3968.465041, 4304.398695, 5175.257103, 5893.964261, 1748]

ras = [159.1742083]
ras = [181.4184]
decs = [46.5416139]
decs = [49.1748]

ras = [112.5627]
decs = [41.7480]

ras = [args.RA]
decs = [args.Dec]
z = args.zs
zl = args.zl
wavlims = [args.lmin, args.lmax]

# z = 0.205
plot_lines_source = args.plotlinessource
plot_lines_lens = args.plotlineslens


for i in range(len(ras)):
    print(i)
    ra, dec = ras[i], decs[i]
    pos = coords.SkyCoord(ra, dec,  unit='deg')

    xid = SDSS.query_region(pos, spectro=True,radius=60*u.arcsec)
    print("xid",xid)
    print(r'len of xid: {len(xid}')
    if len(xid) < 1:
        raise Exception(f'No spectra available for that coordinates {args.RA},{args.Dec}')
    xid = xid[0]
    sp = SDSS.get_spectra(plate=xid['plate'], mjd=xid['mjd'], fiberID=xid['fiberID'])
    flux, model, logwav, ivar = sp[0][1].data['flux'], sp[0][1].data['model'], sp[0][1].data['loglam'], sp[0][1].data['ivar']

    #get ylims
    kappa = np.where((logwav>np.log10(wavlims[0]))&(logwav<np.log10(wavlims[1])))[0]
    ymin, ymax = np.min(flux[kappa]), np.max(flux[kappa])

    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes([0.05, 0.7, 0.9, 0.2])

    #plot the data and best fit model from SDSS
    ax.plot(10.**logwav, flux, color='red', alpha=0.5)
    ax.plot(10.**logwav, model, color='red')
    ax.set_title('RA, DEC='+str(ra)+', '+str(dec))
    ax.set_xlim(args.lmin, args.lmax)
    ax.set_ylim(ymin, ymax)
    textheight = [10,5] * len(linenames)

    if plot_lines_lens:
        for k,line in enumerate(lines):
            #if linenames[k] in ['OII_1', 'OII_2',  'OIII_3', 'OIII_4', 'OIII_5',  'OIII_2', 'Hg', 'Hb', 'Ha']:
            if True:
                ax.vlines(x=line*(1+zl), ymin=0, ymax=1, color='black', ls='--', transform=ax.get_xaxis_transform(), alpha=0.5)
                plt.text(line*(1+zl),textheight[k],linenames[k])
 


    #plot the residuals (data-model)
    ax = plt.axes([0.05, 0.4, 0.9, 0.2])
    ax.plot(10.**logwav, flux-model, color='black')
    ax.set_title('data-model')
    ax.set_xlim(args.lmin, args.lmax)
    ymin, ymax = np.min((flux-model)[kappa]), np.max((flux-model)[kappa])
    ax.set_ylim(ymin, ymax)




    #plot the residual signal-to-noise and a smoothed version of it
    ax = plt.axes([0.05, 0.1, 0.9, 0.2])
    textheight = [6,3] * len(linenames)
    ax.plot(10.**logwav, (flux-model)*(ivar**0.5), color='green', alpha=0.5)
    ax.plot(10.**logwav, gaussian_filter((flux-model)*(ivar**0.5), 2), color='green')
    ax.set_title('(data-model)/noise')
    ax.set_xlim(args.lmin, args.lmax)
    if plot_lines_source:
        for k,line in enumerate(lines):
            if linenames[k] in ['OII_1', 'OII_2',  'OIII_3', 'OIII_4', 'OIII_5',  'OIII_2', 'Hg', 'Hb', 'Ha']:
                ax.vlines(x=line*(1+z), ymin=0, ymax=1, color='black', ls='--', transform=ax.get_xaxis_transform(), alpha=0.5)
                plt.text(line*(1+z),textheight[k],linenames[k])
    plt.show(block=True)
