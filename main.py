# Reuses some info from Stephen Bailey shared on [desi-data 3401] "running fiber assignment on a real target catalog"
import os
import subprocess
from astropy.table import Table
import numpy as np
from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, obsmask, obsconditions
import fitsio
import glob
from desisim.quickcat import quickcat

# target selection
targetfile = "data/dr6_targets.fits"
os.environ['DECALS_PATH'] = '/global/project/projectdirs/cosmo/data/legacysurvey/dr6/'

if not os.path.exists(targetfile):
    cmd = "select_targets {source} {destination}"
    cmd = cmd.format(source=os.environ['DECALS_PATH'], destination=targetfile)
    print(cmd)
    subprocess.call(cmd.split())
    print('Done with target selection')

columns = [
    'TARGETID', 'RA', 'DEC', 'SUBPRIORITY', 'BRICKNAME',
    'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET',
]

#truth file
truthfile = "data/truth.fits"
if not os.path.exists(truthfile):
    import desitarget.mock.mockmaker as mb
    from desitarget.targetmask import desi_mask, bgs_mask, mws_mask

    #targetsfilename = "small_chunk_targets-dr5.0-0.16.2.fits"
    targets = fitsio.read(targetfile, 'TARGETS', columns=columns)
    colnames = list(targets.dtype.names)
    print(colnames)
    nobj = len(targets)
    truth = mb.empty_truth_table(nobj=nobj)
    print(truth.keys())

    for k in colnames:
        if k in truth.keys():
            print(k)
            truth[k][:] = targets[k][:]

    nothing = '          '
    truth['TEMPLATESUBTYPE'] = np.repeat(nothing, nobj)

    masks = ['BGS_ANY', 'ELG', 'LRG', 'QSO', 'STD_FSTAR', 'STD_BRIGHT']
    dict_truespectype = {'BGS_ANY':'GALAXY', 'ELG':'GALAXY', 'LRG':'GALAXY', 'QSO':'QSO', 
                    'STD_FSTAR':'STAR', 'STD_BRIGHT':'STAR'}
    dict_truetemplatetype = {'BGS_ANY':'BGS', 'ELG':'ELG', 'LRG':'LRG', 'QSO':'QSO', 
                        'STD_FSTAR':'STAR', 'STD_BRIGHT':'STAR'}

    for m in masks:
        istype = (targets['DESI_TARGET'] & desi_mask.mask(m))!=0
        print(m, np.count_nonzero(istype))
        truth['TRUESPECTYPE'][istype] = np.repeat(dict_truespectype[m], np.count_nonzero(istype))
        truth['TEMPLATETYPE'][istype] = np.repeat(dict_truetemplatetype[m], np.count_nonzero(istype))
        truth['MOCKID'][istype] = targets['TARGETID'][istype]

    del targets
    print('writing truth')
    truth.write(truthfile, overwrite=True)
    print('done truth')


#generate sky data
from desitarget.mock.sky import random_sky
from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, obsmask, obsconditions
import desimodel.footprint

skyfile = 'data/dense_sky.fits'
if not os.path.exists(skyfile):

    skyra, skydec = random_sky(nside=4096)
    n_p = len(skyra)

    data = Table()
    data['RA'] = skyra
    data['DEC'] = skydec
    data['TARGETID'] = np.int_(1E10+ np.arange(n_p))
    data['DESI_TARGET'] = np.int_(np.ones(n_p) * desi_mask['SKY'])
    data['MWS_TARGET'] = np.int_(np.zeros(n_p))
    data['BGS_TARGET'] = np.int_(np.zeros(n_p))
    data['OBSCONDITIONS'] = np.int_(np.ones(n_p)*(obsconditions['DARK']|obsconditions['BRIGHT']|obsconditions['GRAY']))
    data['BRICKNAME'] = np.chararray(n_p, itemsize=8)
    data['BRICKNAME'][:] = "b0000000"
    data['SUBPRIORITY'] = np.random.random(n_p)
    data.write(skyfile, overwrite=True)
    print('Done creating sky data')

#reloading target data

mtlfile = 'data/mtl.fits'
darkfile = 'data/std-dark.fits'
brightfile = 'data/std-bright.fits'

if (not os.path.exists(mtlfile)) or (not os.path.exists(brightfile)) or (not os.path.exists(darkfile)):
    targetdata = fitsio.read(targetfile, 'TARGETS', columns=columns)
    print('Done reading target data')

#compute MTL
if not os.path.exists(mtlfile):
    import desitarget.mtl
    mtl = desitarget.mtl.make_mtl(targetdata)
    mtl.meta['EXTNAME'] = 'MTL'
    mtl.write(mtlfile)

    #print some stats
    print('MWS_TARGETS: {}'.format(np.count_nonzero(targetdata['MWS_TARGET']!=0)))
    print('BGS_TARGETS: {}'.format(np.count_nonzero(targetdata['BGS_TARGET']!=0)))
    print('DESI_TARGETS: {}'.format(np.count_nonzero(targetdata['DESI_TARGET']!=0)))

#standards
if not os.path.exists(darkfile):
    darkstd = (targetdata['DESI_TARGET'] & desi_mask.mask('STD_FSTAR|STD_WD')) != 0
    darkdata = targetdata[darkstd]
    obscond = np.zeros(len(darkdata), dtype=np.int64)
    darkdata = np.lib.recfunctions.append_fields(
    darkdata, 'OBSCONDITIONS', obscond)
    fitsio.write(darkfile, darkdata, extname='STD')
    print('{} dark standards'.format(np.count_nonzero(darkstd)))
    print('Finished with dark standards')
    
if not os.path.exists(brightfile):
    brightstd = (targetdata['DESI_TARGET'] & desi_mask.mask('STD_BRIGHT')) != 0
    brightdata = targetdata[brightstd]
    obscond = np.zeros(len(brightdata), dtype=np.int64)
    brightdata = np.lib.recfunctions.append_fields(
    brightdata, 'OBSCONDITIONS', obscond)
    fitsio.write(brightfile, brightdata, extname='STD')
    print('Done with bright standards')
    print('{} bright standards'.format(np.count_nonzero(brightstd)))
    
    
#creating tile lists 
import desimodel.io
tiles = desimodel.io.load_tiles()
dx = open('data/dark-tiles.txt', 'w')
bx = open('data/bright-tiles.txt', 'w')
for tileid, program  in zip(tiles['TILEID'], tiles['PROGRAM']):
    if program == 'BRIGHT':
        bx.write(str(tileid)+'\n')
    else:
        dx.write(str(tileid)+'\n')

dx.close()
bx.close()

#run fiberassign
output_bright = 'output/bright/'

if not os.path.exists(output_bright):
    os.makedirs(output_bright)
    cmd = "fiberassign --mtl data/mtl.fits "
    cmd += " --sky data/dense_sky.fits --stdstar data/std-bright.fits "
    cmd += " --footprint /global/common/cori/contrib/desi/desiconda/20170613-1.1.4-spectro/code/desimodel/0.9.0/data/footprint/desi-tiles.fits " 
    cmd += " --surveytiles data/bright-tiles.txt --positioners /global/common/cori/contrib/desi/desiconda/20170613-1.1.4-spectro/code/desimodel/0.9.0/data/focalplane/fiberpos.txt "
    cmd += "--fibstatusfile data/fiberstatus.ecsv --outdir output/bright/"
    print('starting fiberassign for bright tiles')
    subprocess.call(cmd.split())
    print('finished fiberassign for bright tiles')

    
output_bright = 'output/dark/'

if not os.path.exists(output_bright):
    os.makedirs(output_bright)
    cmd = "fiberassign --mtl data/mtl.fits "
    cmd += " --sky data/dense_sky.fits --stdstar data/std-dark.fits "
    cmd += " --footprint /global/common/cori/contrib/desi/desiconda/20170613-1.1.4-spectro/code/desimodel/0.9.0/data/footprint/desi-tiles.fits " 
    cmd += " --surveytiles data/dark-tiles.txt --positioners /global/common/cori/contrib/desi/desiconda/20170613-1.1.4-spectro/code/desimodel/0.9.0/data/focalplane/fiberpos.txt "
    cmd += "--fibstatusfile data/fiberstatus.ecsv --outdir output/dark/"
    print('starting fiberassign for dark tiles')
    subprocess.call(cmd.split())
    print('finished fiberassign for dark tiles')
    
    
zcat_dark = 'data/zcat_dark.fits'
if not os.path.exists(zcat_dark):
    tile_files = glob.glob('output/dark/tile_*.fits')
    print('{} files to gather'.format(len(tile_files)))

    print('reading mtl')
    mtl = Table.read('data/mtl.fits')
    print('reading truth')
    truth = Table.read('data/truth.fits')

    print('making zcat')
    zcat = quickcat(tile_files, mtl, truth, perfect=True)
    print('writing zcat')
    zcat.write(zcat_dark, overwrite=True)
    print('finished zcat')
    
zcat_bright = 'data/zcat_bright.fits'
if not os.path.exists(zcat_dark):
    tile_files = glob.glob('output/bright/tile_*.fits')
    print('{} files to gather'.format(len(tile_files)))

    print('reading mtl')
    mtl = Table.read('data/mtl.fits')
    print('reading truth')
    truth = Table.read('data/truth.fits')

    print('making zcat')
    zcat = quickcat(tile_files, mtl, truth, perfect=True)
    print('writing zcat')
    zcat.write(zcat_bright, overwrite=True)
    print('finished zcat')





