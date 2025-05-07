import os, re, numpy as np
from skimage.io import imread
import argparse, yaml
import tifffile
from joblib import Parallel, delayed
from functools import partial

"""
This script is to adapt ome tiff files from Nikon NIS software to be compatible to DART-FISH analysis pipeline (including rename and split channels) and extract the necessary metadata.
"""

def adaptMultiCh(filename, inputdir, outputdir, regex, metadata_file):
    """
    Convert filename to be compatible to DF analysis workflow
    """
    fov_str = 's' + str(int(regex.match(filename).group('fov'))).zfill(n_imgs_width)
    z_str = 'z' + str(int(regex.match(filename).group('z'))).zfill(n_z_width) 
    multich = tifffile.imread(os.path.join(inputdir, filename))
    if multich.shape[0] == min(multich.shape):
        multich = np.moveaxis(multich, 0, -1)
    for ch in range(min(multich.shape)):
        ch_str = "ch{:02}".format(ch)
        new_name = "{0}_{1}_{2}_{3}.tif".format(rnd, fov_str, z_str, ch_str)
        tifffile.imwrite(os.path.join(outputdir, new_name), 
                        multich[..., ch], imagej=True, photometric = 'minisblack', compression="zlib")
def adaptMetadat(filename, inputdir, outputdir, regex, metadata_file):
    """
    Generate metadata file from raw ome tiff file
    """
    fov_str = 's' + str(int(regex.match(filename).group('fov'))).zfill(n_imgs_width)
    z_str = 'z' + str(int(regex.match(filename).group('z'))).zfill(n_z_width) 
    # print(fov_str,z_str)
    with open(metadata_file,'a') as f:
        if int(regex.match(filename).group('z')) == 0:
            for tag in tifffile.TiffFile(os.path.join(inputdir, filename)).pages[0].tags.values():
                if tag.name == 'ImageWidth':
                    ImageWidth = tag.value
                if tag.name == 'ImageLength':
                    ImageLength = tag.value
                if tag.name == '65326': # pixel size tag in ome tiff
                    PixelSize = tag.value
                if tag.name == '65329': # FOV location tag in ome tiff
                    dXPosition = tag.value[0]
                    dYPosition = tag.value[1]
                    dZPosition = tag.value[2]
            n_seq_width = len(regex.match(filename).group('seq'))
            seq_str = 'Seq' + str(int(regex.match(filename).group('seq'))).zfill(n_seq_width)
            new_seq_str = 'Seq' + str(int(regex.match(filename).group('seq'))+1).zfill(n_seq_width)
            next_file = filename.replace('ZStack0000', 'ZStack0001').replace(seq_str, new_seq_str)
            for tag in tifffile.TiffFile(os.path.join(inputdir, next_file)).pages[0].tags.values():
                if tag.name == '65329':
                    dZPosition_new = tag.value[2]
            ZStepSize = abs(dZPosition_new - dZPosition)
            f.write(f'{fov_str},{dXPosition},{dYPosition},{dZPosition},{ImageWidth},{ImageLength},{PixelSize},{ZStepSize}\n')

# Load parameters
parser = argparse.ArgumentParser()
parser.add_argument('param_file')
args = parser.parse_args()
params = yaml.safe_load(open(args.param_file, "r"))

#Input Data Folder
dir_data_raw = params['dir_nikon_data']

dir_output_aligned = params['dir_data_raw'] # "../0_raw_adapted/"# 
if not os.path.isdir(dir_output_aligned):
    os.makedirs(dir_output_aligned)

rnd_list = params['reg_rounds']
rnd_dict = params['rounds_dir']

regex3d = {key:re.compile(params['rounds_regex3d'][key]) for key in params['rounds_regex3d']}

# running the code
for rnd in rnd_list:
    regex = regex3d[rnd]
    # file the corresponding files with regular expression pattern
    files = [file for file in os.listdir(os.path.join(dir_data_raw, rnd_dict[rnd])) 
                                         if regex.match(file)]
    files = sorted(files)
    
    n_imgs_width = len(str(len(set([regex.match(file).group('fov') for file in files]))))
    n_z_width = len(str(len(set([regex.match(file).group('z') for file in files]))))
    n_seq_width = len(str(len(set([regex.match(file).group('seq') for file in files]))))
    rnd_dir = os.path.join(dir_output_aligned, rnd)
    # create a file to save necessary metadata
    metadata_file = os.path.join(dir_output_aligned, 'metadata.csv')
    with open(metadata_file, 'w') as file:
        file.write('Point,dXPosition,dYPosition,dZPosition,ImageWidth,ImageLength,PixelSize,ZStepSize\n')
    if not os.path.exists(rnd_dir):
        os.mkdir(rnd_dir)
    # loop through files to get metadata for all FOVs
    for file in files:
        adaptMetadat(file, inputdir=os.path.join(dir_data_raw, rnd_dict[rnd]), 
                      outputdir=rnd_dir, regex = regex, metadata_file = metadata_file)
    # adapt files to DART-FISH compatible
    adaptPartial = partial(adaptMultiCh, inputdir=os.path.join(dir_data_raw, rnd_dict[rnd]), 
                           outputdir=rnd_dir, regex=regex, metadata_file = metadata_file)
    Parallel(n_jobs=20, prefer='processes')(delayed(adaptPartial)(file) for file in files)
