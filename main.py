import os, re
import shutil, sys
from datetime import datetime
import yaml
import subprocess as sp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('param_file')
args = parser.parse_args()
params = yaml.safe_load(open(args.param_file, "r"))
param_file = args.param_file

""" Running maximum intensity projection"""
mip_shell = [sys.executable, "file_adapter.py", param_file]
print(mip_shell)
commandOut = sp.run(mip_shell, stderr = sp.PIPE, text = True)

print(commandOut.stderr)

""" Running maximum intensity projection"""
mip_shell = [sys.executable, "rawMaxProjection.py", param_file]
print(mip_shell)
commandOut = sp.run(mip_shell, stderr = sp.PIPE, text = True)

print(commandOut.stderr)

if commandOut.returncode:
    print("Error in MIPing.")
    sys.exit(1)

""" Running image registration"""
align_shell = [sys.executable, "AlignerDriver_2D.py", param_file]
print(align_shell)
commandOut = sp.run(align_shell, stderr = sp.PIPE, text = True)

print(commandOut.stderr)

if commandOut.returncode:
    print("Error in image registration.")
    sys.exit(1)


# background substraction step is usually skipped in my DART-FISH experiment.

# """ Running background subtraction (optional)"""
# if params['background_subtraction']:
#     subt_shell = [sys.executable, "backgroundSubtraction.py", param_file]
#     print(subt_shell)
#     commandOut = sp.run(subt_shell, stderr = sp.PIPE, text = True)
 
#     print(commandOut.stderr)

#     if commandOut.returncode:
#         print("Error in background subtraction:")
#         sys.exit(1)


    
""" Running stitching"""
stitch_shell = [sys.executable, "StitchDriver_2D.py", param_file]
print(stitch_shell)
commandOut = sp.run(stitch_shell, stderr = sp.PIPE, text = True)

print(commandOut.stderr)
if commandOut.returncode:
    print("Error in stitching.")
    sys.exit(1)

""" Run the pre-decoding code"""
pdc_shell = [sys.executable, "findDecodingParams.py", param_file]
print(pdc_shell)
commandOut = sp.run(pdc_shell, stderr = sp.PIPE, text = True)

print(commandOut.stderr)
if commandOut.returncode:
    print("Error in findDecodingParams.py.")
    sys.exit(1)

""" Run decoding"""
dc_shell = ["docker job --type r6i.32xlarge --cpu 15 --memory 500",sys.executable, "decoding_driver.py", param_file]
print(dc_shell)
commandOut = sp.run(dc_shell, stderr = sp.PIPE, text = True)

print(commandOut.stderr)
if commandOut.returncode:
    print("Error in sparse decoding.")
    sys.exit(1)

""" Combinding FOVs"""
comb_shell = [sys.executable, "CombineFOVs.py", param_file]
print(comb_shell)
commandOut = sp.run(comb_shell, stderr = sp.PIPE, text = True)


print(commandOut.stderr)

if commandOut.returncode:
    print("Error in combining FOVs.")
    sys.exit(1)

""" Running segmentation and assignment"""
seg_shell = [sys.executable, "segmentation_driver.py", param_file]
print(seg_shell)
commandOut = sp.run(seg_shell, stderr = sp.PIPE, text = True)

print(commandOut.stderr)
if commandOut.returncode:
    print("Error in segmentation/assignment")
    sys.exit(1)


""" Running QC plots """
qc_shell = [sys.executable, "QC_plots.py", param_file]
print(qc_shell)
commandOut = sp.run(qc_shell, stderr = sp.PIPE, text = True)

print(commandOut.stderr)
if commandOut.returncode:
    print("Error in QC plots.")
    sys.exit(1)    