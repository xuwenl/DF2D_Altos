# Kian: added 201011 after having trouble importing SimpleITK normally
import sys
import SimpleITK as sitk
import tifffile
import numpy as np
import os, re
from os.path import join as pathjoin
from datetime import datetime


class TwoDimensionalAligner():
    def __init__(self, movingImgFiles, refImgFiles, transform = "affine", NumberOfResolutions = 5, MaximumNumberOfIterations = 1000, 
                    NumberOfSpatialSamples = 4000, transParamFile = None, voxelSize = (0.0001, 0.0001, 0.001), AutomaticScale = True, Scales=['2']):
        self.movingImgFiles = movingImgFiles
        self.refImgFiles = refImgFiles
        self.transform = transform
        self.NumberOfResolutions = NumberOfResolutions
        self.MaximumNumberOfIterations = MaximumNumberOfIterations
        self.NumberOfSpatialSamples = NumberOfSpatialSamples
        self.transParamFile = transParamFile
        self.voxelSize = voxelSize
        self.AutomaticScale = AutomaticScale
        self.Scales = Scales
        
    def findTransformParameters(self):   
        """ running elastix on moving and ref images to find the transform parameter map between them """
        self.elastixImageFilter = sitk.ElastixImageFilter() # The basic object to do the transformation
        
        """ Setting the transformation parameters"""        
        parameterMap = self.elastixImageFilter.GetDefaultParameterMap(self.transform) # getting the dafault parameter map for our transformation of interest
        parameterMap['NumberOfHistogramBins'] = ['64'] # a parameter for the image comparison metric, AdvancedMattesMutualInformation, that we are using.
        parameterMap['MaximumNumberOfIterations'] = [str(self.MaximumNumberOfIterations)] # number of iterations per aligning each resolution
        parameterMap['NumberOfResolutions'] = [str(self.NumberOfResolutions)] # number of resolution-decreasing alignments. This is the most critical parameter
        parameterMap['NumberOfSpatialSamples'] = [str(self.NumberOfSpatialSamples)] # number of random samples drawn for image comparison during optimization
        parameterMap['WriteIterationInfo'] = ['true'] # This command writes the report in the current working directory, so we have to move the files later    
        parameterMap['AutomaticTransformInitialization'] = ['true']
        parameterMap['AutomaticTransformInitializationMethod'] = ['CenterOfGravity']
        if self.AutomaticScale: # if the scale of parameters is to be found automatically
            parameterMap['AutomaticScalesEstimation'] = ["true"]
        else:
            parameterMap['AutomaticScalesEstimation'] = ["false"]
            parameterMap['Scales'] = self.Scales  # the rotations we expect (in radians) are almost in the same order of magnitude as the translations we expect (in mm)

        
        self.elastixImageFilter.SetParameterMap(parameterMap) # setting the parameter map to our transformation object
        self.elastixImageFilter.SetMovingImage(self.readMovingImage()) # Setting the origin image, the one we want to transform
        self.elastixImageFilter.SetFixedImage(self.readRefImage()) # Setting the destination/final image
        
        self.elastixImageFilter.LogToFileOn()
        self.elastixImageFilter.SetOutputDirectory("./")

        self.elastixImageFilter.Execute()   # running the transformation
        self.transformParameterMap = self.elastixImageFilter.GetTransformParameterMap() # saving the optimized transformation parameters
        
        if not self.transParamFile is None:
            self.writeParameterFile(reportName = self.transParamFile) # write the parameter map to folder MetaData

    def applyTransformation(self, movingImgFiles_list, outImgFiles_list):
        """ Applies the computed transform to all the moving images within movingImgFiles_list. 
            movingImgFiles_list is a list of lists, each inner list contains one stack.
            outImgFiles_list is a list of lists, each inner list specifies the path to the transformed images.
        """
        self.transformixImageFilter = sitk.TransformixImageFilter() 
        self.transformParameterMap = self.getTransformParameterMap()
        self.transformParameterMap[0]['FinalBSplineInterpolationOrder'] = ['1']
        self.transformixImageFilter.SetTransformParameterMap(self.transformParameterMap) 
        
        for movingImgFiles, outImgFiles in zip(movingImgFiles_list, outImgFiles_list):
            images2D_input = self.readImage(movingImgFiles)
            # images2D_input = sitk.Extract(images3D_input, (images3D_input.GetWidth(), images3D_input.GetHeight(), 0), (0,0,0))
            self.transformixImageFilter.SetMovingImage(images2D_input)
            self.transformixImageFilter.Execute()

            outputImage = self.transformixImageFilter.GetResultImage()
            self.writeImage(sitk.Cast(outputImage, sitk.sitkUInt16), outImgFiles)
            # sitk.WriteImage(sitk.Cast(outputImage, sitk.sitkUInt8),
            #                 outImgFiles)
        
    
    def readMovingImage(self):
        im3d = sitk.ReadImage(self.movingImgFiles)
        im3d.SetSpacing(self.voxelSize)
        im2d = sitk.Extract(im3d, (im3d.GetWidth(), im3d.GetHeight(), 0), (0,0,0))
        return(im2d)
        # im = sitk.DiscreteGaussian(sitk.Cast(im, sitk.sitkFloat32), 3, 32, 0.01, False) # sigma, kernel width, max error, useImageSpacing
        # return(sitk.Cast(im, sitk.sitkUInt8))
        
    def readRefImage(self):
        im3d = sitk.ReadImage(self.refImgFiles)
        im3d.SetSpacing(self.voxelSize)
        im2d = sitk.Extract(im3d, (im3d.GetWidth(), im3d.GetHeight(), 0), (0,0,0))
        return(im2d)
        # im = sitk.DiscreteGaussian(sitk.Cast(im, sitk.sitkFloat32), 3, 32, 0.01, False) # sigma, kernel width, max error, useImageSpacing
        # return(sitk.Cast(im, sitk.sitkUInt8))
        

    def readImage(self, path):
        im3d = sitk.ReadImage(path)
        im3d.SetSpacing(self.voxelSize)
        im2d = sitk.Extract(im3d, (im3d.GetWidth(), im3d.GetHeight(), 0), (0,0,0))
        return(im2d)        
    
    def writeParameterFile(self, reportName):
        self.elastixImageFilter.WriteParameterFile(parameterMap = self.transformParameterMap[0], filename = reportName)
     
    def getTransformParameterMap(self):
        return self.transformParameterMap

    def writeImage(self, img_sitk, outfiles):
        img_np = sitk.GetArrayFromImage(img_sitk).astype(np.uint16)
        tifffile.imwrite(outfiles[0], img_np, imagej=True, photometric = 'minisblack', compression="zlib") # "LZW")
        # for im, out in zip(img_np, outfiles):
        #     # print('X', out)
        #     tifffile.imwrite(out, im, imagej=True, photometric = 'minisblack', compression="zlib") # "LZW")

