#!/usr/bin/env python
# ==================================
# Authors:
# Carlos Brandt - chbrandt@lncc.br
# ==================================

"""Package to make FITS image copy/cuts and update their header"""

##@package image
##@file imcp
#
#
# This package contains functions to cutout smaller regions
# from a original image, paste the image into a bigger array,
# as well as any other function for image copy manipulation.
#
# Executable package: YES
#
# To see the package help message, type:
#
# > python imcp.py --help
#
# 'imcp' has one mandatory argument: the (input) image FITS file.
# If none of the options are given, imcp outputs a file named
# 'output.fits' with 1/4th of the original image area, with the 
# same central point/pixel.
#
# > python imcp.py input_file.fits
#
# The options available make possible to define the central point
# where the cut will be centered, as well as the side lengths for
# output image. These values can be passed in 'pixel' or 'degrees'
# units.



import os
import astropy.io.fits as fits
import numpy as np;
from scipy import ndimage as ndi
import sextractor as se
import math as m
try:
    import cv2
except:
    print "cv2 not found, the code will support images in fits format only"


import sys;
import logging;
import astropy.wcs as wcs
from math import sin, cos ,sqrt, fabs, atan, tan 
# =================================================================================================================
def cutout( img, hdr=None, coord_unit='pixel', xo=0, yo=0, size_unit='pixel', x_size=0, y_size=0, mask=None ):
    """
    Do a snapshot from given fits image.

    cutout( ndarray, ...) -> (ndarray, header)

    The function can deal with 'pixel' and 'degrees' cordinate units. As input, besides the image FITS 
    filename (input_file), it receives a central position (xo,yo) and side lengths (x_size,y_size) for
    output image position and dimensioning. Function returns two values: a numpy.array with resultant 
    image pixel intensities and respective header object.

    In addition, an image mask information can be given to use as interest region selection. 'mask' is
    expected to be a numpy.where like structure (i.e, mask=numpy.where()). If given, returned image will
    have all pixels null except the ones listed in 'mask' parameter.

    If 'xo=0' and 'yo=0', given image (input_file) central pixel will be chosen as (xo,yo).
    If 'x_size=0' and 'y_size=0', half length of each side will be used for output dimensions.


    Input:
     - img        numpy.ndarray : Image array (ndim=2,dtype=float)
     - hdr        fits.header : Image (FITS) header
     - coord_unit           str : Position (xo,yo) units ('pixel','degrees')
     - xo                   int : Centroid (horizontal) position for cutting window
     - yo                   int : Centroid (vertical) position for cutting window
     - size_unit            str : Window sizes (x_size,y_size) units ('pixel','degrees')
     - x_size               int : Horizontal cutting window size
     - y_size               int : Vertical cutting window size
     - mask   (ndarray,ndarray) : Tuple with index arrays (Output like numpy.where())

    Output:
     - (ndarray, header) : Resultant image array and (updated) header instance
    
    ---
    """

    image = img;

    logging.debug("Input parameters (type(image),header,coord_unit,xo,yo,size_unit,x_size,y_size,mask): %s",(type(image),hdr,coord_unit,xo,yo,size_unit,x_size,y_size,mask));

    xo=float(xo);
    yo=float(yo);
    x_size=float(x_size);
    y_size=float(y_size);

    imagem = image;
    
    # Initialize some variables..
    #
    x_diff = 0;   x_edge = 0;
    y_diff = 0;   y_edge = 0;
    x_fin = 0;   x_ini = 0;
    y_fin = 0;   y_ini = 0;


    if ( hdr ):
        dimpix = get_pixelscale( hdr );
        logging.info("Pixel_scale: %s",dimpix);

    y_img_size, x_img_size = imagem.shape;
    logging.debug("Input image shape: %s",(x_img_size,y_img_size));

    # Get the side sides (at least, 1!) and transform for the size_unit if necessary..
    #
    x_cut_size = float(x_size);
    y_cut_size = float(y_size);

    if ( size_unit == 'degrees' ):
        x_cut_size = int(float(x_cut_size)/dimpix);
        y_cut_size = int(float(y_cut_size)/dimpix);

    # And if no side size was given, define a default value correspondig to half of original image..
    #
    if not ( x_size ):
        x_cut_size = int(x_img_size/2);
        logging.warning("'x_size' not given. Using half of image x side(%d)", x_cut_size);

    if not ( y_size ):
        y_cut_size = int(y_img_size/2);
        logging.warning("'y_size' not given. Using half of image y side(%d)", y_cut_size);

    logging.debug("Output image shape: %s",(x_cut_size,y_cut_size));

    # Check requested output.vs.input image sizes..
    #
    if ( x_cut_size == x_img_size and y_cut_size == y_img_size ):
        logging.warning("Requested output sizes are the same as input image. Returning image and header unchanged.");
        return (imagem,hdr);

    # Verify central coordinates values..
    #
    if ( coord_unit == 'pixel' and (xo != 0 and yo != 0) ):
        x_halo = int(float(xo));
        y_halo = int(float(yo));

    elif ( coord_unit == 'degrees' and (xo != 0 and yo != 0) ):
        x_halo = int(radec2xy(hdr,xo,yo)[0]);
        y_halo = int(radec2xy(hdr,xo,yo)[1]);

    elif ( xo == 0 and yo == 0 ):
        x_halo = int(x_img_size/2);
        y_halo = int(y_img_size/2);
        logging.warning("No central coordinates were given for snapshot. Using image central pixel as the cut center.");
    else:
        logging.error("Central positioning is out of valid values.");
        return (False,False);

    logging.debug("Central point for output image: %s",(x_halo,y_halo));

    # Define the images (in/out) slices to be copied..
    #
    x_ini = x_halo - int(x_cut_size/2) #-1;
    x_fin = x_ini + x_cut_size;
    y_ini = y_halo - int(y_cut_size/2) #-1;
    y_fin = y_ini + y_cut_size;

    x_ini_old = max( 0, x_ini );   x_fin_old = min( x_img_size, x_fin );
    y_ini_old = max( 0, y_ini );   y_fin_old = min( y_img_size, y_fin );

    x_ini_new = abs( min( 0, x_ini ));   x_fin_new = x_cut_size - (x_fin - x_fin_old);
    y_ini_new = abs( min( 0, y_ini ));   y_fin_new = y_cut_size - (y_fin - y_fin_old);

    # If header, update pixel<->sky information..
    #
    if ( hdr ):
        hdr = update_coordinates(hdr.copy(), x_ini, y_ini);
        hdr['NAXIS1']=x_cut_size#hdr.update('NAXIS1',x_cut_size);
        hdr['NAXIS2']=y_cut_size#hdr.update('NAXIS2',y_cut_size);

    # Initialize new image, and take all index list..
    #
    #print "inside cutout"
    #print y_cut_size
    #print x_cut_size
    #print imagem.dtype
    imagemnova = np.zeros( (int(y_cut_size),int(x_cut_size)), dtype=imagem.dtype );
    ind_z = np.where(imagemnova == 0);

    # Copy requested image slice..
    #
    imagemnova[ int(y_ini_new):int(y_fin_new), int(x_ini_new):int(x_fin_new) ] = imagem[ int(y_ini_old):int(y_fin_old), int(x_ini_old):int(x_fin_old) ];

    # If 'mask', maintain just "central" object on it..
    #
    if ( mask ):
        msk = ( mask[0]-y_ini, mask[1]-x_ini )

        zip_m = zip( msk[0], msk[1] );
        zip_z = zip( ind_z[0], ind_z[1] );

        L = list(set(zip_z) - set(zip_m));

        try:
            ind_0, ind_1 = zip(*L);
            indx = ( np.array(ind_0), np.array(ind_1) );
            imagemnova[ indx ] = 0;
        except:
            pass;

    return (imagemnova, hdr);
def segmentation(filename,method,smooth=0.0,threspar=2,args={},**kwargs):
    if method=='sextractor':
        segimg,objimg,tbdata,hdr,file_name,area=sex_routine4mediatrix(filename,args=args,**kwargs)
    else:
        segimg,objimg,original_img,hdr,file_name,_ids,area=bit_routine4mediatrix(filename, bit_method=method, args=args, sm=smooth,f=threspar, **kwargs)
    return segimg,objimg,original_img,hdr,file_name,_ids,area     


def verify_kwarg(param_name,default_value,kwargs): 
    if param_name in kwargs.keys():
        param=kwargs[param_name]
    else:
        param=default_value
    return param

def sex_routine4mediatrix(image_file, args={},**kwargs):
    preset=verify_kwarg("preset",'HST_Arcs',kwargs)
    origin_dir=verify_kwarg('origin_dir','',kwargs)
    destination_dir=verify_kwarg('destination_dir','',kwargs)
    save_files=verify_kwarg('save_files',False,kwargs)
    
    if  (origin_dir!=destination_dir):
        os.system("cp "+origin_dir+image_file+" "+destination_dir)
    #preset SLchallengeSpace
    Dsex = se.run_segobj(origin_dir+image_file,params=['NUMBER','X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE', 'ELLIPTICITY','MAG_ISO', 'MAG_AUTO','FLAGS'],args=args,preset='HST_Arcs',quiet=True) 
    file_name=image_file.split("/")[-1]
    file_name=file_name.replace(".fits","")

    if  (origin_dir!=destination_dir) and (save_files==True):
        os.system("cp "+file_name+'_seg.fits'+" "+destination_dir)
        #print "cp "+file_name+'_seg.fits'+" "+destination_dir
        os.system("cp "+file_name+'_obj.fits'+" "+destination_dir)
        os.system("cp "+file_name+'_cat.fit'+" "+destination_dir)
    
    Dsex={'SEGMENTATION': os.getcwd()+'/'+file_name+'_seg.fits', 'OBJECTS': os.getcwd()+'/'+file_name+'_obj.fits', 'CATALOG': os.getcwd()+'/'+file_name+'_cat.fit'}
    
    try:    
        segimg =  fits.getdata(Dsex['SEGMENTATION'],ignore_missing_end=True)
        index_seg=np.where(segimg>1)
        Area=len(index_seg[0])
    except:
        print "SEXTRACTOR FILE IS EMPTY: "+ str(image_file)
        segimg=np.array([])
        Area=-1
    try:
        objimg, hdr = fits.getdata(Dsex['OBJECTS'], header=True)
        tbdata = fits.open(Dsex['CATALOG'])[1].data
    except:
        print "SEXTRACTOR FILE IS EMPTY: "+ str(image_file)
        objimg=np.array([])
        hdr=np.array([])
        tbdata = np.array([])
        Area=-1

    os.system("rm "+file_name+'_seg.fits')
    os.system("rm "+file_name+'_obj.fits')
    os.system("rm "+file_name+'_cat.fit')
    return segimg.copy(),objimg.copy(),tbdata,hdr,file_name,Area

def bit_routine4mediatrix(image_file,bit_method='max',sm=0,f=1,args={},**kwargs):

    opt={'increase': 2, 'relative_increase': True,'connected': False,'object_centered':True, 'out_type':'cutout', 'vmin':0 , 'invert':True ,'out_title': 'Mediatrix Decomposition', 'keys_color': "r" ,'alpha': 1 ,'max_level': 1000, 'near_distance': sqrt(2)/2, 'max_level': 1000, 'method':"brightest", 'world_coord': 'False', 'frac': 0.1}
    opt.update(args)
    
    origin_dir=verify_kwarg('origin_dir','',kwargs)
    destination_dir=verify_kwarg('destination_dir','',kwargs)
    if opt['out_type']=='cutout':
        opt['object_centered']=False

    if  (origin_dir!=destination_dir):
        os.system("cp "+origin_dir+image_file+" "+destination_dir)
    test_string=image_file.rstrip(".jpg")
    test_string=test_string.rstrip(".png")
    test_string=test_string.rstrip(".jpeg")
    test_string=test_string.rstrip(".tiff")
    test_string=test_string.rstrip(".gif")
    if test_string== image_file:
        original_img, hdr = fits.getdata(origin_dir+image_file, header=True)
    else:
        original_img=cv2.imread('testimg.jpg')[:,:,0]
        # if it is a rgb image takes only the first channel
        hdr=None
    #original_img, hdr = fits.getdata(origin_dir+image_file, header=True)
    if sm>0:
        sm=fwhm2std(sm)
        original_img=gauss_convolution_fft(original_img, n_fwhm=4, sigma=sm)


    if bit_method=='max':
        segimg_bool=thres(original_img, frac=f)
    elif bit_method=='simple':
        segimg_bool=thres(original_img, thresh=f)
    else:
        x=76
        y=76
        segimg_bool=thres_std(original_img,f=f)
    file_name=image_file.split("/")[-1]
    file_name=file_name.replace(".fits","")
    segimg=np.zeros(segimg_bool.shape)
    segimg[segimg_bool]=1
    index_seg=np.where(segimg==1)
    Area=len(index_seg[0])
    _ids=[1]
    del segimg_bool
    
    if Area>0:
        objimg,hdr=segstamp(segimg=segimg, objID=_ids[0], objimg=original_img, hdr=hdr, \
        increase=opt['increase'], relative_increase=opt['relative_increase'], \
        connected=opt['connected'], obj_centered=opt['object_centered'])
    else:
        _ids=[]
        objimg=np.zeros(segimg.shape)
    fits.writeto(destination_dir+file_name+"_seg.fits",segimg.astype(float),header=None,  overwrite=True)
    fits.writeto(destination_dir+file_name+"_obj.fits",objimg.astype(float),header=None,  overwrite=True)
    return segimg.copy(),objimg.copy(),original_img.copy(),hdr,file_name,_ids,Area


def thres_std(img, thresh=None, size=9,f=1):
    """
    Segment image using a thresholding algorithm
    """
    if thresh is None:
        thresh = Stats.mode(img)
        thresh += f*np.std(img)
    #thresh =threshold_otsu(img)


    # Take the binary image version "splitted" on the 'thresh' value
    img_bin = (img > thresh)#threshold_adaptive(img, 41, offset=10)#(img > thresh)
    
    # Clean fluctuations
    img_bin = Cleaner.spurious(img_bin) 
    
    # Clean small objects
    img_bin = Cleaner.small(img_bin,size=9)
    #for addarc bkg = 40, for paint arcs bkg=9
    
    return img_bin

def thres(img, thresh=None, frac=0.1, size=9):
    """
    Segment image using a thresholding algorithm
    """
    if thresh is None:
        #print "threshold"
        #print np.max(img)
        #print frac
        thresh = np.max(img)*frac
        #thresh += np.std(img)
    
    # Take the binary image version "splitted" on the 'thresh' value
    img_bin = (img > thresh)
    
    # Clean fluctuations
    img_bin = Cleaner.spurious(img_bin) # namespace correction by debom
    
    # Clean small objects
    img_bin = Cleaner.small(img_bin,size=9) # namespace correction by debom
    
    return img_bin


class Distros:
    """
    Namespace to group functions computing/returning distributions
    """
    @staticmethod
    def histogram(img,nbins=1000,normed=True):
        """
        Return image's histogram
        """
        imhist,bins = np.histogram(img.flatten(),bins=nbins,normed=normed)
        return imhist,bins

    @staticmethod
    def cdf(img,nbins=1000):
        """
        Return the Cumulative Distribution Function of given image
        """
        hist,bins = Distros.histogram(img,nbins)
        imcdf = np.cumsum(hist);
        imcdf/=imcdf[-1];
        return imcdf,bins

class Stats:
    """
    Namespace to group functions computing statistical values
    """
    @staticmethod
    def mode(img):
        """
        Return image's mode value
        """
        _hist,bins = Distros.histogram(img)
        ind = np.argmax(_hist)
        _mode = bins[ind]
        return _mode

class Cleaner:

    @staticmethod
    def spurious(img_bin):
        # And use (MO) binary opening (erosion + dilation) for cleaning spurious Trues
        strct = ndi.generate_binary_structure(2,1)
        img_bin = ndi.binary_opening(img_bin,strct)
        return img_bin

    @staticmethod
    def small(img_bin,size=9):
        # Label each group (Regions==True) of pixels
        regions,nlbl = ndi.label(img_bin)
        for i in xrange(1,nlbl+1):
            inds = np.where(regions==i)
            if inds[0].size < size:
                regions[inds] = 0
        return regions.astype(np.bool)




    
# ---
# ==========================================================================================
def segstamp(segimg, objID, objimg=None, hdr=None, increase=0, relative_increase=False, connected=False, obj_centered=True):
    """
    Identify objects on given images by their IDs and return object images

    segstamp( segimg, objID ... )

    By default, if 'objIDs' is not given, postamp will scan segmentation image 
    'seg_img' for the list of object ID numbers. If 'objIDs' is given, those IDs 
    will be used for object poststamps creation.

    'increase' and 'relative_increase' define whether the poststamps will have a 
    size different from object-only dimensions or just the object (pixels) itself, 
    and, also, if this increasement value is a multiplicative factor (True) or an 
    additive one (False).

    Since a list with object IDs is given, a list with arrays, with each IDentified
    object, is returned.

    Input:
     - segimg : numpy.ndarray(ndim=2,dtype=int)
        Segmented image (e.g, SEx's segmentation image)
        
     - objIDs : [int,]
        List with object IDs of interest in 'segimg'.

     - objimg : numpy.ndarray(ndim=2,dtype=float)
        Objects image (e.g, SEx's objects image)

     - hdr : FITS header instance
        FITS header to be updated for each poststamp
        
     - increase : float
        Value for poststamp resizing (> 0)
        
     - relative_increase : bool
        Is 'increase' a additive value (default,False) or multiplicative one(True)?
        

    Output:
     - (ndarray,header)  : Image (poststamp) array and corresponding header

    ---
    """


    _id = objID;

    ind = create_IDmask(segimg, _id);

    y_min = min( ind[0] );
    x_min = min( ind[1] );

    y_idx = ind[0] - y_min;
    x_idx = ind[1] - x_min;

    y_size = max( y_idx ) + 1;
    x_size = max( x_idx ) + 1;
    
    if obj_centered==False:
        pixels = np.where(segimg >=0)
        y_size=segimg.shape[0]#max(pixels[0])
        x_size=segimg.shape[1]#max(pixels[1]) 
       

    if (connected == True ):		
      ind=separate_disconected(ind, high_area=True)
 

    # Central pixel on original image:
    if obj_centered==True:
        yo = y_size/2 + y_min;
        xo = x_size/2 + x_min;
    else:
        yo = y_size/2; 
        xo = x_size/2; 




    if ( increase != 0 and obj_centered==True ):		
        if (relative_increase == True):
            x_size = x_size*increase;
            y_size = y_size*increase;
        else:
            x_size = x_size + 2*increase;
            y_size = y_size + 2*increase;
    
    if type(objimg)!=type('None'):
        #print "imcp param"
        #print type(objimg)
        #print int(xo)
        #print int(yo)
        #print int(x_size)
        #print int(x_size)
        #print ind        
        image_out, hdr = cutout( objimg, hdr, xo=int(xo), yo=int(yo), x_size=int(x_size), y_size=int(y_size), mask=ind );
    else:
        image_out, hdr = cutout( segimg, hdr, xo=int(xo), yo=int(yo), x_size=int(x_size), y_size=int(y_size), mask=ind );
    
    return ( image_out, hdr );

def create_IDmask(segimg, objID):
    """ 
    Generate mask for each object in given image.
    
    ID in 'objIDs' is used to create a mask for each object (ID) in 'segimg'.
    
    Input:
     - segimg : ndarray(ndim=2,dtype=int)
        Image array with int numbers as object identifiers
     - objIDs : [int,]
        List with IDs for objects inside 'segimg'
        
    Output:
     -> index array (output from numpy.where())
        List of tuples with index arrays in it. The output is a list of "numpy.where()" arrays.
        
    """


    # For each object (id) scan the respective indices for image mask and 'cutout' params
    #
    id = float(objID);
    mask = np.where(segimg == int(id));
    
    return mask;

def radec2xy(hdr,ra,dec):

	"""Transforms sky coordinates (RA and Dec) to pixel coordinates (x and y).
	
	Input:
	- hdr: FITS image header
	- ra <float> : Right ascension value in degrees
	- dec <float>: Declination value in degrees
	
	Output:
	- (x,y) <tuple>: pixel coordinates

	"""

	wcs = wcs.WCS(hdr)
	
	skycrd = np.array([[ra,dec]])

	pixcrd = wcs.wcs_sky2pix(skycrd,1)

	x = pixcrd[0][0]

	y = pixcrd[0][1]

	return (x,y)
	
	


def xy2radec(hdr,x,y):

	""" Transforms pixel coordinates (x and y) to sky coordinates (RA and Dec).
	
	Input:
	- hdr: FITS image header
	- x <float>: x value
	- y <float>: y value
	
	Output:
	- (ra,dec) <tuple>: Right ascension and declination values in degrees 

	"""

	wcs = wcs.WCS(hdr)

	pixcrd = np.array([[x,y]], np.float_)

	skycrd = wcs.wcs_pix2sky(pixcrd,1)

	ra = skycrd[0][0]

	dec = skycrd[0][1]

	return (ra,dec)	

def get_pixelscale(hdr):
    """ Read out header dimpix (pixel_scale) from CD_* (or CDELT*) keys

    Input:
     - hdr  <.Header> : Image header instance

    Output:
     - pixel_scale  float : Image sky-to-pixel scale ("/px)

    ---
    """

    try:
        pixel_scale = hdr['PIXSCALE'];
        return pixel_scale;
    except:
        pass;
        
    try:
        CD1_1 = float(hdr['CD1_1']);
        CD1_2 = float(hdr['CD1_2']);
        CD2_1 = float(hdr['CD2_1']);
        CD2_2 = float(hdr['CD2_2']);
    except:
        CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1;

    if ( CD1_1 == 1 and CD1_2 == 1 and CD2_1 == 1 and CD2_2 == 1 ):
        try:
            CD1_1 = float(hdr['CDELT1']);
            CD1_2 = 0.;
            CD2_1 = 0.;
            CD2_2 = float(hdr['CDELT2']);
        except:
            CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1;
            print "Unable to find the pixel size value";

    try:
        CUNIT1 = str(hdr['CUNIT1']);
        CUNIT2 = str(hdr['CUNIT2']);
    except:
        print >> sys.stderr, "Warning: Unable to find 'CUNIT[1|2]' parameters in header instance for image coordinates unit.";
        print >> sys.stderr, "         Degrees ('deg') is being used as default unit value.";
        CUNIT1 = 'deg';
        CUNIT2 = 'deg';

    conv_factor = 1.;
    CUNIT = CUNIT1
    # Convertion from degrees to arcsec units..
    #
    if ( CUNIT1 != CUNIT2 ):
        print >> sys.stderr, "Error: I do not have support for different per-axis pixel resolutions.";
	print >> sys.stderr, "Error: CUNIT1: %s , CUNIT2: %s" % (CUNIT1,CUNIT2)
	return None

    # Scale is given, using CD1_1, CD1_2, CD2_1, CD2_2 header keys:
    #
    pixel_scale = m.sqrt((CD1_1**2 + CD1_2**2 + CD2_1**2 + CD2_2**2)/2.) #* conv_factor;
    print "Unit: %s/px" % CUNIT

    return pixel_scale;


def update_coordinates(hdr, x_ini, y_ini):
    """Update header information regarding world coordinate system"""


    LTV1 = None;
    LTV2 = None;

    NAXIS1 = int(hdr['NAXIS1']);
    NAXIS2 = int(hdr['NAXIS2']);

    try:
        CRPIX1 = float(hdr['CRPIX1']);
        CRPIX2 = float(hdr['CRPIX2']);
    except:
        CRPIX1 = CRPIX2 = 1.0;

    try:
        CD1_1 = float(hdr['CD1_1']);
        CD1_2 = float(hdr['CD1_2']);
        CD2_1 = float(hdr['CD2_1']);
        CD2_2 = float(hdr['CD2_2']);
    except:
        CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1.0;

    LTV1 = -1*x_ini;
    LTV2 = -1*y_ini;
    CRPIX1 = CRPIX1 + LTV1;
    CRPIX2 = CRPIX2 + LTV2;


    # Add some keys to header..
    #
    WCSDIM = 2
    CDELT1  =  CD1_1;
    CDELT2  =  CD2_2;
    LTM1_1 = 1.0
    LTM2_2 = 1.0
    WAT0_001= 'system=image'
    WAT1_001= 'wtype=tan axtype=ra'
    WAT2_001= 'wtype=tan axtype=dec'


    # Header update..
    #
    if (LTV1 != None) :
        hdr['LTV1']=LTV1##.update('LTV1',LTV1)
        hdr['CRPIX1']=CRPIX1
    if (LTV2 != None) :
        hdr['LTV2']=LTV2
        hdr['CRPIX2']=CRPIX2
    hdr['WCSDIM']=WCSDIM
    hdr['CDELT1']=CDELT1
    hdr['CDELT2']=CDELT2
    hdr['LTM1_1']=LTM1_1
    hdr['LTM2_2']=LTM2_2
    hdr['WAT0_001']=WAT0_001
    hdr['WAT1_001']=WAT1_001
    hdr['WAT2_001']=WAT2_001

    return (hdr);





def separate_disconected(ind, high_area=False):

    """
    From a list of points that belongs to objects it separates in groups of connected objects. If high_area=False return a list which first index represents each object , the second and third its coordinates. if high_area=True  returns the only a list with the coordinates of the object that has greater area.

   
    Input:
     - ind : array like
        the list of points where ind[0] is a list with the first coordinates (e.g. x) and ind[1] the second coordinates (e.g. y)
     - high_area : bool
        enables a criteria to return only the coordinates of the object that has greater area
   
    Output:
     - (nparray)  : if high_area=False  is list which first index represent the  object.  The second  and third index represents  the first and second coordinates. if high_area=True  is list which first and second index represents lists the first and second coordinates

    ---
    """

    p=[]
    for i in range(0,len(ind[0])):
        p_aux=[ind[0][i],ind[1][i]]
        p.append(p_aux)
   
  
    Objects=[]
    objects=[]
    while len(p)>0:
        p_test=[p[0]]
        while len(p_test)>0:
            p_center=p_test[0]
            p_neighbors=[[p_center[0],p_center[1]],[p_center[0]+1,p_center[1]],[p_center[0]+1,p_center[1]+1],[p_center[0]+1,p_center[1]-1],[p_center[0]-1,p_center[1]],[p_center[0]-1,p_center[1]+1],[p_center[0]-1,p_center[1]-1],[p_center[0],p_center[1]-1],[p_center[0],p_center[1]+1]]
            for i in range(0,len(p_neighbors)):
                if (p_neighbors[i] in p) and not(p_neighbors[i] in objects):
                    objects.append(p_neighbors[i])
                    p_test.append(p_neighbors[i])
                    p.remove(p_neighbors[i])  
            p_test.remove(p_test[0])
        Objects.append(objects)
        objects=[]

    if high_area==False:
        return Objects
    else:
        # criteria to select an interest object
   
        Area_max=len(Objects[0])
        id_max=0
        for i in range(0,len(Objects)):
            if len(Objects[i])>Area_max:
                Area_max=len(Objects[i])
                id_max=i
        ind_new=[[-100],[-1000]]
        for i in range(0,len(Objects[id_max])):
            ind_new[0].append(Objects[id_max][i][0])
            ind_new[1].append(Objects[id_max][i][1])
        ind_new[0].remove(-100)
        ind_new[1].remove(-1000)
        ind_n=[np.array(ind_new[0]),np.array(ind_new[1])]
        return ind_n          

