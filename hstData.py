#!/usr/bin/env python

def get(dic, key=None, value=None):

    if dic.has_key('flux'):

        # return the data array with flux
        if key == 'lineID':
            alist = [{'lineID':lineID, 'fileName':fileName, 'image':image, 'xScale':xScale, 'yScale':yScale, 'pixScale':pixScale, 'JD':JD, 'MJD':MJD, 'decYr':decYr, 'Phi':Phi, 'flux':flux} for lineID, fileName, image, xScale, yScale, pixScale, JD, MJD, decYr, Phi, flux in zip(dic['lineID'], dic['fileName'], dic['image'], dic['xScale'], dic['yScale'], dic['pixScale'], dic['JD'], dic['MJD'], dic['decYr'], dic['Phi'], dic['flux']) if lineID == value]
        if key == 'JD':
            alist = [{'lineID':lineID, 'fileName':fileName, 'image':image, 'xScale':xScale, 'yScale':yScale, 'pixScale':pixScale, 'JD':JD, 'MJD':MJD, 'decYr':decYr, 'Phi':Phi, 'flux':flux} for lineID, fileName, image, xScale, yScale, pixScale, JD, MJD, decYr, Phi, flux in zip(dic['lineID'], dic['fileName'], dic['image'], dic['xScale'], dic['yScale'], dic['pixScale'], dic['JD'], dic['MJD'], dic['decYr'], dic['Phi'], dic['flux']) if JD == value]

    elif (dic.has_key('dist')):

        # return the data array with dist, coords & pa
        if key == 'lineID':
            alist = [{'lineID':lineID, 'fileName':fileName, 'image':image, 'xScale':xScale, 'yScale':yScale, 'pixScale':pixScale, 'JD':JD, 'MJD':MJD, 'decYr':decYr, 'Phi':Phi, 'dist':dist, 'coords':coords, 'pa':pa} for lineID, fileName, image, xScale, yScale, pixScale, JD, MJD, decYr, Phi, dist, coords, pa in zip(dic['lineID'], dic['fileName'], dic['image'], dic['xScale'], dic['yScale'], dic['pixScale'], dic['JD'], dic['MJD'], dic['decYr'], dic['Phi'], dic['dist'], dic['coords'], dic['pa']) if lineID == value]
        if key == 'JD':
            alist = [{'lineID':lineID, 'fileName':fileName, 'image':image, 'xScale':xScale, 'yScale':yScale, 'pixScale':pixScale, 'JD':JD, 'MJD':MJD, 'decYr':decYr, 'Phi':Phi, 'dist':dist, 'coords':coords, 'pa':pa} for lineID, fileName, image, xScale, yScale, pixScale, JD, MJD, decYr, Phi, dist, coords, pa in zip(dic['lineID'], dic['fileName'], dic['image'], dic['xScale'], dic['yScale'], dic['pixScale'], dic['JD'], dic['MJD'], dic['decYr'], dic['Phi'], dic['dist'], dic['coords'], dic['pa']) if JD == value]

    elif (dic.has_key('profile')):

        # return the data array with profile & radPos
        if key == 'lineID':
            alist = [{'lineID':lineID, 'fileName':fileName, 'image':image, 'xScale':xScale, 'yScale':yScale, 'pixScale':pixScale, 'JD':JD, 'MJD':MJD, 'decYr':decYr, 'Phi':Phi, 'profile':profile, 'radPos':radPos} for lineID, fileName, image, xScale, yScale, pixScale, JD, MJD, decYr, Phi, profile, radPos in zip(dic['lineID'], dic['fileName'], dic['image'], dic['xScale'], dic['yScale'], dic['pixScale'], dic['JD'], dic['MJD'], dic['decYr'], dic['Phi'], dic['profile'], dic['radPos']) if lineID == value]
        if key == 'JD':
            alist = [{'lineID':lineID, 'fileName':fileName, 'image':image, 'xScale':xScale, 'yScale':yScale, 'pixScale':pixScale, 'JD':JD, 'MJD':MJD, 'decYr':decYr, 'Phi':Phi, 'profile':profile, 'radPos':radPos} for lineID, fileName, image, xScale, yScale, pixScale, JD, MJD, decYr, Phi, profile, radPos in zip(dic['lineID'], dic['fileName'], dic['image'], dic['xScale'], dic['yScale'], dic['pixScale'], dic['JD'], dic['MJD'], dic['decYr'], dic['Phi'], dic['profile'], dic['radPos']) if JD == value]

    else:

        # in this case just return the data array (no update is done)
        if key == 'lineID':
            alist = [{'lineID':lineID, 'fileName':fileName, 'image':image, 'xScale':xScale, 'yScale':yScale, 'pixScale':pixScale, 'JD':JD, 'MJD':MJD, 'decYr':decYr, 'Phi':Phi} for lineID, fileName, image, xScale, yScale, pixScale, JD, MJD, decYr, Phi in zip(dic['lineID'], dic['fileName'], dic['image'], dic['xScale'], dic['yScale'], dic['pixScale'], dic['JD'], dic['MJD'], dic['decYr'], dic['Phi']) if lineID == value]
        if key == 'JD':
            alist = [{'lineID':lineID, 'fileName':fileName, 'image':image, 'xScale':xScale, 'yScale':yScale, 'pixScale':pixScale, 'JD':JD, 'MJD':MJD, 'decYr':decYr, 'Phi':Phi} for lineID, fileName, image, xScale, yScale, pixScale, JD, MJD, decYr, Phi in zip(dic['lineID'], dic['fileName'], dic['image'], dic['xScale'], dic['yScale'], dic['pixScale'], dic['JD'], dic['MJD'], dic['decYr'], dic['Phi']) if JD == value]

    # return a list containing a dict with the data for the selected item
    return alist


def read(dataDir, dataList, velMin, velMax):

    import numpy as np
    import json
    from astropy.io import fits
    from astropy.time import Time
    import astropy.extern.six

    # reading the file with the file names
    # (must be in JSON format)
    data = json.load(open(dataList))

    # storing the file names into one single list
    fitsname = []
    for i, key in enumerate(data.keys()):
        for j, value in enumerate(data.values()[i]):
            fitsname.append(value)

    ntrans = len(data) # number of spectral lines
    nepoch = len((data.items())[0][1]) # number of observations for each line
    dic_id = np.repeat(data.keys(), nepoch) # line ID associated with each image

    JD0, period = 2456874.4, 2022.7 # this must be fixed
    raw_images = []
    images = []
    new_images = []
    dates = []
    tempJD = []
    tempMJD = []
    tempdecYr = []
    tempPhi = []

    # this is the main dictionary with n keys
    dic = {
            "lineID": [],
            "fileName": [],
            "image": [],
            "xScale": [],
            "yScale": [],
            "pixScale": [],
            "JD": [],
            "MJD": [],
            "decYr": [],
            "Phi": []
    }

    for index, item in enumerate(fitsname):

        img0, hdu0 = fits.getdata(dataDir+item, 0, header=True)

        # correcting the orientation of the images
        new_img0 = img0[:,::-1,:]

        # velocity
        vel0 = hdu0.get('VMIN') + hdu0.get('DELV') * np.asarray(range(hdu0.get('NAXIS3')))

        # date-obs
        if (item != 'zero.fits') & (item.find('4815') != -1):
            date_obs = hdu0.get('DATEOBS')
            date = Time(date_obs, scale='utc')
            Phi0 = (date.jd - JD0) / period + 13
            dates.append([date.jd,date.mjd,date.decimalyear,Phi0])
            # this will create a list with nepoch elements
            tempJD.append(date.jd)
            tempMJD.append(date.mjd)
            tempdecYr.append(date.decimalyear)
            tempPhi.append(Phi0)


        # image scale in arcsec
        cols = np.linspace(1,hdu0.get('NAXIS1'),hdu0.get('NAXIS1'))
        rows = np.linspace(1,hdu0.get('NAXIS2'),hdu0.get('NAXIS2'))
        xscale = (hdu0.get('CRPIX1') - cols) * hdu0.get('DELTA')
        yscale = (rows - hdu0.get('CRPIX2')) * hdu0.get('DELTA')

        # integrated spectrum
        # spec0 = np.sum(new_img0, axis=(1,2))

        # entire image @ selected velocity
        raw_sli00 = np.sum(new_img0[(vel0 >= velMin) & (vel0 <= velMax), :, :], axis=0).squeeze()

        # conserving flux *per pixel*
        # (https://montageblog.wordpress.com/2011/06/24/does-montage-conserve-flux-when-changing-the-image-resolution/)
        raw_sli0 = raw_sli00 * (0.1 / hdu0.get('DELTA'))**2

        # contains the entire image [index, name, data cube, xscale, yscale, pixscale, lineID]
        # raw_images.append([index, item, raw_sli0, xscale, yscale, hdu0.get('DELTA'), dic_id[index]])
        # print(dic_id[index], item, raw_sli0.shape, xscale.shape, yscale.shape, hdu0.get('DELTA'))
        dic["lineID"].append(dic_id[index])
        dic["fileName"].append(item)
        dic["image"].append(raw_sli0)
        dic["xScale"].append(xscale)
        dic["yScale"].append(yscale)
        dic["pixScale"].append(hdu0.get('DELTA'))

    # creating a list of ntrans x nepoch elements
    # that repeats itself after nepoch elements
    dic["JD"] = (np.tile(tempJD, ntrans)).tolist()
    dic["MJD"] = (np.tile(tempMJD, ntrans)).tolist()
    dic["decYr"] = (np.tile(tempdecYr, ntrans)).tolist()
    dic["Phi"] = (np.tile(tempPhi, ntrans)).tolist()

    # returns a dictionary with all the info needed
    # dic.keys = {"lineID", "fileName", "image", "xScale", "yScale", "pixScale", "JD", "MJD", "decYr","Phi"}
    return dic


if __name__ == "__main__":

    # provide a list with the name of the FITS files to be read in JSON format:
    # {
    #   "ar3": [
    #     "CUBEArIII7137_9.4_7.5.fits",
    #     "CUBEArIII7137_9.9_7.5.fits",
    #     "zero.fits",
    #     "CUBEArIII7137_11.9_7.5.fits",
    #     "CUBEArIII7137_12.8_7.5.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "CUBEArIII7137.76_2014.4_7.5x7.5-500+500_20kms_global.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "CUBEArIII7137.76_2015.0_7.5x7.5-500+500_20kms_global.fits",
    #     "zero.fits"
    #   ],
    #   "he1": [
    #     "CUBEHeI4714_9.4_7.5.fits",
    #     "CUBEHeI4714_9.9_7.5.fits",
    #     "CUBEHeI4714_10.8_7.5.fits",
    #     "CubEHeI4714_11.9_7.5.fits",
    #     "CUBEHeI4714_12.8_7.5.fits",
    #     "CUBEHeI4714_13.7_7.5.fits",
    #     "CUBEHeI4714.47_2014.2_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEHeI4714.47_2014.4_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEHeI4714.47_2014.6_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEHeI4714.47_2014.7_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEHeI4714.47_2014.9_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEHeI4714.47_2015.0_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEHeI4714.47_2015.2_7.5x7.5-500+500_20kms_global.fits"
    #   ],
    #   "fe3": [
    #     "CUBEFeIII4659_9.4_7.5.fits",
    #     "CUBEFeIII4659_9.9_7.5.fits",
    #     "CUBEFeIII4659_10.8_7.5.fits",
    #     "CUBEFeIII4659_11.9_7.5.fits",
    #     "CUBEFeIII4659_12.8_7.5.fits",
    #     "CUBEFeIII4659_13.7_7.5.fits",
    #     "CUBEFeIII4659.35_2014.2_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeIII4659.35_2014.4_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeIII4659.35_2014.6_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeIII4659.35_2014.7_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeIII4659.35_2014.9_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeIII4659.35_2015.0_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeIII4659.35_2015.2_7.5x7.5-500+500_20kms_global.fits"
    #   ],
    #   "n2": [
    #     "CUBENII5756_9.4_7.5.fits",
    #     "CUBENII5756_9.9_7.5.fits",
    #     "CUBENII5756_10.8_7.5.fits",
    #     "CUBENII5756_11.9_7.5.fits",
    #     "CUBENII5756_12.8_7.5.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "CUBENII5756.19_2014.4_7.5x7.5-500+500_20kms_global.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "CUBENII5756.19_2015.0_7.5x7.5-500+500_20kms_global.fits",
    #     "zero.fits"
    #   ],
    #   "fe2": [
    #     "CUBEFeII4815_9.4_7.5.fits",
    #     "CUBEFeII4815_9.9_7.5.fits",
    #     "CUBEFeII4815_10.8_7.5.fits",
    #     "CUBEFeII4815_11.9_7.5.fits",
    #     "CUBEFeII4815_12.8_7.5.fits",
    #     "CUBEFeII4815_13.7_7.5.fits",
    #     "CUBEFeII4815.88_2014.2_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeII4815.88_2014.4_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeII4815.88_2014.6_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeII4815.88_2014.7_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeII4815.88_2014.9_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeII4815.88_2015.0_7.5x7.5-500+500_20kms_global.fits",
    #     "CUBEFeII4815.88_2015.2_7.5x7.5-500+500_20kms_global.fits"
    #   ],
    #   "ni2": [
    #     "NiII7413_9.4_7.5_20.fits",
    #     "NiII7413_9.9_7.5_20.fits",
    #     "zero.fits",
    #     "NiII7413_11.9_7.5_20.fits",
    #     "NiII7413_12.8_7.5_20.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "CUBENiII7413.65_2014.4_7.5x7.5-500+500_20kms_global.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "zero.fits",
    #     "CUBENiII7413.65_2015.0_7.5x7.5-500+500_20kms_global.fits",
    #     "zero.fits"
    #   ]
    # }
    dataDir = '/Volumes/Kerberos/DATA/ETC/HST/TEDS_CUBE/NEW/'
    dataList = dataDir + 'data_list.json'
    # dataList = '/Volumes/Kerberos/CODES/CLOUDY/MODELS/WEIGELT/PYCLOUDY/data_list.json'
    velMin, velMax = -60, -20

    # providing dataDir & dataList independently
    # allows for having data and list stored
    # either in the same folder or in different places
    read(dataDir, dataList, velMin, velMax)
