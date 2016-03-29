#!/usr/bin/env python

# extract the radial profile along a given direction (radius,PA)
def profile(data, radius, posAngle, slitWidth):

    import numpy as np
    import w_subimg
    import copy
    from skimage import measure
    from skimage import feature

    ntrans = len(np.unique(data['lineID']))
    nepoch = len(np.unique(data['JD']))
    Phi = np.unique(data['Phi'])

    # http://www.python-course.eu/passing_arguments.php
    # https://jeffknupp.com/blog/2012/11/13/is-python-callbyvalue-or-callbyreference-neither/
    dic = copy.deepcopy(data) # <- THIS IS VERY IMPORTANT! KEEP copy.deepcopy()!
    dic['profile'] = [] # new key: line radial profile
    dic['radPos'] = [] # new key: positions where the radial profile was measured

    raw_images = [
                    data['fileName'],
                    data['image'],
                    data['pixScale']
                ]

    for j in range(ntrans):

        for i in range(nepoch):

            index = i + j * nepoch

            # print(j, i, index, raw_images[0][index], coords_row, coords_col)

            if raw_images[0][index] == 'zero.fits':
                raw_images[1][index] = raw_images[1][index] * 0

            img = copy.deepcopy(raw_images[1][index])

            rho = radius # length in arcsec
            # theta=0 @ 6:00; theta=90 @ 3:00; theta=180 @ 12:00
            theta = (posAngle + 180. - 360.) if (posAngle >= 180.) else (posAngle + 180.)
            rho_x = rho * np.cos((theta) / 180. * np.pi)
            rho_y = rho * np.sin((theta) / 180. * np.pi)

            actual_image = img # transf.resize(raw_images[j][2], raw_images[0][2].shape)

            # center of image (same for X and Y)
            x0 = ((img.shape)[0] - 1) / 2
            dx = rho_x / raw_images[2][index]
            dy = rho_y / raw_images[2][index]
            x1 = x0 + dx
            y1 = x0 + dy

            # slit width
            lineWidth = slitWidth / raw_images[2][index]

            xstart, ystart = x0, x0
            xend, yend = x1, y1

            if j == 0:
                xstart0, ystart0, xend0, yend0 = xstart, ystart, xend, yend

            # print(j,(xstart,ystart), (xend,yend))

            profile = measure.profile_line(actual_image, (xstart,ystart), (xend,yend), linewidth=lineWidth)
            xvar = np.linspace(0,rho,profile.shape[0])

            dic["profile"].append(profile)
            dic["radPos"].append(xvar)

    return dic


# measure the centroid of a defined region (for now, a box)
def position(data, coords_row, coords_col):

    import numpy as np
    import w_subimg
    import copy
    from skimage import measure
    from skimage import feature

    ntrans = len(np.unique(data['lineID']))
    nepoch = len(np.unique(data['JD']))
    Phi = np.unique(data['Phi'])

    # http://www.python-course.eu/passing_arguments.php
    # https://jeffknupp.com/blog/2012/11/13/is-python-callbyvalue-or-callbyreference-neither/
    dic = copy.deepcopy(data) # <- THIS IS VERY IMPORTANT! KEEP copy.deepcopy()!
    dic['coords'] = [] # new key: centroid coordinates (x,y)
    dic['dist'] = [] # new key: centroid position
    dic['pa'] = [] # new key: position angle E from N

    raw_images = [
                    data['fileName'],
                    data['image'],
                    data['pixScale']
                ]

    for j in range(ntrans):

        for i in range(nepoch):

            index = i + j * nepoch

            # print(j, i, index, raw_images[0][index], coords_row, coords_col)

            if raw_images[0][index] == 'zero.fits':
                raw_images[1][index] = raw_images[1][index] * 0

            # call the sub-image function
            subimg = w_subimg.main(raw_images[1][index],
                                   raw_images[2][index],
                                   coords_row[0], coords_col[0],
                                   coords_row[1], coords_col[1])

            subim_r, extent = subimg[0:2]
            row1, row2, col1, col2 = subimg[2:]

            # extracting the subimage
            img = subim_r.astype('double')

            mom = measure.moments(img)
            '''
            The following properties can be calculated from raw image moments:
            Area as: m[0, 0].
            Centroid as: {m[0, 1] / m[0, 0], m[1, 0] / m[0, 0]}.
            '''
            centroid = [mom[0,1] / mom[0,0], mom[1,0] / mom[0,0]]
            # centroid2 = get_centroid(img)
            # center of img
            center = (((raw_images[1][index]).shape)[0] - 1) / 2
            cent_col = (center - (col1+centroid[1])) * raw_images[2][index]
            cent_row = (row1+centroid[0] - center) * raw_images[2][index]

            # position angle measured E from N
            dist = np.sqrt(cent_col**2 + cent_row**2)
            pa = (2. * np.pi - np.arctan(cent_col/cent_row)) / np.pi * 180.

            dic["dist"].append(dist)
            dic["coords"].append([cent_row,cent_col])
            dic["pa"].append(pa)

    return dic


# measure the flux inside a defined region (for now, a box)
def flux(data, coords_row, coords_col):

    import numpy as np
    import w_subimg
    import copy

    ntrans = len(np.unique(data['lineID']))
    nepoch = len(np.unique(data['JD']))
    Phi = np.unique(data['Phi'])

    # http://www.python-course.eu/passing_arguments.php
    # https://jeffknupp.com/blog/2012/11/13/is-python-callbyvalue-or-callbyreference-neither/
    dic = copy.deepcopy(data) # <- THIS IS VERY IMPORTANT! KEEP copy.deepcopy()!
    dic['flux'] = [] # new key

    raw_images = [
                    data['fileName'],
                    data['image'],
                    data['pixScale']
                ]

    for j in range(ntrans):

        for i in range(nepoch):

            index = i + j * nepoch

            # print(j, i, index, raw_images[0][index], coords_row, coords_col)

            if raw_images[0][index] == 'zero.fits':
                raw_images[1][index] = raw_images[1][index] * 0

            # call the sub-image function
            subimg = w_subimg.main(raw_images[1][index],
                                   raw_images[2][index],
                                   coords_row[0], coords_col[0],
                                   coords_row[1], coords_col[1])

            subim_r, extent = subimg[0:2]
            row1, row2, col1, col2 = subimg[2:]

            # extracting the subimage
            img = subim_r.astype('double')

            dic["flux"].append(np.sum(img) * (raw_images[2][index] / 0.1)**2)

            # print(dic["flux"][0])

    return dic



if __name__ == "__main__":

    # THIS MODULE CAN BE CALLED IN TWO WAYS:
    # (1) ./hstMeasure.py (in which case the code will run with the defaul parameters listed below)
    # (2) from another script by providing the following arguments:
    #
    # --- from here
    #

    import numpy as np
    import copy
    import hstData, hstMeasure

    dataDir = '/Volumes/Kerberos/DATA/ETC/HST/TEDS_CUBE/NEW/'
    dataList = dataDir + 'data_list.json'
    velMin, velMax = -60, -20
    data = hstData.read(dataDir, dataList, velMin, velMax)

    ntrans = len(np.unique(data['lineID']))
    nepoch = len(np.unique(data['JD']))
    Phi = np.unique(data['Phi'])

    # defining the corners of the area comprising the blobs
    # Weigelt C
    coord1_r_WC, coord2_r_WC =  -0.2, +0.3 # -1, +1
    coord1_c_WC, coord2_c_WC =  -0.5, -0.1 # -1, +1
    # Weigelt D
    coord1_r_WD, coord2_r_WD =  -0.2, -0.8 # -1, +1
    coord1_c_WD, coord2_c_WD =  -0.2, +0.3 # -1, +1
    # concatenating
    coords_r_WC = [coord1_r_WC, coord2_r_WC]
    coords_c_WC = [coord1_c_WC, coord2_c_WC]
    coords_r_WD = [coord1_r_WD, coord2_r_WD]
    coords_c_WD = [coord1_c_WD, coord2_c_WD]

    flux_WC = hstMeasure.flux(data, coords_r_WC, coords_c_WC)
    flux_WD = hstMeasure.flux(data, coords_r_WD, coords_c_WD)

    pos_WC = hstMeasure.position(data, coords_r_WC, coords_c_WC)
    pos_WD = hstMeasure.position(data, coords_r_WD, coords_c_WD)


    # plot the result with something like this:
    # he1_WC = hstData.get(flux_WC, key='lineID', value='he1')
    # he1_WC_flux = [he1_WC[i]['flux'] for i in range(nepoch)]
    # he1_WD = hstData.get(flux_WD, key='lineID', value='he1')
    # he1_WD_flux = [he1_WD[i]['flux'] for i in range(nepoch)]
    # x = [he1_WC[i]['Phi'] for i in range(nepoch)] or x = np.unique(he1_WC['Phi'])
    # plt.plot(x, he1_WC_flux)
    # plt.plot(x, he1_WD_flux)

    ### --- up to here



    ### BELOW IS THE GRAPHIC PART

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from mairan_fix_bbox.fix_bbox import fix_bbox

    # FLUX
    leg_fontsize = 10
    epoch = Phi

    lineID = [u'ar3', u'he1', u'fe3', u'n2', u'fe2', u'ni2']
    line = [u'[Ar\,{\sc iii}]~$\lambda7137$', u'He\,{\sc i}~$\lambda4714$',
            u'[Fe\,{\sc iii}]~$\lambda4659$', u'[N\,{\sc ii}]~$\lambda5756$',
            u'[Fe\,{\sc ii}]~$\lambda4815$', u'[Ni\,{\sc ii}]~$\lambda7413$']

    markers = ['D', 'o']
    colors  = ['k', 'r']

    legentTextID = ['Weigelt C', 'Weigelt D']

    # X and Y for plotting
    ionData_flux1 = np.empty((len(lineID), nepoch))
    ionData_flux2 = np.empty((len(lineID), nepoch))
    varx = np.unique(data['Phi'])

    for i in range(len(lineID)):
        # the order of the list entries is defined by the order in var 'lineID' above
        ionData1 = hstData.get(flux_WC, key='lineID', value=lineID[i]) # retrieving the data for the ion
        ionData_flux1[i,:] = [ionData1[x]['flux'] for x in range(nepoch)] # storing the flux for each epoch

        ionData2 = hstData.get(flux_WD, key='lineID', value=lineID[i]) # retrieving the data for the ion
        ionData_flux2[i,:] = [ionData2[x]['flux'] for x in range(nepoch)] # storing the flux for each epoch



    f, axarr = plt.subplots(3, 2, figsize=(10,7), sharex=True, sharey=True)

    for i in range(len(axarr)):

        nonZero1 = ionData_flux1[2*i] > 0
        nonZero2 = ionData_flux2[2*i+1] > 0

        # left column
        axarr[i,0].plot(varx[nonZero1], np.log10(ionData_flux1[2*i][nonZero1]), color=colors[0], marker=markers[0])
        axarr[i,0].plot(varx[nonZero1], np.log10(ionData_flux2[2*i][nonZero1]), color=colors[1], marker=markers[1])
        # right column
        axarr[i,1].plot(varx[nonZero2], np.log10(ionData_flux1[2*i+1][nonZero2]), color=colors[0], marker=markers[0])
        axarr[i,1].plot(varx[nonZero2], np.log10(ionData_flux2[2*i+1][nonZero2]), color=colors[1], marker=markers[1])

        # ionData1 = hstData.get(flux_WC, key='lineID', value=lineID[])
    # axarr2 = axarr.flatten()
    # axarr2_pos_indx = [0,1,4,5,8,9,2,3,6,7,10,11] # <- what's this?
    # counter = 0
    # for k in range(2):
    #
    #     for i in range(ntrans):
    #
    #         axarr2[axarr2_pos_indx[counter]].xaxis.set_major_locator(MaxNLocator(4))
    #         axarr2[axarr2_pos_indx[counter]].yaxis.set_major_locator(MaxNLocator(4))
    #
    #         for j in range(2):
    #
    #             if k == 0:
    #                 # flux
    #                 flux = flux_array[i*nepoch:(i+1)*nepoch,j]
    #                 # print(flux[flux>0])
    #                 axarr2[axarr2_pos_indx[counter]].set_ylim(np.log10(1e-13),np.log10(1e-10))
    #                 # print(np.min(flux[flux>0]),np.max(flux[flux>0]))
    #                 axarr2[axarr2_pos_indx[counter]].set_title(line[i])
    #                 # axarr2[axarr2_pos_indx[counter]].xaxis.set_major_locator(MaxNLocator(4))
    #                 # axarr2[axarr2_pos_indx[counter]].yaxis.set_major_locator(MaxNLocator(4))
    #                 axarr2[axarr2_pos_indx[counter]].plot(epoch[flux>0], np.log10(flux[flux>0]),
    #                                label='{}'.format(legentTextID[j]),
    #                                linestyle='None',
    #                                marker=markers[j],
    #                                color=colors[j])
    #             else:
    #                 # flux normalized at apastron
    #                 flux = flux_array[i*nepoch:(i+1)*nepoch,j] / flux_array[i*nepoch+3,j] # <-- ERROR HERE
    #                 axarr2[axarr2_pos_indx[counter]].set_ylim(0.0,1.5)
    #                 axarr2[axarr2_pos_indx[counter]].set_title(line[i])
    #                 # axarr2[axarr2_pos_indx[counter]].xaxis.set_major_locator(MaxNLocator(4))
    #                 # axarr2[axarr2_pos_indx[counter]].yaxis.set_major_locator(MaxNLocator(3))
    #                 axarr2[axarr2_pos_indx[counter]].plot(epoch[flux>0], flux[flux>0],
    #                                label='{}'.format(legentTextID[j]),
    #                                linestyle='None',
    #                                marker=markers[j],
    #                                color=colors[j])
    #
    #             axarr2[axarr2_pos_indx[counter]].axvline(x=13, color='0.5', ls='--', marker=' ', zorder=0, lw=1)
    #
    #             if k == 0 and i == 1:
    #                 handles, labels = axarr2[axarr2_pos_indx[counter]].get_legend_handles_labels()
    #                 axarr2[axarr2_pos_indx[counter]].legend(handles, labels, loc='best', ncol=1, numpoints=1, frameon=True, fontsize=leg_fontsize)
    #
    #         counter+=1

    plt.savefig('panel_flux2.eps')
    fix_bbox('panel_flux2.eps')
    plt.close(f)
