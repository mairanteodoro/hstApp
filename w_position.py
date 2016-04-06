#!/usr/bin/env python

# POSITION OF THE BLOBS

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mairan_fix_bbox.fix_bbox import fix_bbox
from skimage import measure
from skimage import feature


def get_centroid(data):

    h,w = np.shape(data)
    x = np.arange(0,w)
    y = np.arange(0,h)

    X,Y = np.meshgrid(x,y)

    cx = np.sum(X * data) / np.sum(data)
    cy = np.sum(Y * data) / np.sum(data)

    return cx,cy


def main(raw_images,ntrans,Phi,coords_r,coords_c):

    import w_subimg

    lines = ntrans * ['line_id']
    epoch = Phi
    nepoch = len(Phi)

    # centroid position for each transition and epoch and blob
    # (middle col refers to the blobs and last col refers to x,y)
    centroid_array = np.zeros((ntrans*nepoch,2,2))
    flux_array = np.zeros((ntrans*nepoch,2))

    for j in range(2):

        coord1_r, coord2_r = coords_r[2*j:2*(j+1)]
        coord1_c, coord2_c = coords_c[2*j:2*(j+1)]
        for i in range(ntrans*nepoch):

            if raw_images[i][1] == 'zero.fits':
                raw_images[i][2] = raw_images[i][2] * 0

            # call the sub-image function
            subimg = w_subimg.main(raw_images[i][2],raw_images[i][5],coord1_r,coord1_c,coord2_r,coord2_c)

            subim_r, extent = subimg[0:2]
            row1, row2, col1, col2 = subimg[2:]

            # extracting the subimage
            img = subim_r.astype('double')

            if raw_images[i][1] != 'zero.fits':
                # pos = feature.blob_dog(img, threshold=1e-14, min_sigma=1.5, max_sigma=2.0, overlap=0.5)
                pos = feature.peak_local_max(img, threshold_abs=1e-14, min_distance=1)
                mom = measure.moments(img)
                '''
                The following properties can be calculated from raw image moments:
                Area as: m[0, 0].
                Centroid as: {m[0, 1] / m[0, 0], m[1, 0] / m[0, 0]}.
                '''
                centroid = [mom[0,1] / mom[0,0], mom[1,0] / mom[0,0]]
                # centroid2 = get_centroid(img)
                # center of img
                center = (((raw_images[i][2]).shape)[0] - 1) / 2
                cent_col = (center - (col1+centroid[1])) * raw_images[i][5]
                cent_row = (row1+centroid[0] - center) * raw_images[i][5]
                centroid_array[i,j,:] = [cent_row,cent_col]
                # *REMOVING* THE FLUX CONSERVATION *PER PIXEL* CORRECTION FACTOR
                flux_array[i,j] = np.sum(img) * (raw_images[i][5] / 0.1)**2
                if pos.shape[0] != 0:
                    pos_col = (center - (col1+pos[:,1]).flatten()) * raw_images[i][5]
                    pos_row = ((row1+pos[:,0]).flatten() - center) * raw_images[i][5]

    return centroid_array



if __name__ == "__main__":

    # THIS MODULE CAN BE CALLED IN TWO WAYS:
    # (1) ./w_flux.py (in which case the code will run with the defaul parameters listed below)
    # (2) from another script by providing the following arguments:
    #
    # --- from here
    #
    # defining the corners of the area comprising the blobs
    # Weigelt C
    coord1_r_WC, coord2_r_WC =  -0.2, +0.3 # -1, +1
    coord1_c_WC, coord2_c_WC =  -0.5, -0.1 # -1, +1
    # Weigelt D
    coord1_r_WD, coord2_r_WD =  -0.2, -0.8 # -1, +1
    coord1_c_WD, coord2_c_WD =  -0.2, +0.3 # -1, +1
    # concatenating
    coords_r = [coord1_r_WC, coord2_r_WC, coord1_r_WD, coord2_r_WD]
    coords_c = [coord1_c_WC, coord2_c_WC, coord1_c_WD, coord2_c_WD]

    import hstData

    velMin, velMax = -60, -20
    data = hstData.read(dataDir, dataList, velMin, velMax)

    raw_images = data[0]
    ntrans = len(data[0]) / len(data[1][0]) # = number of images / number of JD entries
    Phi = data[1][3]
    epoch = Phi
    nepoch = len(Phi)

    centroid_array = main(raw_images,ntrans,Phi,coords_r,coords_c)

    ### --- up to here



    ### BELOW IS THE GRAPHIC PART

    # DISTANCE FROM CENTRAL SOURCE
    leg_fontsize = 10
    ion = [u'Ar$^{2+}$',u'He$^0$',u'Fe$^{2+}$',u'N$^{1+}$',u'Fe$^{1+}$',u'Ni$^{1+}$']
    markers = ['D', '+', 'o', '^', 's', 'x']
    colors  = ['m', 'b', 'g', 'c', 'orange', 'r']

    f, axarr = plt.subplots(2, 2, figsize=(10,9), sharex='col')
    axarr2 = axarr.flatten()

    axarr2[0].set_ylabel('Projected distance from the central source (arcsec)')
    axarr2[0].set_ylim(0.15,0.30)
    axarr2[1].set_ylim(0.30,0.45)

    axarr2[2].set_ylabel('Position angle')
    axarr2[2].set_ylim(260,285)
    axarr2[3].set_ylim(340,365)

    yticks1 = axarr2[2].get_yticks()
    ylabel1 = []
    for x in yticks1:
        ylabel1.append(u'{0}$^\circ$'.format(int(x)))
    axarr2[2].yaxis.set_ticklabels(ylabel1)

    yticks2 = axarr2[3].get_yticks()
    ylabel2 = []
    for x in yticks2:
        ylabel2.append(u'{0}$^\circ$'.format(int(x)))
    axarr2[3].yaxis.set_ticklabels(ylabel2)

    axarr2[0].set_title("Weigelt C")
    axarr2[1].set_title("Weigelt D")

    for j in range(2):

        for i in range(ntrans):

            x = centroid_array[i*13:(i+1)*13,j,1] # cols
            y = centroid_array[i*13:(i+1)*13,j,0] # rows

            dist = np.sqrt(x**2 + y**2)
            pa = (2. * np.pi - np.arctan(x/y)) / np.pi * 180.

            axarr2[j].yaxis.set_major_locator(MaxNLocator(4))
            axarr2[j].plot(epoch[dist>0], dist[dist>0],
                           label='{}'.format(ion[i]),
                           linestyle='None',
                           marker=markers[i],
                           color=colors[i])
            axarr2[j].axvline(x=13, color='0.5', ls='--', marker=' ', zorder=0, lw=1)

            axarr2[j+2].plot(epoch[pa>0], pa[pa>0],
                           label='{}'.format(ion[i]),
                           linestyle='None',
                           marker=markers[i],
                           color=colors[i])
            axarr2[j+2].axvline(x=13, color='0.5', ls='--', marker=' ', zorder=0, lw=1)

            if j == 0:
                handles, labels = axarr2[j].get_legend_handles_labels()
                axarr2[j].legend(handles, labels, loc='best', ncol=2, numpoints=1, frameon=True, fontsize=leg_fontsize)

    plt.savefig('panel_pos2.eps')
    fix_bbox('panel_pos2.eps')
    plt.close(f)
