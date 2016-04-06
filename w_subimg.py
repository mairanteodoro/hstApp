#!/usr/bin/env python

def main(raw_image, pixscale, coord1_r, coord1_c, coord2_r, coord2_c):

    #
    # RETURN A SUB-IMAGE OF A GIVEN IMAGE
    #

    import numpy as np

    # START OF CRITICAL PART
    # coordinates array
    # entire image
    im_noise = raw_image
    r0r, c0r = 0, 0
    # row indexes increase in the same direction as labels
    tick_r = -((im_noise.shape[0]-1)/2. - np.linspace(-0.5,im_noise.shape[0]-0.5,im_noise.shape[0])) * pixscale
    # row indexes increase in the opposite direction of increasing labels
    tick_c = ((im_noise.shape[1]-1)/2. - np.linspace(-0.5,im_noise.shape[1]-0.5,im_noise.shape[1])) * pixscale
    # the coordinates of the center of each pixel are stored here:
    tick_r_array = np.around((tick_r[0] + pixscale/2) + (np.array(range(r0r,r0r+tick_r.shape[0])) - r0r) * pixscale, decimals=2)
    tick_c_array = np.around((tick_c[0] - pixscale/2) - (np.array(range(c0r,c0r+tick_c.shape[0])) - c0r) * pixscale, decimals=2)

    # print('Extracting the coordinates betwee {:+0.2f} and {:+0.2f}:'.format(coord1,coord2))
    # print(tick_r_array[np.where((tick_r_array >= coord1) & (tick_r_array <= coord2))])

    # - difference between the coord arrays and the desired coords
    # - round up to the 2 decimal place
    dif0 = np.around(tick_r_array - coord1_r, decimals=2)
    dif1 = np.around(tick_r_array - coord2_r, decimals=2)
    dif2 = np.around(tick_c_array - coord1_c, decimals=2)
    dif3 = np.around(tick_c_array - coord2_c, decimals=2)

    # - now, the absolute value of the difference
    abs_dif0 = abs(dif0)
    abs_dif1 = abs(dif1)
    abs_dif2 = abs(dif2)
    abs_dif3 = abs(dif3)

    # now, find out where the difference is minimum
    # print('Sel indexes:')
    # selecting the indexes corresponding to the coords
    indx_r0 = (np.where(dif0 == np.min(abs_dif0)))[0]
    # print(indx_r0)
    indx_r1 = (np.where(dif1 == np.min(abs_dif1)))[0] #+ 1.0
    # print(indx_r1)
    indx_c0 = (np.where(dif2 == np.min(abs_dif2)))[0]
    # print(indx_c0)
    indx_c1 = (np.where(dif3 == np.min(abs_dif3)))[0] #+ 1.0
    # print(indx_c1)
    # np.around([0.37, 1.64], decimals=1)
    # print(np.around(abs(tick_c_array - coord2_c),decimals=2), np.around(np.min(abs(tick_c_array - coord2_c)),decimals=2))
    # remmember to add 1 to include the entire array
    # indx_r0, indx_r1 = indx_r[0][0], indx_r[0][-1] + 1
    # indx_c0, indx_c1 = indx_c[0][0], indx_c[0][-1] + 1

    row1, row2 = int(((sorted([indx_r0,indx_r1]))[0]).squeeze()), int(((sorted([indx_r0,indx_r1]))[-1]).squeeze())
    col1, col2 = int(((sorted([indx_c0,indx_c1]))[0]).squeeze()), int(((sorted([indx_c0,indx_c1]))[-1]).squeeze())
    # print(row1,row2,col1,col2)
    # print([tick_c_array[col1]+0.1/2.,tick_c_array[col2]-0.1/2.,tick_r_array[row2]+0.1/2.,tick_r_array[row1]-0.1/2.])
    subim_r = im_noise[row1:row2+1,col1:col2+1]
    # extent_r = [tick_cr[0],tick_cr[-1],tick_rr[-1],tick_rr[0]]
    extent_r = [tick_c_array[col1]+pixscale/2.,tick_c_array[col2]-pixscale/2.,
                tick_r_array[row2]+pixscale/2.,tick_r_array[row1]-pixscale/2.]

    return [subim_r, extent_r, row1, row2, col1, col2]

    # END OF CRITICAL PART

if __name__ == "__main__":

    main(raw_image,pixscale,coord1_r,coord1_c,coord2_r,coord2_c)
