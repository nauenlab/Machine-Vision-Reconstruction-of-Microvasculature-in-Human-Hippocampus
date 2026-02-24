import math
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow.keras.models import Model, load_model
from tensorflow_addons.losses import SigmoidFocalCrossEntropy
import skimage.morphology as morphology
from scipy.ndimage import filters, distance_transform_edt
from skimage.feature import peak_local_max
from skimage.morphology import label, remove_small_objects
from skimage.segmentation import find_boundaries, watershed
from skimage import io
import csv

def image_to_tile(stack, tile_size=128, overlap=0):
    """
    Cuts the pre-processed z-stack into 128*128 subtiles and removes blank tiles

    Parameters
    ----------
    stack :
        the pre-processed microscopy image
    tile_size : double
        tile size on x and y; default to 128
    overlap : double between [0,1]
        degree of overlap between two subtiles on both x and y direction; default to no overlap

    Returns
    -------
    tiles :
        subtiles
    """

#    tilenum_row = math.floor((stack.shape[1]-tile_size) / tile_size*(1-overlap))
#    tilenum_col = math.floor((stack.shape[2]-tile_size) / tile_size*(1-overlap))
# update 1/13/22
    tilenum_row = math.floor(stack.shape[1] / tile_size * (1 - overlap))
    tilenum_col = math.floor(stack.shape[2] / tile_size*(1-overlap))
    tiles = np.zeros([stack.shape[0]*tilenum_row*tilenum_col, tile_size, tile_size, 1])

    # add another dimension to stack
    stack = np.reshape(stack, (stack.shape[0], stack.shape[1], stack.shape[2], 1))

    counter = 0
    for i in range(stack.shape[0]): # loop through layers
        for j in range(tilenum_row):  # loop through row
            for k in range(tilenum_col): # loop through column
                curr_tile = stack[i,
                            j*tile_size*(1-overlap) : j*tile_size*(1-overlap)+tile_size,
                            k*tile_size*(1-overlap) : k*tile_size*(1-overlap)+tile_size, :]
                tiles[counter,:,:, :] = curr_tile
                counter = counter + 1

    return tiles

def unet2d_prediction(tiles, model_fname, img_shape, threshold, tile_size=128, overlap=0):
    """
    Makes binary segmentation on cut tiles using trained weights

    Parameters
    ----------
    tiles : ndarray
        the cut tiles, normalized to within [0,1]
    model_fname : string
        file name of the unet2d model and trained weights
    img_shape : array
        shape of the original image (layer, row, column)
    threshold : double
        binarization threshould after prediction
    tile_size : double
        size of each tile in xy direction; default to 128
    overlap : double within [0,1]
        degree of overlap when cutting the tiles; default to 0

    Returns
    -------
    preds : ndarray
        prediction
    """

    # unet2d prediction
    print("[INFO] Loading best model:")
    model = load_model(model_fname)
    print("[INFO] Using the best model to make predictions on testing set:")
    preds_tiles = model.predict(tiles, verbose=1)
    preds_tiles = preds_tiles > threshold

    # stitch prediction and blank tiles together
    print("[INFO] stitching tiles")
    preds = np.zeros((img_shape[0], img_shape[1], img_shape[2], 1))
    counter = 0
#    tilenum_row = math.floor((stack.shape[1]-tile_size) / tile_size*(1-overlap))
#    tilenum_col = math.floor((stack.shape[2]-tile_size) / tile_size*(1-overlap))
# update 1/13/22
    tilenum_row = math.floor(img_shape[1] / tile_size * (1 - overlap))
    tilenum_col = math.floor(img_shape[2] / tile_size * (1 - overlap))
    for i in range(img_shape[0]): # loop through layers
        for j in range(tilenum_row):  # loop through row
            for k in range(tilenum_col): # loop through column
                if counter >= preds_tiles.shape[0]:
                    break
                preds[i,
                    j * tile_size * (1 - overlap): j * tile_size * (1 - overlap) + tile_size,
                    k * tile_size * (1 - overlap): k * tile_size * (1 - overlap) + tile_size, :] = preds_tiles[counter]
                counter = counter + 1

    return preds

def unet2d_predNoThreshold(tiles, model_fname, img_shape, tile_size=128, overlap=0):
    """
    Makes prediction (0-1; possibility, no threshold) on cut tiles using trained weights

    Parameters
    ----------
    tiles : ndarray
        the cut tiles (normalized to within [0,1])
    model_fname : string
        file name of the unet2d model and trained weights
    img_shape : array
        shape of the original image (layer, row, column)
    tile_size : double
        size of each tile in xy direction; default to 128
    overlap : double within [0,1]
        degree of overlap when cutting the tiles; default to 0

    Returns
    -------
    preds : ndarray
        prediction
    """

    # unet2d prediction
    print("[INFO] Loading best model:")
    model = load_model(model_fname)
    print("[INFO] Using the best model to make predictions on testing set:")
    preds_tiles = model.predict(tiles, verbose=1)

    # stitch prediction and blank tiles together
    print("[INFO] stitching tiles")
    preds = np.zeros((img_shape[0], img_shape[1], img_shape[2], 1))
    counter = 0
#    tilenum_row = math.floor((img_shape[1] - tile_size) / tile_size * (1 - overlap))
#    tilenum_col = math.floor((img_shape[2] - tile_size) / tile_size * (1 - overlap))
# Edit 1/13/22
    tilenum_row = math.floor((img_shape[1]) / tile_size * (1 - overlap))
    tilenum_col = math.floor((img_shape[2]) / tile_size * (1 - overlap))
    for i in range(img_shape[0]): # loop through layers
        for j in range(tilenum_row):  # loop through row
            for k in range(tilenum_col): # loop through column
                if counter >= preds_tiles.shape[0]:
                    break
                preds[i,
                    j * tile_size * (1 - overlap): j * tile_size * (1 - overlap) + tile_size,
                    k * tile_size * (1 - overlap): k * tile_size * (1 - overlap) + tile_size, :] = preds_tiles[counter]
                counter = counter + 1

    return preds

def green_bubble_removal(bn_mask, min_size, assign_diff_labels=False):
    """
        Removes green bubbles/nuclei on capillary segmentation outputs (assign either value of
        0 or 2), based on size exclusion.

        Parameters
        ----------
        bn_mask : ndarray
            the binary mask after capillary segmentation
        min_size : double
            the minimum number of voxels below which this structure will be discarded; default to 7500
        assign_diff_labels : bool
            if true, kept structures will be labeled as 2, and removed structures will be labeled as as 1
            if false, kept will be 1, removed (and background) will both be 0

        Returns
        -------
        output_mask : ndarray
            the output mask after green bubble removal
        """

    # get all connected regions in 3d
    labeled_mask, label_num = label(bn_mask, return_num=True)

    # remove small regions
    output_mask = remove_small_objects(labeled_mask, min_size=min_size, connectivity=3)
    output_mask[output_mask > 0] = 1

    # optionally, not remove small regions, but assigned differnt labels to then
    if assign_diff_labels:
        output_mask = bn_mask + output_mask

    output_mask = output_mask.astype(np.float64)

    return output_mask

def watershed_2d(bn_mask, sigma, min_distance):
    """
    Segment cells in each layer of the 3D image by 2D _watershed
    Author: Chentao Wen
    Edits:
        Added erosion after taking distance map
        Added handle to manipulate degree of gaussian blur

    Parameters
    ----------
    bn_mask :
        the binary image of cell region and background (predicted by 2D U-net)
    sigma :
        the degree of gaussian blur
    min_distance :
        the minimum cell distance allowed in the result

    Returns
    -------
    bn_output :
        binary image (cell/bg) removing boundaries detected by _watershed
    boundary :
        image of cell boundaries
    """

    print("[INFO] performing 2d watershed")

    bn_mask = np.squeeze(bn_mask)
    boundary = np.zeros(bn_mask.shape, dtype='bool')

    for i in range(bn_mask.shape[0]):
        # convert binary segmentation into logic data type
        bn_image = bn_mask[i, :, :] > 0.5*255
        # calculate distance map (euclidean distance)
        dist = distance_transform_edt(bn_image, sampling=[1, 1])
        # erode the distance map to make dendrites disappear
        dist_eroded = morphology.erosion(dist)
        # apply gaussian filter to blur out details
        dist_smooth = filters.gaussian_filter(dist_eroded, sigma=sigma, mode='constant')

        # find local maximum in the dist map within neighborhood of min_distance
        local_maxi = peak_local_max(dist_smooth, min_distance=min_distance, indices=False)
        # local_maxi = peak_local_max(dist_smooth, footprint=np.ones((7,7)), indices=False)
        # assign one label for all connected local maximum
        markers, num_labels = morphology.label(local_maxi, return_num=True)
        # 2d seeded watershed
        labels_ws = watershed(-dist_smooth, markers, mask=bn_image)
        # find the boundary of the watershed image
        labels_bd = find_boundaries(labels_ws, connectivity=2, mode='outer', background=0)
        # save to boundary
        boundary[i, :, :] = labels_bd

        print(i)

    # remove boundary from watershed image
    bn_output = bn_mask > 0.5
    bn_output[boundary == 1] = 0

    return bn_output, boundary

def watershed_3d(image_watershed2d, samplingrate, sigma_3d, min_distance, min_size):
    """
    Segment cells by 3D _watershed
    Author: Chentao Wen
    Edits:
        Added 3d morphological erosion
        Added handle to manipulate the degree of gaussian blur
        Got rid of methods, min_size, cell_num --> no structure is removed here

    Parameters
    ----------
    image_watershed2d :
        the binary image (cell/bg) obtained by watershed_2d
    samplingrate : list
        resolution in x, y, and z axis to calculate 3D distance
    sigma_3d : list
        the degree of gaussian blur on 3 dimensions
    min_size :
        minimum size of neurons (in pixels)

    Returns
    -------
    labels_wo_bd :
        label image of cells removing boundaries (set to 0)
    labels_clear :
        label image of cells before removing boundaries

    Notes
    -----
    For peak_local_max function, exclude_border=0 is important. Without it, the function will exclude the cells
    within bottom/top layers (<=min_distance layers)
    """
    # 3d distance map
    print("3d distance map")
    dist = distance_transform_edt(image_watershed2d, sampling=samplingrate)
    # erode the distance map to make dendrites disappear
    print("erode")
    dist_eroded = morphology.erosion(dist)
    # gaussian filter to blur out distance map
    print("gauss")
    dist_smooth = filters.gaussian_filter(dist_eroded, sigma_3d, mode='constant')
    # find local max in distance map (local max separated by min_distance
    print("local max")
    local_maxi = peak_local_max(dist_smooth, min_distance=min_distance, exclude_border=0, indices=False)
    # condense local max if they are connected + label local max
    print("condense local max")
    markers = morphology.label(local_maxi)
    # 3d watershed
    print("3d watershed")
    labels_ws = watershed(-dist_smooth, markers, mask=image_watershed2d)
    # assign labels
    labels_clear = remove_small_objects(labels_ws, min_size=min_size, connectivity=3)

    labels_bd = find_boundaries(labels_clear, connectivity=3, mode='outer', background=0)
    labels_wo_bd = labels_clear.copy()
    labels_wo_bd[labels_bd == 1] = 0
    labels_wo_bd = remove_small_objects(labels_wo_bd, min_size=min_size, connectivity=3)

    return labels_wo_bd, labels_clear

def neuron_features(image_watershed3d, px_res):
    """
    This function takes the serial-numbered object instance segmentation to calculate neuron features
    Parameters
    ----------
    image_watershed3d :
        Object instance segmentation mask in 3D, where each object is numbered
    px_res : double
        microns per pixel

    Returns
    -------
    bbox :
        list; (z_min, z_max, x_min, x_max, y_min, y_max) for each object instance
    centroid :
        list; the z,x,y location of the centroid pixel
    obj_volume :
        list; volume of each object instance in um^3
    neurons :
        list; m*3 ndarray with z,x,y location of each neuron pixel, and m == object volume
    neuron_label :
        list; the label of a neuron that satisfies vol requirement
    """

    # total number of cells identified
    num = np.max(image_watershed3d)

    # initialize return values
    bbox = [None] * num
    centroid = [None] * num
    obj_volume = [None] * num
    neurons = [None] * num
    neuron_label = [None] * num

    # Go through each object instance (only analyze those that satisfy size requirement)
    for i in range(num):
        # extract pixel (x,y,x) location with label==i+1
        neuron = np.argwhere(image_watershed3d == i + 1)

        if neuron.shape[0] > 0:
            # neuron pixel location info
            neurons[i] = neuron

            # bbox info
            bbox[i] = [np.min(neuron[:, 0]),
                       np.max(neuron[:, 0]),
                       np.min(neuron[:, 1]),
                       np.max(neuron[:, 1]),
                       np.min(neuron[:, 2]),
                       np.max(neuron[:, 2])]

            # centroid info
            centroid[i] = [np.mean(neuron[:, 0]),
                           np.mean(neuron[:, 1]),
                           np.mean(neuron[:, 2])]

            # object volume info
            obj_volume[i] = neuron.size * px_res

            # label of that neuron
            neuron_label[i] = i+1

            print("Now processing neuron #" + str(i))

    return bbox, centroid, obj_volume, neurons, neuron_label