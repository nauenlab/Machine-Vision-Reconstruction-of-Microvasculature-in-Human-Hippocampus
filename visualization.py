import numpy as np
import matplotlib.pyplot as plt
import tifffile as tiff
import hippoReconstruction as hippo

# files
img_fname = ["good_B5_CA1_img.tif",
             "greenBubbles_C6_CA2_img.tif",
             "shadowyNeurons_B5_CA3_img.tif",
             "specklyBkg_C1_CA1_img.tif"]
model_fname = ["condition1_best-model.h5",
               "condition2_with_best-model.h5",
               "condition2_without_best-model.h5"]

for i in range(len(model_fname)):
    for j in range(len(img_fname)):
        # load testing images
        print("[INFO] loading testing images")
        img = tiff.imread(img_fname[j]) # tiff encoding: (depth, x, y)
        img = img / 255     # normalize 8-bit image to within [0,1]
        img = img[0:100,:,:]    # only taking the first 100 layers for the sake of time

        # Making predictions on testing images
        print("[INFO] cutting zstack images to tiles")
        tiles = hippo.image_to_tile(img)
        print("[INFO] network making binary segmentation prediction")
        pred = hippo.unet2d_prediction(tiles, model_fname[i], img.shape, 0.5) # (depth,x,y)
        pred = np.transpose(pred, (1,2,0,3))[:,:,:,0]  # drop last dimension; (x,y,depth)

        # save binary prediction output
        pred_fname = img_fname[i][0:-7] + model_fname[j][0:-14] + "_pred.tif"
        tiff.imwrite(pred_fname, np.transpose(pred, (2,0,1)).astype(np.uint8) * 255)  # (depth, x, y)