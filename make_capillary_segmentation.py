from skimage import io
import numpy as np
import hippoReconstruction as hippo
import nibabel as nib
import pandas as pd

gmodel_fname = "condition2_without_279341_best-model.h5"
meta = pd.read_csv("make_capillary_segmentation_meta.csv").values.tolist()
thresh = 0.47
meta = meta[0:4]

for i in range(len(meta)):
    print("[INFO] now processing: " + meta[i][0])
    print("[INFO] importing zstack green channel images")
    gstack = io.imread(meta[i][0])    # tiff encoding: [depth, x, y]
    gstack = gstack / np.max(gstack)   # normalize input files to within [0,1]

    print("[INFO] trimming unwanted layers; total layers trimmed: ", meta[i][1])
    gstack = gstack[int(meta[i][1]):,:,:]

    print("[INFO] cutting zstack green channel images to tiles")
    gtiles = hippo.image_to_tile(gstack)
    print("[QC] image dataset should be normalized to within [0,1]. Maximum value now: ", np.max(gtiles))

    print("[INFO] network making binary segmentation prediction on capillaries")
    gpred = hippo.unet2d_prediction(gtiles, gmodel_fname, gstack.shape, threshold=thresh)
    gpred = np.transpose(gpred, (2,1,0,3))[:,:,:,0]  # [depth,x,y,1] --> [y,x,depth]

    print("[INFO] removing green bubbles")
    gpred = hippo.green_bubble_removal(gpred, min_size=7000, assign_diff_labels=False)

    nib.save(nib.Nifti1Image(gpred, np.eye(4)), meta[i][0][9:-8] + "processedSeg.nii.gz")
