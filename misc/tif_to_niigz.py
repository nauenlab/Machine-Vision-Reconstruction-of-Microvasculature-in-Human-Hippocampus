import numpy as np
import tifffile as tiff
import nibabel as nib

# files
img_fname = ["good_B5_CA1_condition2_with_pred.tif",
             "good_B5_CA1_condition2_without_pred.tif",
             "good_B5_CA1_condition3_all_pred.tif",
             "greenBubbles_C6_CA2_condition2_with_pred.tif",
             "greenBubbles_C6_CA2_condition2_without_pred.tif",
             "greenBubbles_C6_CA2_condition3_all_pred.tif",
             "shadowyNeurons_B5_CA3_condition2_with_pred.tif",
             "shadowyNeurons_B5_CA3_condition2_without_pred.tif",
             "shadowyNeurons_B5_CA3_condition3_all_pred.tif",
             "specklyBkg_C1_CA1_condition2_with_pred.tif",
             "specklyBkg_C1_CA1_condition2_without_pred.tif",
             "specklyBkg_C1_CA1_condition3_all_pred.tif"]

for i in range(len(img_fname)):
    img = tiff.imread(img_fname[i])
    img = img/255
    nib.save(nib.Nifti1Image(np.transpose(img, (2, 1, 0)), np.eye(4)), img_fname[i][0:-3] + "nii.gz")