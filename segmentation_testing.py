import numpy as np
import os
import nibabel as nib
import matplotlib.pyplot as plt
import tifffile as tiff
import hippoReconstruction as hippo
from sklearn.metrics import precision_recall_curve, auc
from tensorflow.keras.models import load_model
import pandas as pd

output_prefix = "condition2_with_"
dataset_exp_condition = "dataset_fname_condition2_with.csv"
data_layers = 43   # hardcoded; number of 1024x1024 layers in the input dataset
testing_set_size = 455
model_fname = "condition2_without_279341_best-model.h5"

# import data
print("[INFO] loading input dataset")
fname = pd.read_csv(dataset_exp_condition).values.tolist()
label = np.zeros((data_layers,1024,1024))
im = np.zeros((data_layers,1024,1024))
counter = 0
for i in range(len(fname)):
    # get current labels
    curr_nifti = nib.load(fname[i][0] + '_seg.nii.gz')
    curr_label = curr_nifti.get_fdata() # nii.gz encoding: (y, x, depth)
    curr_label = np.transpose(curr_label, (2,1,0))  # resize into: (depth, x, y)
    curr_label = curr_label > 0.5   # ensure that it is truly binary segmentation
    # get current images
    curr_im = tiff.imread(fname[i][0] + '.tif')  # tiff encoding: (depth, x, y)
    curr_im = curr_im / np.max(curr_im)   # normalize into within [0,1]
    # write into pre-allocated matrices
    label[counter:counter+curr_label.shape[0], :, :] = curr_label
    im[counter:counter+curr_im.shape[0], :, :] = curr_im
    # update counter
    counter = counter + curr_label.shape[0]

# assemble data into testing set
print("[INFO] cutting zstack images and ground truth labels to tiles")
im_tiles = hippo.image_to_tile(im)
label_tiles = hippo.image_to_tile(label)
index = np.arange(len(im_tiles))
np.random.shuffle(index)
index = index[0:testing_set_size]    # randomly select data pairs from images with green bubbles
testing_img = im_tiles[index]
testing_gt = label_tiles[index]
print("[QC] image dataset should be normalized to within [0,1]. Maximum value now: ", np.max(testing_img))
print("Testing dataset size: ", len(testing_img))

# Making predictions on testing images
print("Now evaluating the best model on held-out testing set: \n")
best_model = load_model(model_fname)
# metrics
best_model.evaluate(testing_img, testing_gt, verbose=1)
# PR curve and AUC
testing_pred = best_model.predict(testing_img, verbose=1)
precision, recall, thresh = precision_recall_curve(testing_gt.flatten(), testing_pred.flatten())
auc = auc(recall, precision)
print("UNet AUC: ", auc)
optimal_idx = np.argmax(precision + recall)
optimal_threshold = thresh[optimal_idx]
opti_precision = precision[optimal_idx]
opti_recall = recall[optimal_idx]
print("Optimal threshold: ", optimal_threshold)

plt.figure()
plt.plot(recall, precision, marker='.', label='UNet')
plt.plot(opti_recall, opti_precision, 'go', label="Operating Point")
plt.xlabel('Recall (sensitivity): TP/(TP+FN)')
plt.ylabel('Precision: TP/(TP+FP)')
plt.savefig("without_model_tested_on_with_image_PR_curve.png")