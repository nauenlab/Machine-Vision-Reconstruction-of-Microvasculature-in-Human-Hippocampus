# packages for I/O
import os
import nibabel as nib
import tifffile as tiff
import matplotlib.pyplot as plt
import pandas as pd
# packages for matrix operations
import numpy as np
import math
# packages for ML
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_curve, auc
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Input, BatchNormalization, Dropout, Conv2D, \
    Conv2DTranspose, MaxPooling2D, concatenate, add
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
from tensorflow_addons.losses import SigmoidFocalCrossEntropy
from tensorflow.keras.losses import Reduction
from tensorflow.keras.optimizers.schedules import ExponentialDecay
plt.style.use("ggplot")



###============ Experiment setup ============###
#init_lr = 3e-04    # initial learning rate
epoch_num = 100     # max number of training epochs
output_prefix = "condition2_without_"
batch_size = 32
dataset_exp_condition = "dataset_fname_condition2_without.csv"
data_layers = 71   # hardcoded; number of 1024x1024 layers in the input dataset



###============ Import data ============###
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
    curr_im = tiff.imread(fname[i][0] + '.tif')  # tiff encoding: (depth, x, y); uint16 encoding
    curr_im = curr_im / (np.max(curr_im))   # normalize into within [0,1]
    # write into pre-allocated matrices
    label[counter:counter+curr_label.shape[0], :, :] = curr_label
    im[counter:counter+curr_im.shape[0], :, :] = curr_im
    # update counter
    counter = counter + curr_label.shape[0]

# Reshape dataset for ML
im = np.reshape(im, im.shape + (1,))
label = np.reshape(label, label.shape + (1,))

# Set parameters for images
num_slides = im.shape[0]
raw_width = im.shape[1]
raw_height = im.shape[2]

# Set parameters for subtiles
im_width = 128
im_height = 128
border = 5
dataset_size = int(num_slides * math.floor(raw_width/im_width) * math.floor(raw_height/im_height))

# Initialize the lists to store images and labels
X = np.zeros((dataset_size, im_height, im_width, 1), dtype=np.float32)
Y = np.zeros((dataset_size, im_height, im_width, 1), dtype=np.float32)

n = 0   # starting location of storing image/label into X and Y
for k in range(num_slides):
    for i in range(math.floor(raw_width/im_width)):
        for j in range(math.floor(raw_height/im_height)):
            X[n, :, :, :] = im[k, i*im_width:(i+1)*im_width, j*im_height:(j+1)*im_height]
            Y[n, :, :, :] = label[k, i*im_width:(i+1)*im_width, j*im_height:(j+1)*im_height]
            n = n+1

# Split train:validation:test = 70:20:10
print("[QC] image dataset should be normalized to within [0,1]. Maximum value now: ", np.max(X))
X_train_val, X_test, Y_train_val, Y_test = train_test_split(X, Y, test_size=0.1, random_state=42, shuffle=True)
X_train, X_val, Y_train, Y_val = train_test_split(X_train_val, Y_train_val, test_size=0.22, random_state=42, shuffle=True)
print("All images loaded")
print("Training dataset size: ", len(X_train))
print("Validation dataset size: ", len(X_val))
print("Held-out testing dataset size: ", len(X_test))



###============ ML architecture ============###
# Define one convolution block (== two convolution layers)
def conv2d_block(input_tensor, n_filters, kernel_size = 3, batchnorm = True):
    
    # first layer
    x = Conv2D(filters = n_filters, kernel_size = (kernel_size, kernel_size), activation='relu', kernel_initializer='he_normal', padding='same')(input_tensor)
    if batchnorm:
        x = BatchNormalization()(x)
    
    # second layer
    x = Conv2D(filters = n_filters, kernel_size = (kernel_size, kernel_size), activation='relu', kernel_initializer='he_normal', padding='same')(x)
    if batchnorm:
        x = BatchNormalization()(x)
    return x

# Define overall UNet architecture
def get_unet(X, n_filters=16, dropout=0.1, batchnorm=True):

    # Retrieve the size of inputs
    X_shape = X.shape
    inputs = tf.keras.Input(shape=(X_shape[1],X_shape[2], X_shape[3]))

    # contracting path
    conv1 = conv2d_block(inputs, n_filters * 1, kernel_size = 3, batchnorm = batchnorm)
    pool1 = MaxPooling2D((2,2))(conv1)
    dropout1 = Dropout(dropout)(pool1)

    conv2 = conv2d_block(dropout1, n_filters * 2, kernel_size = 3, batchnorm = batchnorm)
    pool2 = MaxPooling2D((2,2))(conv2)
    dropout2 = Dropout(dropout)(pool2)

    conv3 = conv2d_block(dropout2, n_filters * 4, kernel_size = 3, batchnorm = batchnorm)
    pool3 = MaxPooling2D((2,2))(conv3)
    dropout3 = Dropout(dropout)(pool3)

    conv4 = conv2d_block(dropout3, n_filters * 8, kernel_size = 3, batchnorm = batchnorm)
    pool4 = MaxPooling2D((2,2))(conv4)
    dropout4 = Dropout(dropout)(pool4)

    # Horizontal path
    conv5 = conv2d_block(dropout4, n_filters * 16, kernel_size = 3, batchnorm = batchnorm)

    # Expansive path
    convtrans6 = Conv2DTranspose(n_filters * 8, (3, 3), strides = (2, 2), padding = 'same')(conv5)
    conc6 = concatenate([convtrans6, conv4])
    dropout6 = Dropout(dropout)(conc6)
    conv6 = conv2d_block(dropout6, n_filters * 8, kernel_size = 3, batchnorm = batchnorm)

    convtrans7 = Conv2DTranspose(n_filters * 4, (3, 3), strides = (2, 2), padding = 'same')(conv6)
    conc7 = concatenate([convtrans7, conv3])
    dropout7 = Dropout(dropout)(conc7)
    conv7 = conv2d_block(dropout7, n_filters * 4, kernel_size = 3, batchnorm = batchnorm)

    convtrans8 = Conv2DTranspose(n_filters * 2, (3, 3), strides = (2, 2), padding = 'same')(conv7)
    conc8 = concatenate([convtrans8, conv2])
    dropout8 = Dropout(dropout)(conc8)
    conv8 = conv2d_block(dropout8, n_filters * 2, kernel_size = 3, batchnorm = batchnorm)

    convtrans9 = Conv2DTranspose(n_filters * 1, (3, 3), strides = (2, 2), padding = 'same')(conv8)
    conc9 = concatenate([convtrans9, conv1])
    dropout9 = Dropout(dropout)(conc9)
    conv9 = conv2d_block(dropout9, n_filters * 1, kernel_size = 3, batchnorm = batchnorm)

    # Final block
    outputs = Conv2D(1, (1, 1), activation='sigmoid')(conv9)

    # Organize the model
    model = Model(inputs=inputs, outputs=outputs)
    
    return model


# Fix random seed for reproducibility
# The Answer to the Ultimate Question of Life, the Universe, and Everything is 42,
# quote from The Hitchhiker's Guide to the Galaxy.
seed = 42
np.random.seed(seed)

# Clear any previous model
tf.keras.backend.clear_session()

# Initialize the model
input_img = Input((im_height, im_width, 1), name='img')
model = get_unet(input_img, n_filters=16, dropout=0.05, batchnorm=True)
model.compile(
    optimizer='Adam',
    loss=SigmoidFocalCrossEntropy(alpha=0.5, gamma=2, reduction=Reduction.SUM_OVER_BATCH_SIZE),
    metrics=[tf.keras.metrics.BinaryAccuracy(),
             tf.keras.metrics.TruePositives(),
             tf.keras.metrics.TrueNegatives(),
             tf.keras.metrics.FalsePositives(),
             tf.keras.metrics.FalseNegatives()]
)

callbacks = [
    EarlyStopping(patience=15, verbose=1),
    ReduceLROnPlateau(factor=0.1, patience=5, min_lr=0.00001, verbose=1),
    ModelCheckpoint(output_prefix+"best-model.h5", verbose=1, save_best_only=True)
]



###============ ML training ============###
# Apply the network model on input dataset
results = model.fit(X_train, Y_train, batch_size=batch_size, epochs=epoch_num, callbacks=callbacks, validation_data=(X_val, Y_val))

# Plot out the learning epoch
plt.figure(figsize=(8, 8))
plt.title("Learning curve")
plt.plot(results.history["loss"], label="loss")
plt.plot(results.history["val_loss"], label="val_loss")
plt.plot(np.argmin(results.history["val_loss"]), np.min(results.history["val_loss"]), marker="x", color="r",
         label="best model")
plt.xlabel("Epochs")
plt.ylabel("log_loss")
plt.legend()
plt.savefig(output_prefix + "learning_epoch.pdf")



###============ ML performance evaluation ============###
# Evaluate on validation set (this must be equals to the best log_loss)
best_model = load_model(output_prefix + "best-model.h5")
print("Now evaluating the best model on validation set: \n")
best_model.evaluate(X_val, Y_val, verbose=1)

# Evaluation on held-out testing set
print("Now evaluating the best model on held-out testing set: \n")
# metrics
best_model.evaluate(X_test, Y_test, verbose=1)
# PR curve and AUC
Y_pred_nothresh = best_model.predict(X_test, verbose=1)
precision, recall, thresh = precision_recall_curve(Y_test.flatten(), Y_pred_nothresh.flatten())
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
plt.savefig(output_prefix + 'PR_curve.png')
