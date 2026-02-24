# Machine-Vision-Reconstruction-of-Microvasculature-in-Human-Hippocampus

========
Step 1: image preprocessing
- System: FIJI
- Steps:
	1. Open the .czi file in FIJI
	2. “Split channels” to split the green and red images into two separate images
	3. Save metadata
	4. Convert red and green into 8-bit images
	5. Process on green: “image subtraction” between raw green channel image and processed red channel —> “brightness/contrast” manual adjustment —> “subtract background” (50 pixels) —> “despeckle” —> save file as "xx_CAx_G_prep.tif"


========
Step 2: capillary binary segmentation
- System: Python (JHPCE)
- Codes:
	-- make_capillary_segmentation.py: the main execution file
	-- hippoReconstruction.py: all functions to be called on
	-- condition2_without_279341_best-model.h5: model and weights
	-- make_capillary_segmentation_meta.csv: file with images sizes and etc
	-- visualization.py: see some sample images
- Files:
	-- "xx_CAx_G_prep.tif"
- Steps:
	1. Load and normalize preprocessed G images
	2. Trim unwanted layers
	3. Cut images into tiles
	4. CNN makes binary segmentation
	5. Remove green bubbles
	6. Stitch and save file as "xx_CAx_G_processedSeg.nii.gz"


========
Step 3: manually correct capillary binary segmentation
- System: ITK-SNAP (workstation)
- Files:
	-- "xx_CAx_G_prep.tif"
	-- "xx_CAx_G_processedSeg.nii.gz"
- Steps: 
	1. Load the processed G channel and the capillary binary segmentation in ITK-SNAP
	2. Manually correct any incorrect object instance segmentation by painting over existing labels
	3. Save file as "xx_CAx_G_correctedSeg.nii.gz" --> note that this file does NOT remove the layers affected by coverslip compression!! These layers are removed in step 6 (vasculature feature extraction)


========
Step 4: capillary centerline extraction
- System: matlab (workstation)
- Codes:
	-- DeepVessPostProcessingOptimize_14June2022.m
	-- Associated source codes in "private" folder
- Files:
	-- "xx_CAx_G_correctedSeg.nii.gz"
- Steps: 
	1. Run DeepVessPostProcessingOptimize_14June2022.m
	2. Save Skel, C, and V as "xx_CAx_initialCline.mat" and "xx_CAx_initialCline.nii.gz"


======== 
Step 5: correct capillary centerline connection points
- System: ITK-SNAP
- Codes:
	-- cline_dilation.m: dilate centerline for visualization purposes
	-- post_fixing_centerline_MAIN.m
	-- Associated source codes in "private" folder
- Files:
	-- "xx_CAx_initialCline.nii.gz"
- Steps:
	1. Run cline_dilation.m
	2. Load the processed G channel and the dilated centerline in ITK-SNAP
	3. Manually correct any disconnected points on the centerline map, save file as "xx_CAx_initialCline_corrected.nii.gz"
	4. Run post_fixing_centerline_MAIN.m
	5. Save file as "xx_CAx_finalCline.nii.gz" and "xx_CAx_finalSkel.mat"


========
Step 6: vasculature feature extraction
- System: matlab
- Codes:
	-- vasculature_feature_extraction_19Dec2023.m
- Files:
	-- "xx_CAx_finalCline.nii.gz"
	-- "xx_CAx_finalSkel.mat"
	-- "xx_CAx_G_correctedSeg.nii.gz" --> note that this file does NOT have layers affected by coverslip compression removed!! These layers are removed in this step, during file import
- Steps:
	1. Run vasculature_feature_extraction_19Dec2023.m
	2. Save numerical output in the overall output sheet "outputs_numerical.xlsx"
	3. Save all output in the workspace separately as "xx_CAx_outputs.mat"











