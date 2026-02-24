path = "/Volumes/LaCie/NauenLab/files_24cases/A4/";
zstack = "A4_CA3";

C = niftiread(path + zstack + "_finalCline.nii.gz");

skelDilateR = 3;
skel = imdilate(C, strel('sphere', skelDilateR));

niftiwrite(skel, path + "_dilated", "Compressed", true);