%{
This function extracts vasculature features.
Inputs:
    - C:
        3d array, contains centerline of extracted vessel. vessels contain
        serial labeling. 
        Source file: "[caseNumber]_[CARegion]_finalCline.nii.gz"
    - V:
        3d array, binary segmentation of vessels. 
        Source file: "[caseNumber]_[CARegion]_G_correctedSeg.nii.gz"
    - N:
        3d array, binary segmentation of neuron center
        Source file: "[caseNumber]_[CARegion]_R_seg.nii.gz"
    - area:
        Float, cross sectional imaged area in um^2
        Default to 340.08*340.08
    - res:
        3x1 array [x, y, z], um per voxel width
        Default to [0.332, 0.332, 1.205]

Outputs:
== numerical outputs ==
    out_vvf: vessel volume fraction
    out_vld: vessel length density
    out_neurond: neuron density
    out_nodevd: node volume density
    out_nodeld: node length density
    out_vsl: vessel segment partitioning
== vector outputs ==
    prep_seg_length: length [um] of each vessel segment
    out_tort: tortuosity of each vessel segment
    out_NCdist: each neuron's distance [um] to its nearest vessel centerline
    out_NVdist: each neuron's distance [um] to its nearest vessel
        segmentation (currently skipped because it takes too long)
== image outputs ==
    vessel centerline projection
    neuron centroid projection
    the first two projections overlay on top of each other

Potentially:
    [] neuron distance to metabolic minimum
    [] neuron distribution in relation to metabolic concentration
%}

%% Set computation variables
fpath = "/Users/yunonnebai/NauenLab/2023_data_analysis/dataset/images/";
output_fname = "/Users/yunonnebai/NauenLab/2023_data_analysis/dataset/outputs_numerical_NEWNAME.xlsx"; % sheet names: common, tortuosity, segment_length, neuron_centerline_distance
case_num = ["A1","A2","A3","A4","A5","A6", ...
            "B1","B2","B3","B4","B5","B6", ...
            "C1","C2","C3","C4","C5","C6", ...
            "D1","D2","D3","D4","D5","D6"];
CA_num = ["CA1","CA2","CA3","CA4"];
output_row_num = string(2:97);
coverslip_compression = [15	40	60	0	10	10	0	10	0	0	0	45	0	10	0	0	0	0	0	0	20	60	0	60	0	0	0	0	0	0	0	60	0	0	0	0	20	10	40	30	55	0	0	0	0	0	0	0	0	0	0	0	20	0	0	0	25	0	0	0	30	0	0	0	0	0	0	0	0	0	0	50	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	10	35	0	0	0	0	0	0];
case_num_neuron = ["A4", "B2", "B3", "C1", "C2", "C5"];
coverslip_compression_neuron = [0 10 0 0 0 0 0 60 0 0 0 0 0 0 0 0 20 0 0 0 0 0 0 0];
res = [0.332, 0.332, 1.205];
area = 340.08*304.08;

%% Main: vasculature analysis
counter = 0;
for i=1:length(case_num)
    for j=1:length(CA_num)
        counter = counter+1;
        % file input
        fname = case_num(i) + "_" + CA_num(j)
        C = niftiread(fpath + fname + "_finalCline.nii.gz");
        V = niftiread(fpath + fname + "_G_correctedSeg.nii.gz");
        V = V(:,:,coverslip_compression(counter)+1 : end);
        
        % Ensure that the segmentation matrices do not have negative values
        if (min(C,[],'all') ~= 0)   
            C(C<0) = 0;
            "Input centerline matrix contains negative values; now normalized to 0"
        end
        
        % Metric variable calculation
        prep_Cseg = seg_centerline(C);
        [prep_seg_length, out_tort, out_weighted_tort] = segment_based_calculations(prep_Cseg, size(C), res);
        out_nvv = normalized_vess_vol(V);
        out_nvl = nomalized_vess_leng(prep_seg_length, size(C), area, res);
        out_avd = avg_vess_diameter(V, prep_seg_length, res);
        out_nodevol = node_per_vol(C, size(C), area, res);
        out_nodeleng = node_per_leng(C, prep_seg_length);
        out_vsp = vessel_segment_partitioning(prep_seg_length);
        
        % Save output metric variable results
        save(fname+"_outputs.mat");
        writematrix(out_nvv,output_fname,'Sheet','common','Range',"J"+output_row_num(counter));
        writematrix(out_nvl,output_fname,'Sheet','common','Range',"K"+output_row_num(counter));
        writematrix(out_avd,output_fname,'Sheet','common','Range',"L"+output_row_num(counter));
        writematrix(out_nodevol,output_fname,'Sheet','common','Range',"M"+output_row_num(counter));
        writematrix(out_nodeleng,output_fname,'Sheet','common','Range',"N"+output_row_num(counter));
        writematrix(out_vsp,output_fname,'Sheet','common','Range',"O"+output_row_num(counter));
        writematrix(out_tort,output_fname,'Sheet','tortuosity','Range',"E"+output_row_num(counter));
        writematrix(out_weighted_tort,output_fname,'Sheet','weighted_tortuosity','Range',"E"+output_row_num(counter));
        writematrix(prep_seg_length,output_fname,'Sheet','segment_length','Range',"E"+output_row_num(counter));
    end
end

%% Main: neuron analysis
counter = 0;
for i=1:length(case_num_neuron)
    for j=1:length(CA_num)
        counter = counter+1;
        % file input
        fname = case_num_neuron(i) + "_" + CA_num(j)
        C = niftiread(fpath + fname + "_finalCline.nii.gz");
        N = niftiread(fpath + fname + "_R_segmentation.nii.gz");
        N = N(:,:,coverslip_compression_neuron(counter)+1 : end);
        N = N>0;    % binarize neuron segmentation matrix
        V = niftiread(fpath + fname + "_G_correctedSeg.nii.gz");
        V = V(:,:,coverslip_compression(counter)+1 : end);
        
        % Ensure that the segmentation matrices do not have negative values
        if (min(C,[],'all') ~= 0)   
            C(C<0) = 0;
            "Input centerline matrix contains negative values; now normalized to 0"
        end

        % prepartory variable calculation
        prep_Ncentroid = neuron_centroid(N);

        % output metric variable calculation
        out_neurond = neuron_density(N, area, res);
        out_NVdist = distance_to_vessel(V, N, res);
        out_NNdist = distance_to_neuron(N, res);

        % Save output metric variable results
        save(fname+"_outputs_neuron.mat");
        writematrix(out_neurond,output_fname,'Sheet','common','Range',"P"+output_row_num(counter));
        writematrix(out_NVdist,output_fname,'Sheet','NV_distance','Range',"E"+output_row_num(counter));
        writematrix(out_NNdist,output_fname,'Sheet','NN_distance','Range',"E"+output_row_num(counter));
    end
end


%% Functions
function out = dist_um(pt1, pt2, res)
%{
    Calculate the eudlidean distance between two points in um (corrected 
    for resolution)
    input: 
        Coordinate of the two points [which row, which column, which layer]
        Resolution [x,y,z]
    output: distance between the two points in um
%}
    % get distance in three directions
    x = res(1) * (pt1(1)-pt2(1));
    y = res(2) * (pt1(2)-pt2(2));
    z = res(3) * (pt1(3)-pt2(3));
    % get hypotenus
    out = sqrt(x^2 + y^2 + z^2);
end
function out = path_um(path, res)
%{
    Calculate the length of path (curve) in um (corrected for resolution)
    input:
        - path: coordinates of all points on the path (n by 3 array)
    output: 
        - out: (float) path length in um
%}
    out = 0;    % initialize output variable
    for i=1:(size(path,1)-1)    % loop through all intervals
        % get coordinates of two neighboring points
        temp_pt1 = path(i,:);   
        temp_pt2 = path(i+1,:);
        dist_um(temp_pt1, temp_pt2, res);   % distance between the two points in um
        out = out + dist_um(temp_pt1, temp_pt2, res);
    end
end
function out = seg_centerline(C)
%{ 
    Remove the nodes from centerline and output only disconnected segments
    input: 
        C: vessel centerline
    output: 
        out: segment centerline (struct, index of each segment formatted as 
            [which row, which column, which layer])
%}
    % find nodes
    Cconv = convn(C, ones(3,3,3), 'same');
    Cbifur = (Cconv.*C)>3;
    % remove nodes to extract segments only
    Cseg = C - Cbifur;
    out = bwconncomp(Cseg);
end
function out = neuron_centroid(N)
%{
    Compute the centroid voxel from each segmented neruon "blob"
    input:
        - N: binarized neuron segmentation matrix (logical)
    output
        - out: index of each neuron centroid voxel, formatted as 
               [which row, which column, which layer])
%}
    % get centroid location, formatted as [x,y,z]
    stat = regionprops3(N, "Centroid"); % calculate centroid voxel
    stat = stat.Centroid;    % extract voxel location from table
    stat(all(isnan(stat),2),:) = [];  % remove nan values
    % change output format to [which row, which column, which layer]
    out = [stat(:,2), stat(:,1), stat(:,3)];  
end
function [out_sl, out_tort, out_wtort] = segment_based_calculations(Cseg, C_shape, res)
%{
    Compute the length of each vessel segments in um
    input:
        Cseg: centerline of each segment, with nodes removed
        C_shape: shape of centerline matrix ([x,y,z])
        res: um per voxel width ([x,y,z])
    output: 
        out: vector storing the length of each segment in um ([1,n])
%}
    
    % Initialize variables
    out_sl = zeros(1,size(Cseg.PixelIdxList, 2));
    out_tort = zeros(1, length(out_sl));
    endpts_dist = zeros(1, length(out_sl));

    for i=1:length(out_sl)  % loop through each segment
        % gather all points within this segment
        seg_temp = Cseg.PixelIdxList{i};    % gather all points on that segment
        [x,y,z] = ind2sub(C_shape, seg_temp);   % turn point index into 3d coordinate
        points = [x y z];
        
        % initialize variables to prepare for reordering points
        num_points = size(points, 1);
        ordered_points = zeros(num_points, 3);
        ordered_points(1, :) = points(1, :);
        remaining_points = 1:num_points;
        remaining_points(1) = [];

        % find the closest points and order them
        for j=2:num_points
            last_point = ordered_points(j-1, :);
            distances = sqrt(sum((points(remaining_points, :) - last_point).^2, 2));
            [~, closest_idx] = min(distances);
            ordered_points(j, :) = points(remaining_points(closest_idx), :);
            remaining_points(closest_idx) = [];
        end
        
        % calculate segment length in um
        out_sl(i) = path_um(ordered_points, res);
        % calculate tortuosity
        endpts_dist(i) = dist_um(ordered_points(1,:), ordered_points(end,:), res); 
        out_tort(i) = out_sl(i) / endpts_dist(i);
    end
    
    % calculate weighted tortuosity
    total_path_length = sum(out_sl);
    out_wtort = out_tort .* (out_sl / total_path_length);
end
function out = normalized_vess_vol(V)
%{
    Vessel area fraction =  vessel volume / total image volume
    purpose: spatial density of vessel that incorporates vessel thickness
    input: 
        - V: vessel binary segmentation
    output: 
        - (float) vessel area fraction [unitless]
    Note: no need to account for resolution because it's canceled
%}
    out = sum(V,'all') / (size(V,1) * size(V,2) * size(V,3));
end
function out = nomalized_vess_leng(seg_length, C_shape, area, res)
%{
    vessel length density = centerline length / total image volume
    purpose: spatial vessel density; oxygenation/nutrient delivery dysfunc
    input: 
        seg_length: length of each vessel segment in um
        C_shape: shape of centerline matrix ([x,y,z])
        area: cross sectional area [num of pixels]
        res: um per voxel width ([x,y,z])
    output: (float) vessel length density [1/um^2]
%}

    % compute total segment length
    path_leng = sum(seg_length);
    
    % calculate length density
    out = path_leng / (area*C_shape(3)*res(1)*res(2)*res(3));
            
end
function out = avg_vess_diameter(V, seg_length, res)
%{
    average vessel diameter of all featured vessels
    input: 
        V: vessel binary segmentation
        seg_length: length of each vessel segment in um
        C_shape: shape of centerline matrix ([x,y,z])
        area: cross sectional area [num of pixels]
        res: um per voxel width ([x,y,z])
    output: (float) average vessel diameter [um]
%} 
    vol_um = sum(V,'all') * res(1) * res(2) * res(3);   % total vessel vol in um^3
    length_um = sum(seg_length);    % total vess length in um
    out = sqrt(vol_um / (pi * length_um)) * 2;  % average vessel diameter in um
end
function out = node_per_vol(C, C_shape, area, res)
%{
    node volume density = number of nodes / total image volume
    purpose: interconnectedness of network; resistance to occlusion or
    blockage; flow dynamics
    input: 
        C: vessel centerline
        C_shape: shape of centerline matrix ([x,y,z])
        area: cross sectional area [num of pixels]
        res: um per voxel width ([x,y,z])
    output: (float) node density [count / um]
%}
    % find nodes
    Cconv = convn(C, ones(3,3,3), 'same');
    Cbifur = (Cconv.*C)>3;  % points on cline with more than 3 connected voxels
    CC = bwconncomp(Cbifur);    % identify clusters of points
    num_nodes = size(CC.PixelIdxList, 2);   % each cluster represents one point
    
    % node volume density
    out = num_nodes / (area*C_shape(3)*res(1)*res(2)*res(3));
end
function out = node_per_leng(C, seg_length)
%{
    node length density = number of nodes / centerline length
    purpose: interconnectedness of network; resistance to occlusion or
    blockage; flow dynamics
    input: 
        C: vessel centerline
        seg_length: length of each vessel segment in um
    output: (float) node density [count / um]
%}
    % find nodes
    Cconv = convn(C, ones(3,3,3), 'same');
    Cbifur = (Cconv.*C)>3;  % points on cline with more than 3 connected voxels
    CC = bwconncomp(Cbifur);    % identify clusters of points
    num_nodes = size(CC.PixelIdxList, 2);   % each cluster represents one point
    
    % compute total segment length
    path_leng = sum(seg_length);
    
    % node density
    out = num_nodes / path_leng;
end
function out = vessel_segment_partitioning(seg_length)
%{
    # of segments / segment centerline length in um
    purpose: interconnectedness of network; resistance to occulusion; flow
    dynamics
    input: 
        seg_length: length of each vessel segment in um
    output: (float) vessel segment partitioning (count/um)
%}
    
    % compute total segment length
    path_leng = sum(seg_length);
    
    % calculate partitioning
    num_seg = length(seg_length);
    out = num_seg / path_leng;

end
function out = neuron_density(N, area, res)
%{ 
    neuron density = number of neurons / total image volume
    purpose: energy demand
    input: 
        N: neuron segmentation
        area: image cross sectional area (number of pixels)
        res: um per voxel width ([x,y,z])
    output: neuron density [count / um^3]
%}
    Ncentroid = neuron_centroid(N);
    out = length(Ncentroid) / 3 / (area*res(1)*res(2)*res(3));   % /3 because stat gives x,y,z data
end
function out = distance_to_vessel(V, N, res)
%{
  Euclidean distance between neuron's centroid and its nearest vessel  
  (by vol segmentation) in um.
- A cube-shaped neighborhood is selected with the centroid being the
  origin, extending 40um in radius directions. Analysis is only
  performed within this neighborhood. 
- If vessels are found within the neighborhood, return the euclidean distance. 
- If no vessel is found, return NaN.
    inputs: 
        V: vessel segmentation matrix
        N: neuron segmentation matrix
        res: um per voxel width ([x,y,z])
    output: (1d array) each neuron centroid's distance to shortest vessel centerline [um]
%}
    
    % setting param
    padwidth = 120; % buffer on 3 dir added to C, V, N
    radius_xy = 120;   % ~40um; radius on x/y of analyzed neighborhood from centroid
    radius_z = 33;   % ~40um; radius on z of analyzed neighborhood from centroid

    % padding, so neighborhood slicing doesn't go out of bound
    padsize = [padwidth, padwidth, padwidth];
    V_padded = padarray(V, padsize);
    N_padded = padarray(N, padsize);
    Ncentroid = neuron_centroid(N_padded);

    % slicing neighborhood, so the analysis doesn't take forever
    out = nan(1, length(Ncentroid)); % initialize output matrix
    for i=1:length(Ncentroid)
        % select neighborhood
        nbhood_ctr = floor(Ncentroid(i,:));
        nbhood_V = V_padded(nbhood_ctr(1)-radius_xy:nbhood_ctr(1)+radius_xy, ...
                            nbhood_ctr(2)-radius_xy:nbhood_ctr(2)+radius_xy, ...
                            nbhood_ctr(3)-radius_z:nbhood_ctr(3)+radius_z);
        if max(nbhood_V,[],'all')~=0
            vess = find(nbhood_V);
            dist = zeros(length(vess), 1); % temp array to store distances
            for j = 1:length(vess)
                % coordinate of the jth cline voxel
                [x, y, z] = ind2sub(size(nbhood_V), vess(j));
                % dist between the j-th cline voxel and the i-th centroid
                dist(j) = dist_um([radius_xy+1,radius_xy+1,radius_z+1], [x,y,z], res);
            end
            if min(dist)<=40
                out(i) = min(dist);
            end
        end
    end
end
function out = distance_to_neuron(N, res)
%{
  Euclidean distance between each neuron centroid's and its closest neuron in um.
    inputs: 
        N: neuron segmentation matrix
        res: um per voxel width ([x,y,z])
    output: (1d array) each neuron centroid's distance to nearest neuron centroid [um]
%}
    Ncentroid = neuron_centroid(N);
    out = nan(1, length(Ncentroid)); % initialize output matrix
    for i=1:length(Ncentroid)
        dist = nan(length(Ncentroid), 1); % temp array to store distances
        for j = 1:length(Ncentroid)
            dist(j) = dist_um(Ncentroid(i,:), Ncentroid(j,:), res);
        end
        out(i) = min(dist(dist>0));
    end
end
