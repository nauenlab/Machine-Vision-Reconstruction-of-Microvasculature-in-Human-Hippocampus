%% For all files
% file import
path = "/Users/yunonnebai/NauenLab/2022_7/July25PostFixingCenterlining/";
zstack = "D3_CA4";
C = niftiread(path + zstack + "_initialCline_corrected.nii.gz");

% post fixing centerlining
[Skel, C_new] = post_fixing_centerlining(C);

% save
niftiwrite(C_new, path + zstack + "_finalCline", "Compressed", true);


%% For testing/visualization only
% find where nodes are
C_binary = C_new>0;
Cconv = convn(C_binary, ones(3,3,3), 'same');
Cbifur = (Cconv.*C_binary)>3;

% colorcode bifurcation points on the cline map
C_nodemap = C_binary + Cbifur;

% dilate
C_dilated = imdilate(C_new, strel('sphere', 4));
C_nodemap = imdilate(C_nodemap, strel('sphere', 4));

% save
niftiwrite(C_dilated, path + zstack + "_finalCline_dilated", "Compressed", true);
niftiwrite(C_nodemap, path + zstack + "_finalCline_nodemap", "Compressed", true);

% find the number of bifurcations
CC = bwconncomp(Cbifur);
num_nodes = size(CC.PixelIdxList, 2)


%% Functions
function [Skel, C_new] = post_fixing_centerlining(C)
    % parameters
    deadEndLimit = 40; % to remove dead end hairs
    smallVesselLimit = 20; % to remove short end-to-end vessels
    
    % binarize
    C = C>0;
    
    % skeletonize
    skel = Skeleton3D(C);
    
    % noise removal
    change = 1;
    while change
        % remove short dead end vessels
        change = skel;
        C = convn(skel, ones(3, 3, 3), 'same');
        C2 = convn(skel .* (C==2), ones(3, 3, 3),'same') .* skel;   % find ends of vessels
        CC1 = skel .* (C==3);   % find middle segment of vessels (no branches)
        CC0 = bwconncomp(CC1);  % returns all connected elements in CC1
        for i = 1:CC0.NumObjects
            if numel(CC0.PixelIdxList{i}) < deadEndLimit && ...
                    any(C2(CC0.PixelIdxList{i}([1, end])))  % if it is an "end vessel" and shorter than limit
                skel(CC0.PixelIdxList{i}) = 0;  % then remove that vessel
            end
        end
        % remove single pixel connected to node
        C = convn(skel, ones(3, 3, 3),'same');
        C2 = convn(C>3, ones(3, 3, 3), 'same');     % find only branches
        CC0 = (skel>0) .* (C==2) .* C2;     % find pixels that's both an endpoint and next to a node
        skel(CC0>0) = 0;
        % remove isolated pixel
        C = convn(skel, ones(3, 3, 3), 'same');
        skel(and(skel, C==1)) = 0;  % find single isolated pixel
        
        %%%%%% YB's edit: added this section
        % remove isolated double-voxel edges
        C = convn(skel, ones(3, 3, 3), 'same');
        CC1 = skel .* (C==2);   % find all endpoints
        CC0 = bwconncomp(CC1);  % find all connected endpoints (double voxel edges)
        for i=1:CC0.NumObjects
            if numel(CC0.PixelIdxList{i})== 2
                skel(CC0.PixelIdxList{i}) = 0;  % remove double-voxel edges
            end
        end
        %%%%%%
        
        % remove single pixle loop
        C = convn(skel, ones(3, 3, 3), 'same');
        C2 = convn(skel .* (C>3), ones(3, 3, 3), 'same') .* skel;
        CC1 = skel .* (C==3) .* (C2==2);    % voxels that are both an edge and next to a node
        for k = find(CC1)'
            CC0 = bwconncomp(skel);
            skel(k) = 0;    % removing the voxel shouldn't change the number of edges
            CC1 = bwconncomp(skel);
            if CC0.NumObjects ~= CC1.NumObjects
                skel(k) = 1;    
            end
        end
        change = ~isempty(find(change ~= skel, 1));
    end
    
    % generate 
    CC0 = bwconncomp(skel);
    Skel=cell(3, CC0.NumObjects);   % info summary of each vessel
    C_new = zeros(size(skel), 'single');  % to create a 3d array where voxels value of a vessel == its serial number
    counter = 1;
    for i=1:CC0.NumObjects
        if numel(CC0.PixelIdxList{i}) > smallVesselLimit    % remove vessels that are too small
            C_new(CC0.PixelIdxList{i}) = counter;
            [x,y,z] = ind2sub(size(skel), CC0.PixelIdxList{i});
            Skel{1,counter} = [x,y,z];    % location info of all voxels in this vessel
            Skel{2,counter} = 5;  % idk what this is
            Skel{3,counter} = counter;  % serial number of this vessel
            counter = counter + 1;
        end
    end
    keep = any(~cellfun('isempty',Skel),1);
    Skel = Skel(:,keep);

end