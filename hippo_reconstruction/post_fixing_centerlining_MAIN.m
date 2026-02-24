%% For all files
zstack = ["A1_CA1", "A1_CA2", "A1_CA3", "A1_CA4", ...
          "A2_CA1", "A2_CA2", "A2_CA3", "A2_CA4", ...
          "A3_CA1", "A3_CA2", "A3_CA3", "A3_CA4", ...
          "A4_CA1", "A4_CA2", "A4_CA3", "A4_CA4", ...
          "A5_CA1", "A5_CA2", "A5_CA3", "A5_CA4", ...
          "A6_CA1", "A6_CA2", "A6_CA3", "A6_CA4", ...
          "B1_CA1", "B1_CA2", "B1_CA3", "B1_CA4", ...
          "B2_CA1", "B2_CA2", "B2_CA3", "B2_CA4", ...
          "B3_CA1", "B3_CA2", "B3_CA3", "B3_CA4", ...
          "B4_CA1", "B4_CA2", "B4_CA3", "B4_CA4", ...
          "B5_CA1", "B5_CA2", "B5_CA3", "B5_CA4", ...
          "B6_CA1", "B6_CA2", "B6_CA3", "B6_CA4", ...
          "C1_CA1", "C1_CA2", "C1_CA3", "C1_CA4", ...
          "C2_CA1", "C2_CA2", "C2_CA3", "C2_CA4", ...
          "C3_CA1", "C3_CA2", "C3_CA3", "C3_CA4", ...
          "C4_CA1", "C4_CA2", "C4_CA3", "C4_CA4", ...
          "C5_CA1", "C5_CA2", "C5_CA3", "C5_CA4", ...
          "C6_CA1", "C6_CA2", "C6_CA3", "C6_CA4", ...
          "D1_CA1", "D1_CA2", "D1_CA3", "D1_CA4", ...
          "D2_CA1", "D2_CA2", "D2_CA3", "D2_CA4", ...
          "D3_CA1", "D3_CA2", "D3_CA3", "D3_CA4", ...
          "D4_CA1", "D4_CA2", "D4_CA3", "D4_CA4", ...
          "D5_CA1", "D5_CA2", "D5_CA3", "D5_CA4", ...
          "D6_CA1", "D6_CA2", "D6_CA3", "D6_CA4"];

for i=1:length(zstack)
    % file import
    zstack_curr = zstack(i);
    C = niftiread(zstack_curr + "_initialCline_corrected.nii.gz");
    "Now processing: " + zstack_curr
    % post fixing centerlining
    [Skel, C_new] = post_fixing_centerlining(C);
    % save
    niftiwrite(C_new, "output_" + zstack_curr + "_finalCline", "Compressed", true);
    save("output_" + zstack_curr + "_finalSkel.mat", 'Skel')

end




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