% Copyright 2017-2018, Mohammad Haft-Javaherian. (mh973@cornell.edu)
%   References:
%   -----------
%   [1] Haft-Javaherian, M; Fang, L.; Muse, V.; Schaffer, C.B.; Nishimura, 
%       N.; & Sabuncu, M. R. (2018) Deep convolutional neural networks for 
%       segmenting 3D in vivo multiphoton images of vasculature in 
%       Alzheimer disease mouse models. *arXiv preprint, arXiv*:1801.00880.

clear all;
sample_all = ["D6_CA2", "D6_CA3", "D6_CA4"];
zstart_all = [1, 1, 1];
path = "C:\Users\Saturn\Desktop\yunong bai\centerlining\";

for i=1:length(sample_all)
    
    sample = sample_all(i);
    zstart = zstart_all(i);
    
    % Grab the segmentation
    "Grabing the segmentation"
    input_fname = path + sample + "_G_correctedSeg.nii.gz";
    V_all = niftiread(input_fname);
    V_all = V_all(:,:,zstart:end);
    
    % Post process
    "Starting centerlining post processing"
    [Skel, C, V] = DeepVess_function(V_all);
    
    % Save centerline result
    "Saving centerlining result"
    output_fname = path + sample + "_cline";
    niftiwrite(C, output_fname, "Compressed", true);
    save(path + sample + "_cline.mat", 'C', 'V', "Skel");
end


%% DeepVess function (under optimization)
function [Skel, C, V] = DeepVess_function(V_in)
    %% Default parameters
    pad_size = 10; % padding to make room for dilation
    skelDilateR = 5; % to dilate skeleton for smoothness 
    vDilateR = 1; % to dilate segmentation to improve connectivity
    boxFiltW = 3; % to smooth the segmentation and skeleton results
    
    %%%%%% YB's edit: changed deadEndLimit; added endToEndLimit
%    deadEndLimit = 11; % to remove dead end hairs
    deadEndLimit = 40; % to remove dead end hairs
    smallVesselLimit = 20; % to remove short end-to-end vessels
    %%%%%%
    
    diameterLimit = 25; % to remove large vessels
    elongationLimit = 1; % elongationLimit = Diameter / Length

    %% Crude vessel centerline
    "Extracting crude vessel centerline"
    % Normalize
    V = V_in / max(max(max(V_in)));
    % Averaging to smoothen the boundary
    V = imboxfilt3(V,3) > 0.5; 
    % Fill holes
    V = single(imfill(imboxfilt3(single(V),3)>0.5,'hole'));
%    V = imfill(imboxfilt3(single(V),3)>0.5,'hole');
    % Add a rim of zeros around the xy plane
    V1 = padarray(V, [pad_size, pad_size, pad_size]);
    % Extract skeleton
    skel = Skeleton3D(V1);

    %% Refined vessel centerline
    "refining vessel centerline"
    % make tubes of similar radius + smoothen the tubes
    skel = imdilate(skel, strel('sphere', skelDilateR));
    skel = imboxfilt3(single(skel), boxFiltW)>0.5;
    % fill holes
    for i = 1:size(skel, 3)
        skel(:, :, i) = imfill(skel(:, :, i), 'holes');
    end
    % smoothen + fill holes again
    skel = single(imfill(imboxfilt3(single(skel), boxFiltW)>0.5, 'hole'));
    % skeletonize from the tubes
    skel = Skeleton3D(skel);
    % apply the tube skeleton to the padded vessels --> so that skeleton
    % doesn't run in areas outside of vessel volume
    skel = skel .* imdilate(V1, strel('sphere', vDilateR));     

    %% Noise removal
    "removing noise"
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
        
        %%%%%% YB's edit: disabled because it removes short, non-end edges
%         % remove double pixle loop
%         C = convn(skel, ones(3, 3, 3), 'same');
%         C2 = convn(skel .* (C==2), ones(3, 3, 3), 'same') .* skel;  % convolve around endpoints
%         CC1 = skel .* (C==3);   % find edges only
%         CC0 = bwconncomp(CC1);
%         for i=1:CC0.NumObjects
%             if numel(CC0.PixelIdxList{i})<3 
%                 skel(CC0.PixelIdxList{i}) = 0;  % remove edges less than 3 voxels long
%             end
%         end    
        %%%%%%
        
        change = ~isempty(find(change ~= skel, 1));
    end

    % remove pads: if any pixel in the pad is nonzero, then the pixel next to 
    % the pad is nonzero
    skel(pad_size+1, :, :) =  any(skel(1:pad_size+1, :, :), 1);
    skel(end-pad_size, :, :) =  any(skel((end-pad_size):end, :, :), 1);
    skel(:, pad_size+1, :) =  any(skel(:, 1:pad_size+1, :), 2);
    skel(:, end-pad_size, :) =  any(skel(:, (end-pad_size):end, :), 2);
    skel(:, :, pad_size+1) =  any(skel(:, :, 1:pad_size+1), 3);
    skel(:, :, end-pad_size) =  any(skel(:, :, (end-pad_size):end), 3);

    skel = skel((pad_size+1):(end-pad_size), (pad_size+1):(end-pad_size), ...
        (pad_size+1):(end-pad_size));
    V1 = V1((pad_size+1):(end-pad_size), (pad_size+1):(end-pad_size), ...
        (pad_size+1):(end-pad_size));

    %% remove rounded object and big vessels
    "removing rounded objects and big vessels"
    C = convn(skel, ones(3,3,3), 'same');
    
    %%%%%% YB's edit: incorporate endpoints to final results
%    CC1 = (skel>0) .* (C==3);   % find the middle of vessles (barring end points)
    CC1 = (skel>0) .* ((C==3)+(C==2));  % find edges and endpoints
    %%%%%%

    %%%%%% YB's edit: added this part to store the branch points separately
    skelBifur = (skel>0) .* (C>3);  % find all bifurcation points
    %%%%%%

    CC0 = bwconncomp(CC1);  % find all connected elements
    D = zeros(CC0.NumObjects,1);
    L = zeros(CC0.NumObjects,1);
    bwDist=bwdist(~V1); % Euclidean distance transform (distance from the nearest foreground voxel)
    skelLarge = zeros(size(skel));  % array to store large vessels to be removed
    for i=1:CC0.NumObjects
        a = bwDist(CC0.PixelIdxList{i});
        a = 2 * max(a(a>0)) * (sum(a>0)>(0.5*numel(a)));    % max diameter if sum of radii is longer than half length
        if isempty(a)
            a=nan;
        end
        D(i)= a;
        L(i) = numel(CC0.PixelIdxList{i});
        if L(i)<(elongationLimit*D(i)) || D(i)> diameterLimit ...
                || isnan(D(i)) || D(i)==0
            CC1(CC0.PixelIdxList{i}) = 0;   % remove large/round/no vessel associated
            skelLarge(CC0.PixelIdxList{i}) = 1;   % store voxel indices of removed large edges
        end
    end
    
    %%%%%% YB's edit: add back the bifurcation points
    Creg = convn(CC1, ones(3,3,3), 'same');             % convolve regular vessels
    CC0 = bwconncomp(skelBifur); % group bifurcations together
    CC2 = (Creg>0) .* (skelBifur>0); % vertices adjacent to regular vessels
    skelBifurFinal = zeros(size(skelBifur));
    for i=1:CC0.NumObjects
        if any(CC2(CC0.PixelIdxList{i})>0)    % if any vertex voxel is adjacent to regular vessels
            skelBifurFinal(CC0.PixelIdxList{i}) = 1; % then add the entire vertice group
        end
    end
    CC1 = CC1 + skelBifurFinal; % add "real" vertices back to centerline
    %%%%%%

    %%%%%% YB's edit: disabled D and L (they are not used)
    CC0 = bwconncomp(CC1);
%     D = zeros(CC0.NumObjects,1);
%     L = zeros(CC0.NumObjects,1);
    for i=1:CC0.NumObjects
%         a = bwDist(CC0.PixelIdxList{i});
%         D(i)=2 * median(a(a>0));
%         L(i) = numel(CC0.PixelIdxList{i});
        a=CC0.PixelIdxList{i};
        CC0.PixelIdxList{i}=a(V1(a)>0);   % remove voxels that are not part of orig vessel
    end
    %%%%%%
    

    %% Generate the skeleton outputs
   "generating skeleton outputs"
    Skel=cell(3, CC0.NumObjects);   % info summary of each vessel
    %%%%%% YB's edit: solve "label1 mystery"; remove small vessels
%    C = skel;   % to create a 3d array where voxels value of a vessel == its serial number + 1
    C = zeros(size(skel), 'single');  % to create a 3d array where voxels value of a vessel == its serial number
    counter = 1;
    for i=1:CC0.NumObjects
%        C(CC0.PixelIdxList{i}) = i+1;
        if numel(CC0.PixelIdxList{i}) > smallVesselLimit    % remove vessels that are too small
            C(CC0.PixelIdxList{i}) = counter;
            [x,y,z] = ind2sub(size(skel), CC0.PixelIdxList{i});
            Skel{1,counter} = [x,y,z];    % location info of all voxels in this vessel
            Skel{2,counter} = 5;  % idk what this is
            Skel{3,counter} = counter;  % serial number of this vessel
            counter = counter + 1;
        end
    end
    % remove empty elements in the cell
    keep = any(~cellfun('isempty',Skel),1);
    Skel = Skel(:,keep);
    %%%%%%%

    %%%%%%% YB's edit: disabled because this maneuver will cause vessels to disappear
%     % Fix the path of centerlines in Skel{1,:} to have straight centerlines,
%     %   other wise they will be zig-zagged
%     for i = 1:size(Skel, 2)
%         temp = Skel{1, i};
%         if size(temp, 1) < 4
%             continue
%         end
%         [Di, I] = pdist2(temp, temp, 'euclidean', 'Smallest', 3);
%         A = zeros(size(temp, 1));
%         for j = 1:size(temp, 1)
%             for k = 2:3
%                 if Di(k, j) <= sqrt(3)
%                     A(j,I(k, j)) = 1;
%                 end
%             end
%         end
%         A = single((A + A') > 0);
%         BGobj = biograph(A);
%         dist = allshortestpaths(BGobj, 'Directed', false);
%         dist(isinf(dist)) = 0;
%         [~, j] = max(dist(:));
%         [S, T] = ind2sub(size(A), j(1));
%         [~, path, ~] = shortestpath(BGobj, S, T);
%         if ~isempty(path)
%             Skel{1, i} = temp(path, :);
%         end
%     end

end


%% Functions (YB added this section)
function V = load_tiff(fname)
    info = imfinfo(fname);
    V = zeros([size(imread(fname)) length(info)]);
    numPage = length(info);
    for k=1:numPage
        V(:,:,k) = imread(fname, k);
    end
end

function out = toSubBlocks(tiff, block_height, block_width, block_depth)
    num_blocks = floor(size(tiff,1)/block_height) * floor(size(tiff,2)/block_width) * floor(size(tiff,3)/block_depth);
    out = zeros(num_blocks, block_height, block_width, block_depth);
    id = 1;
    for i=1:floor(size(tiff,1)/block_height)
        for j=1:floor(size(tiff,2)/block_width)
            for k=1:floor(size(tiff,3)/block_depth)
                out(id,:,:,:) = tiff(block_height*(i-1)+1:block_height*i, block_width*(j-1)+1:block_width*j, block_depth*(k-1)+1:block_depth*k);
                id = id+1;
            end
        end
    end
end

function C = skel2C(skel, block_height, block_width, block_depth)
    C = zeros(block_height, block_width, block_depth);
    for j=1:size(skel, 2)
        for k=1:size(skel{1,j}, 1)
            C(skel{1,j}(k,1), skel{1,j}(k,2), skel{1,j}(k,3))=1;
        end
    end
end

function view3D(sub_vol, sub_cline)
    obj = sub_vol*0.1 + sub_cline*0.9;
    intensity = [0 20 255];
    alpha = [0 0.001 1];
    color = ([0 0 0; 0 255 0; 255 0 0]) ./ 255;
    queryPoints = linspace(min(intensity),max(intensity),256);
    alphamap = interp1(intensity,alpha,queryPoints)';
    colormap = interp1(intensity,color,queryPoints);
    volshow(squeeze(obj), 'Colormap', colormap, 'Alphamap', alphamap);
end