output_fname = "/Users/yunonnebai/NauenLab/2023_data_analysis/dataset/outputs_numerical_NEWNAME.xlsx"; % sheet names: common, tortuosity, segment_length, neuron_centerline_distance
nvdist = readmatrix(output_fname, sheet=5);

for i=1:size(nvdist,1)
    for j=1:size(nvdist,2)
        if nvdist(i,j) > 40
            nvdist(i,j) = NaN;
        end
    end
end