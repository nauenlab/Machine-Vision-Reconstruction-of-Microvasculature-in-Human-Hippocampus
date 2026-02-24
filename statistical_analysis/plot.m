clear all
in_fpath = "/Users/yunonnebai/NauenLab/2023_data_analysis/dataset/images/";
out_fpath = "/Users/yunonnebai/NauenLab/2023_data_analysis/q_20230822/png_numVessVox_v3/";
case_num = ["A1","A2","A3","A4","A5","A6", ...
            "B1","B2","B3","B4","B5","B6", ...
            "C1","C2","C3","C4","C5","C6", ...
            "D1","D2","D3","D4","D5","D6"];
CA_num = ["CA1","CA2","CA3","CA4"];
remove_layer = [15	40	60	0	10	10	0	10	0	0	0	45	0	10	0	0	0	0	0	0	20	60	0	60	0	0	0	0	0	0	0	60	0	0	0	0	20	10	40	30	55	0	0	0	0	0	0	0	0	0	0	0	20	0	0	0	25	0	0	0	30	0	0	0	0	0	0	0	0	0	0	50	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	10	35	0	0	0	0	0	0];
res = [0.332, 0.332, 1.205];

%% Supp figure 4: maximal projection
counter = 0;
for i=1:length(case_num)
    for j=1:length(CA_num)
        counter = counter+1;
        case_num(i) + CA_num(j)

        % image prep
        V = niftiread(in_fpath + case_num(i) + "_" + CA_num(j) + "_G_correctedSeg.nii.gz");
%        V = V(:,:,remove_layer(counter)+1:end);
        cline = niftiread(in_fpath + case_num(i) + "_" + CA_num(j) + "_finalCline.nii.gz");
        V_scaled = V;
        vess_vol = zeros(1,size(V,3));
        for k=1:size(V,3)
            V_scaled(:,:,k) = V(:,:,k) * k;
            vess_vol(end-k+1) = sum(V(:,:,k),'all');
        end
        cline(cline>0) = 1;

        % max proj plot
        diff_img = uint16(max(cline,[],3)) .* (max(V_scaled,[],'all') - max(V_scaled,[],3) + 1);
        MIPimg = imagesc(max(V_scaled,[],3) + diff_img);
        cmap = jet(256);
        cmap(1,:) = 0;
        cmap(256,:) = 1;
        colormap(cmap);
        cb = colorbar('XTick', linspace(1,size(V,3),8), ...
                      'XTickLabel', floor(size(V,3)*1.205 - linspace(1,size(V,3)*1.205,8)));
        cb.Label.String = "Depth [um]";
        cb.Label.Rotation = 270;
        cb.Label.VerticalAlignment = "bottom";
        cb.Label.FontSize = 15;
        axis square
        
        title(case_num(i) + " " + CA_num(j), 'FontSize', 25)
        xticks(linspace(0,1024,10))
        xticklabels(num2cell(floor(linspace(0,1024*0.332,10))));
        yticks(linspace(0,1024,10))
        yticklabels(num2cell(floor(linspace(0,1024*0.332,10))));
        xlabel("x location [um]")
        ylabel("y location [um]")
        set(gcf, 'Position', [100 100 700 600])
        ax = gca();
        ax.FontSize = 15;
        saveas(gcf, out_fpath+case_num(i)+"_"+CA_num(j)+"_maxProj.png");
%        print(out_fpath+case_num(i)+"_"+CA_num(j)+"_maxProj", '-dpdf', '-bestfit')
    end
end


%% Supp figure 3: vessel voxel by depth
counter = 0;
for i=1:length(case_num)
%for i=1
    figure(1)
    ymax = zeros(4,1);
    for j=1:length(CA_num)
        counter = counter+1;
        case_num(i)+CA_num(j)
        V = niftiread(in_fpath + case_num(i) + "_" + CA_num(j) + "_G_correctedSeg.nii.gz");
        V = V(:,:,remove_layer(counter)+1:end);
        vess_vol = zeros(1,size(V,3));
        for k=1:size(V,3)
            vess_vol(end-k+1) = sum(V(:,:,k),'all') * res(1)*res(2)*res(3);
        end
        hold on
        plot(vess_vol, 'LineWidth', 1.5)
        ymax(j) = max(vess_vol);
    end
    xlim([0,249]);
    xticks(linspace(0,249,10));
    xticklabels(floor(linspace(0,249*1.205,10)));
    yticks(linspace(0,max(ymax),10))
    yticklabels(num2cell(floor(linspace(0,max(ymax),10))));
    xlabel("Depth [um]") 
    ylabel("Vessel volume at a given depth [um^3]")
    colororder([1 0 0; 1 0.65 0; 0 1 1; 0 0 1])
    legend({'CA1', 'CA2', 'CA3', 'CA4'})
    axis square
    ax = gca();
    ax.FontSize = 15;
    set(gcf, 'Position', [100 100 700 600])
    saveas(gcf, out_fpath+case_num(i)+"_numVessVox.png");
    clf()
end

