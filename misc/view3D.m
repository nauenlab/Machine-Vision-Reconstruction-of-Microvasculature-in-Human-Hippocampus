C(C>0) = 1;
V = single(V);
view3D(V, C);

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