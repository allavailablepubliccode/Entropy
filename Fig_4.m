clear;clc;close all;
thicken = 3;
figure
colormap gray
set(gcf,'position',[1 385 1024 152])
for mm = 1:5
    subplot(1,5,mm)
    load(['~/Dropbox/Two_Photon/M' num2str(mm) '/Natural_Movies/regions.mat'])
    map = padarray(map, [thicken+2 thicken+2], 0, 'both');
    map = bwmorph(map,'remove');
    map = imdilate(map, strel('square', thicken));
    map = double(map);
    map(map == 0) = nan;
    h = imagesc(map);
    set(h, 'AlphaData', ~isnan(map));
    axis off
    drawnow
end

