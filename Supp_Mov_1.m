clear;clc;close all;
load('~/Dropbox/Two_Photon/M3/Natural_Movies/movie.mat')
playdata(movie(:,:,764:784),1)

function playdata(data, waitt)
figure;
set(gcf,'position',[531 135 494 402])
colormap viridis
data(data==0) = nan;
for t=1:size(data,3)
    p = data(:,:,t);
    h = imagesc(p);
    set(h, 'AlphaData', ~isnan(p))
    title([num2str(t),'/',num2str(size(data,3))]);
    axis off
    drawnow;
    if t==1 || t==size(data,3)
        pause(waitt*2)
    else
        pause(waitt)
    end
end
end