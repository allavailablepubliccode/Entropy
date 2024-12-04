clear;clc;close all;
tpoints = 200;

for ii = 1:5
    clear movie
    load(['~/Dropbox/Two_Photon/M' num2str(ii) '/Natural_Movies/movie.mat'])
    movie(movie==0)=nan;
    allmovies{ii} = movie(:,:,1:tpoints);
end

playdata(allmovies,tpoints)

function playdata(allmovies,tpoints)
figure;
set(gcf,'position',[238 64 787 473])
colormap viridis
for t=1:tpoints
    for ii = 1:5
        subplot(2,3,ii)
        p = allmovies{ii}(:,:,t);
        h = imagesc(p);
        set(h, 'AlphaData', ~isnan(p))
        axis off
    end
    subplot(2,3,6)
    title(['t = ' num2str(t),'/',num2str(tpoints)], 'Position', [0.5 0.5], 'FontSize', 10);
    axis off
    drawnow;
end
end