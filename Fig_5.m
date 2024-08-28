clear;clc;close all;
load('~/Dropbox/Two_Photon/M3/Natural_Movies/movie_short.mat')
% movie_short is time points 765:2:783 from movie
movie(movie==0) = nan;
figure
set(gcf, 'Position',[1 277 1024 260]);
for ii = 1:size(movie,3)
    subplot(2, 5, ii)
    p = movie(:,:,ii)*10^5;
    h = imagesc(p);
    set(h, 'AlphaData', ~isnan(p))
    axis off
    colorbar
    drawnow
end