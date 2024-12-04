clear;clc;close all;
load('movie_short')
smooth_no   = 3;
bar_space   = 3;
colors{1} = 'k';
colors{2} = 'r';

c = 0;
for ii = [1 size(movie,3)]
    c = c + 1;
    p           = movie(:,:,ii);
    p           = sum(p);
    p           = smoothfn(p,smooth_no);
    p           = p(1:bar_space:end);
    p           = p - min(p);
    min_px(c)   = min(p);
    max_px(c)   = max(p);
end
min_px           = min(min_px);
max_px           = max(max_px);

figure
set(gcf, 'Position',[297 265 728 272]);
subplot(1,2,1)
c = 0;
hold on
for ii = [1 size(movie,3)]
    c = c + 1;
    p = movie(:,:,ii);
    p = sum(p);
    p = smoothfn(p,smooth_no);
    p = p(1:bar_space:end);
    p = p - min(p);
    bar(p, 'FaceColor', 'none', 'EdgeColor', colors{c})
    alpha 0.1
end
ylim([min_px max_px])
axis tight
hold off
drawnow

c = 0;
clear min_px max_px
for ii = [1 size(movie,3)]
    c = c + 1;
    p           = movie(:,:,ii);
    p           = sum(p,2);
    p           = smoothfn(p,smooth_no);
    p           = p(1:bar_space:end);
    p           = p - min(p);
    min_px(c)   = min(p);
    max_px(c)   = max(p);
end
min_px           = min(min_px);
max_px           = max(max_px);

subplot(1,2,2)
c = 0;
hold on
for ii = [1 size(movie,3)]
    c = c + 1;
    p = movie(:,:,ii);
    p = sum(p,2);
    p = smoothfn(p,smooth_no);
    p = p(1:bar_space:end);
    p = p - min(p);
    bar(p, 'FaceColor', 'none', 'EdgeColor', colors{c})
    alpha 0.1
end
camroll(90)
set(gca, 'XDir','reverse')
ylim([min_px max_px])
axis tight
hold off
drawnow

function Sdata = smoothfn(origdata,n)
Sdata = origdata;
for ii = 1:n
    Sdata = smoothdata(Sdata);
end
end