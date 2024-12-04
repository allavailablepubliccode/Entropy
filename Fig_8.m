clear;clc;close all;
load('tbls.mat','NconsS_m','consS_m')
erode_no    = 3;
erode_it    = 2;

figure
colormap viridis
set(gcf,'position',[1 305 1024 232])
c = 1;
for M = 1:5

    C = consS_m(:,M);
    NC = NconsS_m(:,M);

    C = C/max(abs(C));
    NC = NC/max(abs(NC));

    load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/regions.mat'],'map')
    eroded_map = customErode(double(logical(map)),erode_no, erode_it);
    eroded_regions = eroded_map .* map;

    clear ind
    for ii = 1:max(map(:))
        ind{ii} = find(map == ii);
    end

    subplot(2,5,c)
    mapdisp = nan(size(map));
    for jj = 1:numel(C)
        mapdisp(ind{jj}) = C(jj);
    end
    h = imagesc(mapdisp);
    set(h, 'AlphaData', ~isnan(mapdisp))
    axis off
    clim([-1 1])
    colorbar
    drawnow

    c = c + 5;

    subplot(2,5,c)
    mapdisp = nan(size(map));
    for jj = 1:numel(NC)
        mapdisp(ind{jj}) = NC(jj);
    end
    h = imagesc(mapdisp);
    set(h, 'AlphaData', ~isnan(mapdisp))
    axis off
    clim([-1 1])
    colorbar
    drawnow

    c = c - 4;

end

function erodedImage = customErode(image,N,N2)
erodedImage = image;
for ii = 1:N2
    se                     = ones(N);
    [imageRows, imageCols] = size(erodedImage);
    [seRows, seCols]       = size(se);
    padSizeRows            = floor(seRows / 2);
    padSizeCols            = floor(seCols / 2);
    paddedImage            = ones(imageRows + 2 * padSizeRows, imageCols + 2 * padSizeCols);
    paddedImage(padSizeRows + 1:end - padSizeRows, padSizeCols + 1:end - padSizeCols) = erodedImage;
    erodedImage            = ones(imageRows, imageCols);
    for i = 1:imageRows
        for j = 1:imageCols
            neighborhood          = paddedImage(i:i + seRows - 1, j:j + seCols - 1);
            if all(neighborhood(se == 1) == 1)
                erodedImage(i,j) = 1;
            else
                erodedImage(i,j) = 0;
            end
        end
    end
    erodedImage(1,:)   = 0;
    erodedImage(end,:) = 0;
    erodedImage(:,1)   = 0;
    erodedImage(:,end) = 0;
end
end