clear;clc;close all;
load('tbls.mat', 'consS_x_m', 'consS_y_m', 'consS_x_s', 'consS_y_s')
erode_no    = 3;
erode_it    = 2;
scalefact2 = 100;
errorfactor = 10;
scalefacterror = 10;

figure
set(gcf,'position',[1 1 1024 536])
colormap gray
for M = 1:5
    subplot(2,3,M)

    load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/regions.mat'],'map')

    eroded_regions2 = double(logical(map));
    eroded_regions2(eroded_regions2==0) = nan;

    scalefact = max(sqrt(consS_x_m(:,M).^2 + consS_y_m(:,M).^2));

    imagesc2(eroded_regions2)
    hold on
    for jj = 1:6

        mx = consS_x_m(jj,M)/scalefact;
        my = consS_y_m(jj,M)/scalefact;

        sx = consS_x_s(jj,M)/scalefact;
        sy = consS_y_s(jj,M)/scalefact;

        mx = mx*scalefact2;
        my = my*scalefact2;
        sx = sx*scalefact2*scalefacterror;
        sy = sy*scalefact2*scalefacterror;

        [row,col] = find(map == jj);
        n = numel(row);

        row = round(mean(row));
        col = round(mean(col));

        clear center
        center(1) = col;
        center(2) = row;

        quiver(center(1), center(2), mx, my, 0, 'r', 'LineWidth', 1, 'MaxHeadSize', 0.1);
        rectangle('Position', [center(1) + mx - sx, center(2) + my - sy, 2*sx, 2*sy], 'Curvature', [1, 1], 'EdgeColor', 'g', 'LineStyle', '-','LineWidth',1);

    end
    hold off
    axis off
    drawnow

end


function h = imagesc2 ( img_data )

h = imagesc(img_data);
axis image off

if ndims( img_data ) == 2
    set(h, 'AlphaData', ~isnan(img_data))
elseif ndims( img_data ) == 3
    set(h, 'AlphaData', ~isnan(img_data(:, :, 1)))
end

if nargout < 1
    clear h
end
end

function factor = CIcalc(CI,n)
alpha = (100-CI)/100;
if n > 30
    critical_value = norminv(1 - alpha/2);
else
    critical_value = tinv(1 - alpha/2, n - 1);
end
factor = critical_value;
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