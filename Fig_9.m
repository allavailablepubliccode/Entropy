clear;clc;close all;
erode_no    = 3;
erode_it    = 2;

figure
set(gcf,'position',[1 417 1024 120])
c = 1;
for M = 1:5

    disp(['mouse M' num2str(M)])

    load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/regions.mat'],'map')
    load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/movie.mat'],'movie')

    eroded_map = customErode(double(logical(map)),erode_no, erode_it);
    eroded_regions = eroded_map .* map;

    % x = 1:size(movie,2);
    % y = 1:size(movie,1);
    % for ii = 1:size(movie,3)
    %     for jj = 1:max(map(:))
    %         tmp = eroded_regions;
    %         tmp(tmp ~= jj) = 0;
    %         tmp = double(logical(tmp));
    %         tmp = tmp .* movie(:,:,ii);
    %         tmp = tmp / sum(tmp(:));
    %
    %         px          = sum(tmp);
    %         py          = sum(tmp,2);
    %         mu_x(ii,jj) = sum(px .* x);
    %         mu_y(ii,jj) = sum(py .* y');
    %         D_x(ii,jj)  = sum(px .* x.^2);
    %         D_y(ii,jj)  = sum(py .* y'.^2);
    %
    %     end
    % end
    % D_x          = mean(diff(D_x - mu_x.^2));
    % D_y          = mean(diff(D_y - mu_y.^2));
    % mu_x         = mean(diff(mu_x));
    % mu_y         = mean(diff(mu_y));
    load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/coeffs.mat'],'D_x','D_y','mu_x','mu_y')

    clear ind ind2
    for ii = 1:max(map(:))
        ind{ii}  = find(eroded_regions == ii);
        ind2{ii} = find(map == ii);
    end

    % clear dS_dt dconsS_dt dconsS_dt_x dconsS_dt_y
    % for ii = 1:size(movie,3)-1
    %     p          = movie(:,:,ii);
    %     logp       = log2(p);
    %     logp(p==0) = 0;
    %     sig        = - p .* logp;
    %     dsigdx     = [zeros(1,size(sig,1))', diff(sig,1,2)]                        .* eroded_map;
    %     dsigdy     = [zeros(1,size(sig,2))', diff(sig,1,1)']'                      .* eroded_map;
    %     d2sigdx2   = [zeros(1,size(sig,1))', diff(sig,2,2), zeros(1,size(sig,1))'] .* eroded_map;
    %     d2sigdy2   = [zeros(1,size(sig,2)); diff(sig,2,1); zeros(1,size(sig,2))]   .* eroded_map;
    %     dlogpdx    = [zeros(1,size(sig,1))', diff(logp,1,2)]                       .* eroded_map;
    %     dlogpdy    = [zeros(1,size(sig,2)); diff(logp,1,1)]                        .* eroded_map;
    %     for jj = 1:numel(ind)
    %         dconsS_dt(ii,jj)   = mean(D_x(ii,jj)*d2sigdx2(ind{jj}) + D_y(ii,jj)*d2sigdy2(ind{jj}) - mu_x(ii,jj)*dsigdx(ind{jj}) - mu_y(ii,jj)*dsigdy(ind{jj}));
    %     end
    % end
    % dconsS_dt = mean(dconsS_dt);
    load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/dS.mat'],'dconsS_dt')

    subplot(1,5,c)
    mapdisp = nan(size(map));
    for jj = 1:numel(ind2)
        mapdisp(ind2{jj}) = dconsS_dt(jj)*10^7;
    end
    h = imagesc(mapdisp);
    set(h, 'AlphaData', ~isnan(mapdisp))
    axis off
    colorbar
    drawnow

    c = c + 1;

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