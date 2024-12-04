clear;clc;close all;
erode_no    = 3;
erode_it    = 2;

load('tbls.mat','NconsS_m', 'NconsS_s', 'consS_m', 'consS_s')

m = [NconsS_m'; nan(1,6); consS_m'];
s = [NconsS_s'; nan(1,6); consS_s'];

m = m * 10^8;
s = s * 10^8;

c = 0;
for ii = 7:11
    c = c + 1;
    Cd = distinguishable_means(m(ii,:), s(ii,:));
    kk(c) = sum(Cd(:));
end

[row,col] = find(isnan(m));
row = unique(row);

resultall = strings(23,6);

for ii = 1:size(m,1)
    for jj = 1:size(m,2)
        if isnan(m(ii,jj))
            resultall(ii,jj) = sprintf('', m(ii,jj), s(ii,jj));
        else
            resultall(ii,jj) = sprintf('%.2f ± %.2f', m(ii,jj), s(ii,jj));
        end
    end
end

% for M = 1:5
% 
%     disp(['mouse M' num2str(M)])
% 
%     load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/regions.mat'],'map')
%     load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/movie.mat'],'movie')
% 
%     eroded_map = customErode(double(logical(map)),erode_no, erode_it);
%     eroded_regions = eroded_map .* map;
% 
%     load(['~/Dropbox/Two_Photon/M' num2str(M) '/Natural_Movies/coeffs.mat'],'D_x','D_y','mu_x','mu_y')
% 
%     clear ind ind2
%     for ii = 1:max(map(:))
%         ind{ii}  = find(eroded_regions == ii);
%         ind2{ii} = find(map == ii);
%     end
% 
%     for jj = 1:numel(ind)
% 
%         clear dNconsS_dt dconsS_dt dconsS_x_dt dconsS_y_dt
%         for ii = 1:size(movie,3)-1
% 
%             p          = movie(:,:,ii);
%             logp       = log2(p);
%             logp(p==0) = 0;
%             sig        = - p .* logp;
%             dsigdx     = [zeros(1,size(sig,1))', diff(sig,1,2)]                        .* eroded_map;
%             dsigdy     = [zeros(1,size(sig,2))', diff(sig,1,1)']'                      .* eroded_map;
%             d2sigdx2   = [zeros(1,size(sig,1))', diff(sig,2,2), zeros(1,size(sig,1))'] .* eroded_map;
%             d2sigdy2   = [zeros(1,size(sig,2)); diff(sig,2,1); zeros(1,size(sig,2))]   .* eroded_map;
%             dlogpdx    = [zeros(1,size(sig,1))', diff(logp,1,2)]                       .* eroded_map;
%             dlogpdy    = [zeros(1,size(sig,2)); diff(logp,1,1)]                        .* eroded_map;
% 
%             dNconsS_dt(:,ii)  = D_x(ii,jj)*p(ind{jj}).*dlogpdx(ind{jj}).^2 + D_y(ii,jj)*p(ind{jj}).*dlogpdy(ind{jj}).^2;
%             dconsS_dt(:,ii)   = D_x(ii,jj)*d2sigdx2(ind{jj}) + D_y(ii,jj)*d2sigdy2(ind{jj}) - mu_x(ii,jj)*dsigdx(ind{jj}) - mu_y(ii,jj)*dsigdy(ind{jj});
%             dconsS_x_dt(:,ii) = D_x(ii,jj)*d2sigdx2(ind{jj}) - mu_x(ii,jj)*dsigdx(ind{jj});
%             dconsS_y_dt(:,ii) = D_y(ii,jj)*d2sigdy2(ind{jj}) - mu_y(ii,jj)*dsigdy(ind{jj});
% 
%         end
% 
%         dNconsS_dt  = cumsum(dNconsS_dt,2);
%         dconsS_dt   = cumsum(dconsS_dt,2);
%         dconsS_x_dt = cumsum(dconsS_x_dt,2);
%         dconsS_y_dt = cumsum(dconsS_y_dt,2);
% 
%         dNconsS_dt  = smooth_signal_2D(dNconsS_dt);
%         dconsS_dt   = smooth_signal_2D(dconsS_dt);
%         dconsS_x_dt = smooth_signal_2D(dconsS_x_dt);
%         dconsS_y_dt = smooth_signal_2D(dconsS_y_dt);
% 
%         dNconsS_dt  = diff(dNconsS_dt')';
%         dconsS_dt   = diff(dconsS_dt')';
%         dconsS_x_dt = diff(dconsS_x_dt')';
%         dconsS_y_dt = diff(dconsS_y_dt')';
% 
%         dNconsS_dt  = dNconsS_dt(:);
%         dconsS_dt   = dconsS_dt(:);
%         dconsS_x_dt = dconsS_x_dt(:);
%         dconsS_y_dt = dconsS_y_dt(:);
% 
%         NconsS_m(jj,M) = mean(dNconsS_dt);
%         NconsS_s(jj,M) = std(dNconsS_dt)/sqrt(numel(dNconsS_dt));
% 
%         consS_m(jj,M) = mean(dconsS_dt);
%         consS_s(jj,M) = std(dconsS_dt)/sqrt(numel(dconsS_dt));
% 
%         consS_x_m(jj,M) = mean(dconsS_x_dt);
%         consS_x_s(jj,M) = std(dconsS_x_dt)/sqrt(numel(dconsS_x_dt));
% 
%         consS_y_m(jj,M) = mean(dconsS_y_dt);
%         consS_y_s(jj,M) = std(dconsS_y_dt)/sqrt(numel(dconsS_y_dt));
% 
%     end
% end
% 
% save('tbls.mat', 'NconsS_m', 'NconsS_s', 'consS_m', 'consS_s','consS_x_m', 'consS_x_s', 'consS_y_m', 'consS_y_s')

function smoothed_data = smooth_signal_2D(data)
    % Ensure data is 2D
    [num_signals, num_timepoints] = size(data);
    
    % Parameters
    lowpass_cutoff = 0.01;   % Low-pass filter cutoff frequency (normalized)
    moving_avg_window = 20;  % Window size for moving average
    
    % Step 1: Temporal Smoothing (Second Dimension)
    smoothed_data = zeros(size(data));
    for i = 1:num_signals
        signal = data(i, :);
        [b, a] = butter(2, lowpass_cutoff, 'low');  % 2nd-order Butterworth filter
        smoothed_signal = filtfilt(b, a, signal);
        smoothed_data(i, :) = smoothdata(smoothed_signal, 'movmean', moving_avg_window);
    end
    
    % Step 2: Spatial Smoothing (First Dimension)
    for t = 1:num_timepoints
        smoothed_data(:, t) = smoothdata(smoothed_data(:, t), 'movmean', moving_avg_window);
    end
end

function distinguishableMatrix = distinguishable_means(means, standardErrors)

% Ensure means and standardErrors are column vectors
means = means(:);
standardErrors = standardErrors(:);

% Number of elements
n = length(means);

% Initialize the output matrix
distinguishableMatrix = zeros(n, n);

% Calculate critical difference for each pair
for i = 1:n
    for j = 1:n
        % Calculate the absolute difference of the means
        mean_diff = abs(means(i) - means(j));

        % Calculate the combined standard error for the pair
        combined_SE = sqrt(standardErrors(i)^2 + standardErrors(j)^2);

        % Compare against the critical value for 95% confidence (1.96 * combined_SE)
        if mean_diff > 1.96 * combined_SE
            distinguishableMatrix(i, j) = 1; % Distinguishable
        else
            distinguishableMatrix(i, j) = 0; % Not distinguishable
        end
    end
end
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


