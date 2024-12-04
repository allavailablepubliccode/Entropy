clear; clc; close all;

% folder where data is stored, downloaded from:
% https://figshare.com/articles/dataset/Data_from_Functional_Parcellation_of_Mouse_Visual_Cortex_Using_Statistical_Techniques_Reveals_Response-Dependent_Clustering_of_Cortical_Processing_Areas_/13476522/1
basedir = '~/Dropbox/Two_Photon/';
erode_no    = 3;
erode_it    = 2;

cortical_regions = {'AL' 'AM' 'LM' 'PM' 'RL' 'V1'};

for mouse = 1:5

    disp(['mouse ' num2str(mouse) '/5, re-structuring images'])
    clear movie map
    for jj = 1:numel(cortical_regions)

        movie_tmp = load([basedir 'M' num2str(mouse) '/Natural_Movies/data.mat'],[cortical_regions{jj}]);
        movie_tmp = movie_tmp.([cortical_regions{jj}]);

        map_tmp = load([basedir 'M' num2str(mouse) '/Natural_Movies/map.mat'],[cortical_regions{jj} '_MAP']);
        map_tmp = map_tmp.([cortical_regions{jj} '_MAP']);

        for kk = 1:size(map_tmp,1)
            movie(map_tmp(kk,1), map_tmp(kk,2), :) = movie_tmp(kk,:);
            map(map_tmp(kk,1), map_tmp(kk,2)) = jj;
        end
    end

    disp(['mouse ' num2str(mouse) '/5, masking images'])
    for jj = 1:size(movie,3)
        movie(:,:,jj) = movie(:,:,jj) .* logical(map);
    end

    disp(['mouse ' num2str(mouse) '/5, cropping images'])
    summed_movie = sum(movie,3);
    summed_movie_1 = sum(summed_movie,1);
    summed_movie_1 = find(summed_movie_1 == 0);
    summed_movie_2 = sum(summed_movie,2);
    summed_movie_2 = find(summed_movie_2 == 0);

    movie(:,summed_movie_1,:) = [];
    movie(summed_movie_2,:,:) = [];

    map(:,summed_movie_1) = [];
    map(summed_movie_2,:) = [];

    disp(['mouse ' num2str(mouse) '/5, normalizing images'])
    for jj = 1:size(movie,3)
        tmp = movie(:,:,jj);
        tmp(tmp == 0) = nan;
        tmp = tmp + abs(min(tmp(:)));
        tmp(isnan(tmp)) = 0;
        tmp = tmp / sum(tmp(:));
        movie(:,:,jj) = tmp;
    end

    disp(['mouse ' num2str(mouse) '/5, saving images'])
    save([basedir 'M' num2str(mouse) '/Natural_Movies/movie.mat'],'movie')
    save([basedir 'M' num2str(mouse) '/Natural_Movies/regions.mat'],'map')

end