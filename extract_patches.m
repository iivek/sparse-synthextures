% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% This file is part of the code available at
% https://github.com/iivek/sparse-synthextures
% which comes under GPL-3.0 license.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear selection
% Extracting image patches from the images
for selection = 1:6
    switch selection
        case(1)
            % Monochrome desert + pyramids
            filename = './datasets/260058.jpg'
            data = imread(filename);
            % Defining rectangle of interest.
            top = 200;
            left = 1;
            bottom = 321;
            right = 257;
            patchSize = [16, 16];
            step = [2,2];   % stepsize between patches
        case(2)
            % Cherries
            filename = './datasets/9221235455_02570f8348_o.jpg'
            data = imread(filename);
            data = imresize(data,0.04);
            imagesc(data); colormap bone
            % Defining rectangle of interest.
            top = 1;
            left = 1;
            bottom = 110;
            right = 180;
            patchSize = [16, 16];
            channels = size(data,3);
            step = [2,2];   % stepsize between patches
        case(3)
            % monochrome fingerprint
            filename = './datasets/fingerprint.png'
            data = imread(filename);
            % data = rgb2gray(data);
            data = imresize(data,0.4);
            imagesc(data); colormap bone
            % Defining rectangle of interest.
            top = 1;
            left = 1;
            bottom = 140;
            right = 140;
            patchSize = [16, 16];
            channels = size(data,3);
            step = [1,1];   % stepsize between patches
        case(4)
            % RGB desert
            filename = './datasets/christmas-desert-scene-for-e-mail.jpg'
            data = imread(filename);
            data = imresize(data,0.1);
            % Defining ROI
            top = 150;
            left = 1;
            bottom = 200;
            right = 75;
            patchSize = [16, 16];
            channels = size(data,3);
            step = [1,1];   % stepsize between patches
        case(5)
            filename = './datasets/leaf-nature-green-spring-122429.jpeg'
            data = imread(filename);
            data = imresize(data,0.1);
            % Defining ROI
            top = 40;
            left = 300;
            bottom = 140;
            right = 400;
            patchSize = [16, 16];
            channels = size(data,3);
            step = [1,1];   % stepsize between patches
        case(6)
            filename = './datasets/Serach_Mystic_Place_Forest_camouflage_tent_2010.jpg'
            data = imread(filename);
            data = imresize(data,0.2);
            % Defining ROI
            top = 250;
            left = 350;
            bottom = 350;
            right = 450;
            patchSize = [16, 16];
            channels = size(data,3);
            step = [1,1];   % stepsize between patches
    end


    channels = size(data,3);
    leftPatch = left:step(1):right-patchSize(2);
    topPatch = top:step(2):bottom-patchSize(1);
    numPatches = numel(topPatch)*numel(leftPatch);

    X = zeros(prod(patchSize)*channels,numPatches);
    cnt = 1;
    for row=1:numel(topPatch)
        for col=1:numel(leftPatch)
            % top-left point of the patch
            l = leftPatch(col); 
            t = topPatch(row);
            patch = data(t:t+patchSize-1,l:l+patchSize(1)-1,:);
            X(:,cnt) = patch(:);
            cnt = cnt+1;
%             bla=data;
%             bla(t:t+patchSize-1,l:l+patchSize(1)-1,:) = 0;
%             imagesc(bla); pause
        end
    end
    
    selection_ = selection; % TODO
    clear selection
    save(strcat('patches', num2str(selection_)))
    selection = selection_;

end


% draw roi over image
out = data;
imagesc(out(top:bottom,left:right,:))