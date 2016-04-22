% %% Dataset 1
% %
% % Extracting image patches from the images
% filename = './datasets/260058.jpg'
% data = imread(filename);
% 
% % Defining rectangle of interest.
% top = 200;
% left = 1;
% bottom = 321;
% right = 257;
% patchSize = [16, 16];
% step = [2,2];   % stepsize between patches
% 
% leftPatch = left:step(1):right-patchSize(2);
% topPatch = top:step(2):bottom-patchSize(1);
% numPatches = numel(topPatch)*numel(leftPatch);
% 
% X = zeros(prod(patchSize),numPatches);
% cnt = 1;
% for row=1:numel(topPatch)
%     for col=1:numel(leftPatch)
%         % top-left point of the patch
%         l = leftPatch(col); 
%         t = topPatch(row);
%         patch = data(t:t+patchSize-1,l:l+patchSize(1)-1);
%         X(:,cnt) = patch(:);
%         cnt = cnt+1;
% %         % Visualize patches one by one
% %         bla=data;
% %         bla(t:t+patchSize-1,l:l+patchSize(1)-1) = 0;
% %         imagesc(bla); pause
% %         %/Visualize patches one by one
%     end
% end
%
% %/ Dataset 1


%% Dataset 2
%
% Extracting image patches from the images
filename = './datasets/9221235455_02570f8348_o.jpg'
data = imread(filename);
% data = rgb2gray(data);
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
%         % Visualize patches one by one
%         bla=data;
%         bla(t:t+patchSize-1,l:l+patchSize(1)-1) = 0;
%         imagesc(bla); pause
%         %/Visualize patches one by one
    end
end
%
%/ Dataset 2