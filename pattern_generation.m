% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% This file is part of the code available at
% https://github.com/iivek/sparse-synthextures
% which comes under GPL-3.0 license.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%  TODO: the update rules can be paralelized
%

% Algorithm parameters
iterations = 5;
epsilon = 1e-12;
beta = 2;   % 2 is a must for now
assert(beta==2, 'Has to be 2 - Euclidean distance.')
lambda = 1.0; % sparsity penalization

% Parameters related to the dataset
height = size(T,1);     % number of elements in a patch.
factors = size(T,2);    % number of factors
width = size(S,1);   % number of activation vectors (that is, the number of patches)

% Initialization of the fixed patches, a sort of seeds that are taken from
% the original image
seedPatches = cell(numel(toConsider),1);
randomized = randperm(size(X,2));
for loc = find(~toConsider)'
    randomPatch = X(:,randomized(loc));
%     randomPatch = T*V(:,randomized(loc));
%     randomPatch = zeros(size(randomPatch));
%     randomPatch = randomPatch';
%     randomPatch(1:4) = 1;
%     randomPatch(17:21) = 1;
%     randomPatch = randomPatch';
    seedPatches{loc} = randomPatch(:);  % take a random patch
end

% Initialization. TODO: Robust and automatic initial values
V = rand(factors,width)/20;

% Not letting activations become larger than they are in the
% original decomposition
activationLimits = max(V,[],2);

obj = Inf;
for i=1:iterations    
    %
    % Unit l1 norm on the activations.
    % In this case l1 penalization won't work. But we can penalize ratio of
    % l1 and l2 to promote sparsity (direct connections with Hoyer's sparsity)
    %
    % Traversing patch by patch, while keeping other patches fixed.
    % Objective function
    
    for patch = 1:width
        if( ~toConsider(patch) )
            % This means that we have fixed this patch and it should not
            % get optimized.
            continue
        end
        overlaps = find(~cellfun('isempty',S(patch,:,1)));
        % Traverse patches that 'patch' overlaps with. There the divergence
        % needs to be taken into account and the gradients need to be
        % accumulated (we accumulate in 'positive' and negative')
        negative = 0;
        positive = 0;
        v1 = V(:,patch);
        for overlap = overlaps
            S1 = S{patch,overlap,1};
            S2 = S{patch,overlap,2};
%             % First step - sanity check:
%             % draw the patches that overlap with patch 'patch', in the
%             % global image space            
%             inWholeImage = zeros(imageSize);
%             inWholeImage(templateInfo(2,1,patch):templateInfo(2,1,patch)+patchSize(1)-1,...
%                 templateInfo(2,2,patch):templateInfo(2,2,patch)+patchSize(2)-1) = ...
%                 reshape(S1*ones(size(T))*V(:,patch),patchSize);
%             imagesc(inWholeImage); pause
%             % Seems alright :)
            T1 = S1*T;
            if( toConsider(overlap) )
                y = S2*T*V(:,overlap);
            else
               % The overlapping patch is a fixed one, not a part of the optimization 
                y = S2*seedPatches{overlap};
            end
            negative = negative...
                + T1'*((beta-1).*y.*(T1*v1).^(beta-2)+y.^(beta-1));
            positive = positive...
                + beta.*T1'*(T1*v1);
        end
        % l1/l2
        negative = negative + lambda.*v1.*norm(v1,1)./norm(v1,2).^3;
        positive = positive + lambda./norm(v1,2);
        % l1/l2 ends
%         l1-l2
%         negative = negative + lambda.*v1./norm(v1,2);
%         positive = positive + lambda;
%         % l1-l2 ends
        numerator = negative; % + repmat(sum(v1.*positive,1),[factors,1]);
        denominator = positive; % + repmat(sum(v1.*negative,1), [factors,1]); 
        v1 = v1.*numerator./denominator;
        overflow = v1>activationLimits;
        v1(overflow) = activationLimits(overflow);
%         v1 = v1./repmat(sum(v1,1),[factors,1]);    % unit L1 norm accross cols
        V(:,patch) = v1;
    end
    
    curr = 0;   % Accumulator for objective function
    for patch = 1:width
        if( ~toConsider(patch) )          
            continue
        end
        overlaps = find(~cellfun('isempty',S(patch,:,1)));
        divPart = 0;
        for overlap = overlaps
            S1 = S{patch,overlap,1};
            S2 = S{patch,overlap,2};
            T1 = S1*T;
            if(toConsider(overlap))
                y = S2*T*V(:,overlap);
            else
               % The overlapping patch is a fixed one, not a part of the optimization 
                y = S2*seedPatches{overlap};
            end
            v1 = V(:,patch);
            % symmetric stuff
            divPart = divPart...
                + 1./(beta-1).*sum( (T1*v1-y).*((T1*v1).^(beta-1)-y.^(beta-1)) );       
        end
        divPart = divPart/2;  % dividing by two because for each pair symmetric divergence has been added twice
        curr = curr + divPart + lambda*norm(v1,1)./norm(v1,2); % l1/l2
        % curr = curr + divPart + lambda*(norm(v1,1)-norm(v1,2)); % l1-l2
    end    
    
    
    %
    figure(1)
    obj = [obj curr];
    if obj(end)>obj(end-1)
        'broken'
        obj(end)-obj(end-1)
        %        break
    end
    subplot(3,1,1); imagesc(T*V); title(num2str(i));
    subplot(3,1,2); imagesc(V);
    subplot(3,1,3); plot(obj);
    
    theWhole = zeros([imageSize, channels]);
    for patch = 1:size(templates,2)
        if( ~toConsider(patch) )
             y = seedPatches{patch};
        else
            y = T*V(:,patch);
        end
        theWhole(templateInfo(2,1,patch):templateInfo(2,1,patch)+patchSize(1)-1,...
            templateInfo(2,2,patch):templateInfo(2,2,patch)+patchSize(2)-1,:) = ...
            reshape(y,[patchSize,channels]);
%         y = (y - min(y(:))) ./ ( max(y(:))-min(y(:)) );
%         if patch == 1
%             figure(11); imagesc(reshape(y,[patchSize,channels]));
%         else
%             figure(12);  imagesc(reshape(y,[patchSize,channels]));
%         end
%         drawnow
    end
    theWhole = (theWhole - min(theWhole(:))) ./ ( max(theWhole(:))-min(theWhole(:)) );
    figure(2);
    imagesc(theWhole(:,:,:)); colormap bone
    
    theWhole = zeros([imageSize, channels]);
    for patch = 1:size(templates,2)
        if( ~toConsider(patch) )
            continue
        else
            y = T*V(:,patch);
        end
        theWhole(templateInfo(2,1,patch):templateInfo(2,1,patch)+patchSize(1)-1,...
            templateInfo(2,2,patch):templateInfo(2,2,patch)+patchSize(2)-1,:) = ...
            reshape(y,[patchSize,channels]);
%        
    end
    theWhole = (theWhole - min(theWhole(:))) ./ ( max(theWhole(:))-min(theWhole(:)) );
    figure(3);
    imagesc(theWhole(:,:,:)); colormap bone

    drawnow
end

theWhole = zeros([imageSize, channels]);
for patch = 1:size(templates,2)
    theWhole(templateInfo(2,1,patch):templateInfo(2,1,patch)+patchSize(1)-1,...
        templateInfo(2,2,patch):templateInfo(2,2,patch)+patchSize(2)-1,:) = ...
        max(reshape(T*V(:,patch),[patchSize channels]),...
            theWhole(templateInfo(2,1,patch):templateInfo(2,1,patch)+patchSize(1)-1,...
            templateInfo(2,2,patch):templateInfo(2,2,patch)+patchSize(2)-1,:));
end
theWhole = (theWhole - min(theWhole(:))) ./ ( max(theWhole(:))-min(theWhole(:)) );
   

% Cleaning up the gradients where the patches overlap using median
% filtering
img = theWhole;
[H,W,C] = size(img); img = double(img); 
gx = zeros(H,W,channels); gy = zeros(H,W,channels); j = 1:H-1; k = 1:W-1;
gx(j,k,:) = (img(j,k+1,:) - img(j,k,:)); gy(j,k,:) = (img(j+1,k,:) - img(j,k,:));
gluedAt = [repeatEvery(2):repeatEvery(2):size(gx,2), patchSize(2):repeatEvery(2):size(gx,2)];
for i=1:channels
    filtered(:,:,i) = medfilt2(gx(:,:,i),[3,3]);
    gx(:,gluedAt,i) = filtered(:,gluedAt,i);
end
for i=1:channels
    filtered(:,:,i) = medfilt2(gy(:,:,i),[3,3]);
    gy(gluedAt,:,i) = filtered(gluedAt,:,i);
end

% Reconstruct image from gradients
for i=1:channels
    img_rec(:,:,i) = poisson_solver_function(gx(:,:,i),gy(:,:,i),img(:,:,i));
end
figure;imagesc(img);colormap gray;colorbar;title('Image')
img_rec = (img_rec - min(img_rec(:))) ./ ( max(img_rec(:))-min(img_rec(:)) );
figure;imagesc(img_rec);colormap gray;colorbar;title('Reconstructed');
% figure;imagesc(abs(img_rec-img));colormap gray;colorbar;title('Abs error'); 



%% Below are things for way2 of pattern generation
% %% Let's reconstruct the entire image. Happy happy happy joy joy joy.
% theWhole = zeros(imageSize);
% for patch = 1:theOtherImageBeginsAt-1
%     theWhole(templateInfo(2,1,patch):templateInfo(2,1,patch)+patchSize(1)-1,...
%         templateInfo(2,2,patch):templateInfo(2,2,patch)+patchSize(2)-1) = ...
%         reshape(T*V(:,patch),patchSize);
% end
% 
% theWhole2 = zeros(imageSize);
% for patch = theOtherImageBeginsAt:size(templates,2)
%     theWhole2(templateInfo(2,1,patch):templateInfo(2,1,patch)+patchSize(1)-1,...
%         templateInfo(2,2,patch):templateInfo(2,2,patch)+patchSize(2)-1) = ...
%         reshape(T*V(:,patch),patchSize);
% end
% 
% 
% figure; imagesc(theWhole(patchSize(1)+1:end-ceil(patchSize(1)/2)-1,patchSize(2)+1:end-ceil(patchSize(2)/2))-1); colormap bone
% figure; imagesc(theWhole2(ceil(patchSize(1)/2)+1:end-patchSize(1)-1,ceil(patchSize(2)/2)+1:end-patchSize(2)-1)); colormap bone
% figure; imagesc(theWhole+theWhole2); colormap bone
% 
% figure
% difference = sqrt((theWhole-theWhole2).^2);
% imagesc((difference(patchSize(1)+1:end-patchSize(1)-1,patchSize(2)+1:end-patchSize(2)-1)))


