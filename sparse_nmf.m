%% Ad-hoc script which does nonnegative matrix factorization with sparsity
%  regularization, to find the dictionary (in this case the left matrix,
%  with dictionary elements organised as columns).
%  I am using a variant where the left matrix is l1 (sum-to-one) normalized
%  accross cols, and 1l/l2 sparsity constraint is imposed on the right
%  matrix (the activation matrix).
%  Basically any other sparse (nonnegative) dictionary learning method should work
%
%  Note: kind of converges slowly - initialization?
%  Or use any other method for sparse (nonnegative) dictionary learning


% X is input matrix
extract_patches;
% X = (X - repmat(min(X,[],1),[size(X,1),1])) ./ ( repmat(max(X,[],1)-min(X,[],1),[size(X,1),1]) );
X = (X - min(X(:))) ./ ( max(X(:))-min(X(:)) );


% Algorithm parameters
iterations = 1000;
epsilon = 1e-12;
beta = 2;
lambda = 0.75; % sparsity penalization

% Parameters related to the dataset
[height width] = size(X);
factors = 256;
% Initialization
V = rand(factors,width);
% V = V./repmat(sum(V,1),[factors,1]);    % unit l1 norm accross cols
T = rand(height,factors);
T = T./repmat(sqrt(sum(T.^2,1)),[height,1]);    % unit norm accross cols
    
% T = X/V;
% T(T<0) = epsilon;

R = T*V;
obj = Inf;
for i=1:iterations    

    %% Update T
    piece1 = (R.^(beta-2).*X)*V';
    piece2 = R.^(beta-1)*V';
    numerator = piece1 + T.*repmat(sum(T.*piece2,1),[height,1]);
    denominator = piece2 + T.*repmat(sum(T.*piece1,1), [height,1]);
    denominator = max(denominator,epsilon);
    T = T.*numerator./denominator;
    T = T./repmat(sqrt(sum(T.^2,1)),[height,1]);    % unit norm accross cols
    
    R = T*V;
    
    %% Update V
    piece1 = T'*(R.^(beta-2).*X);
    piece2 = T'*(R.^(beta-1));
    negative = piece1 + lambda.*V.*repmat( sum(V,1)./sqrt(sum(V.^2,1)).^3, [factors,1]);
    positive = piece2 + lambda./repmat(sqrt(sum(V.^2,1)), [factors,1]); 
    numerator = negative; %+ repmat(sum(V.*positive,1),[factors,1]);
    denominator = positive; %+ repmat(sum(V.*negative,1), [factors,1]); 
    V = V.*numerator./denominator;
%     V = V./repmat(sum(V,1),[factors,1]);    % unit L1 norm accross cols
    
    R = T*V;
    
    %% Calculate objective function
    switch beta
        case 0  % IS divergence
            curr = sum(sum( X./R-log(X./R)))-1*height*width;
        case 1  % generalized KL-divergence
            curr = sum(sum( X.*(log(X)-log(R)) + (R-X) ));
        otherwise
            curr = 1./(beta.*(beta-1)).*sum(sum( X.^beta-R.^beta...
                -beta.*R.^(beta-1).*(X-R) ));
    end
    curr = curr + lambda.*sum( sum(V,1)./sqrt(sum(V.^2,1)) );
%     curr = curr + lambda*sum(V(:));%./norm(V,2);
    
    %    
    obj = [obj curr];
    if obj(end)>obj(end-1)
        'Convergence violated.'
        obj(end)-obj(end-1)
        %        break
    end
    figure(1)
    subplot(2,1,1); imagesc(data);
        %% Let's show the reconstruction
        %  Variables from extract_patches need to be intact
        %
        recon = zeros(size(data));
        cnt = 1;
        for row=1:numel(topPatch)
            for col=1:numel(leftPatch)
                % top-left point of the patch
                l = leftPatch(col); 
                t = topPatch(row);
                recon(t:t+patchSize(1)-1,l:l+patchSize(2)-1,:) = reshape(T*V(:,cnt),[patchSize,channels]);
                cnt = cnt+1;
        %         bla=data;
        %         bla(t:t+patchSize(1)-1,l:l+patchSize(2)-1) = 0;
        %         imagesc(bla); pause
            end
        end
        recon = (recon - min(recon(:))) ./ ( max(recon(:))-min(recon(:)) );

    subplot(2,1,2); imagesc(recon); title(num2str(i));
    
    figure(2)
    subplot(3,1,1); imagesc(V); title(num2str(i));
    subplot(3,1,2); hist(V(:),20);
    subplot(3,1,3); plot(obj);
    drawnow
   
end

save checkpoint

%% Let's show the reconstruction
%  Variables from extract_patches need to be intact
%
recon = zeros(size(data));
cnt = 1;
for row=1:numel(topPatch)
    for col=1:numel(leftPatch)
        % top-left point of the patch
        l = leftPatch(col); 
        t = topPatch(row);
        recon(t:t+patchSize(1)-1,l:l+patchSize(2)-1,:) = reshape(T*V(:,cnt),[patchSize,channels]);
        cnt = cnt+1;
%         bla=data;
%         bla(t:t+patchSize(1)-1,l:l+patchSize(2)-1) = 0;
%         imagesc(bla); pause
    end
end
recon = (recon - min(recon(:))) ./ ( max(recon(:))-min(recon(:)) );



%% Let's visualize the patches
percol = round(sqrt(factors));
perrow = ceil(factors/sqrt(factors));

leftPatch = 1:patchSize(2):percol.*patchSize(2);
topPatch = 1:patchSize(1):perrow.*patchSize(1);
numPatches = numel(topPatch)*numel(leftPatch);
suchNice = zeros(perrow.*patchSize(1),perrow.*patchSize(1),channels);

%X = zeros(prod(patchSize),numPatches);
cnt = 1;
for row=1:numel(topPatch)
    for col=1:numel(leftPatch)
        % top-left point of the patch
        l = leftPatch(col); 
        t = topPatch(row);
        suchNice(t:t+patchSize-1,l:l+patchSize(1)-1,:) = reshape(T(:,cnt)./max(T(:,cnt)),[patchSize,channels]);
        cnt = cnt+1;
    end
end

figure; imagesc(recon), colormap bone
figure; imagesc(suchNice); colormap bone
