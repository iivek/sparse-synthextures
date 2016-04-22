% % So here goes. From an overcomplete dictionary derived from image patches,
% % I want to generate new patches, but an entire image with minimal blocky
% % artifacts.
% 
% % Convention. Entire image is divided into overlapping blocks. Indexing
% % starts from top left to bottom left, then towards the right.
% patchDim = 5;
% % Patches are square.
% % Patches are vectorized.
% % Here is an example of how indexing goes for 'patchdim'=5. 
% %
% % x1  x6  x11 x16 x21
% % x2  x7  x12 x17 x22
% % x3  x8  x13 x18 x23
% % x4  x9  x14 x19 x24
% % x5  x10 x15 x20 x25
% %
% % As an illustration, consider a problem where we have two overlapping
% % patches. 
% % Objective: min[ sum_i sum_j D(S_ij_1*T*v_i, S_ij_2*T*v_j) ]. Note that V are
%   variables. The objective comes down to symetric KL divergence, (in case
%   of KL divergence)
% 
% % We'll be seeing the forumlation S*T*V. T will be 
% 
% T = [1, 0, 0;
%      1, 0, 0;
%      0, 1, 0;
%      0, 0, 1;
%      0, 0, 1;]
% v1 = [0 1 3]';
% v2 = [3 0 1]';
% % Suppose the two vectors overlap where the 3s are. This means that the
% % following hold
% S1 = [  0 0 0 0 0;
%         0 0 0 0 0;
%         0 0 0 0 0;
%         0 0 0 1 0;
%         0 0 0 0 1]
% S2 = [  0 0 0 0 0;
%         0 0 0 0 0;
%         0 0 0 0 0;
%         1 0 0 0 0;
%         0 1 0 0 0]
% S1*T*v1 == S2*T*v2

% %% What follows is okay and works but i'm commenting it
% % Now let's generate mapping suitable for image patches. For the convention
% % of how patches are indexed and how the patches are vectorized see above.
% patchSize = [16 16];    % height and width
% repeatEvery = [8,8];    % Determines location of the patches which will compose the entire image. 
% imageSize = [128, 128];
% assert(all(mod(imageSize,patchSize)==0), 'Please choose imageSize(i) as a multiple of PatchSize(i).');
% %
% % For every patch, determining where the other patches overlap with it.
% % Deliberately using for loops to migrate the code more easily to C
% rowIdx = 1:repeatEvery(2):imageSize(2)-patchSize(2)+1;
% colIdx = 1:repeatEvery(1):imageSize(1)-patchSize(1)+1;
% templates = sparse(imageSize(1)*imageSize(2),numel(rowIdx)*numel(colIdx));
% templateInfo = zeros(2,2,numel(rowIdx)*numel(colIdx));
% %
% globalCounter = 1;
% colCounter = 1;
% for col=colIdx
%     rowCounter = 1;
%     for row=rowIdx
%         template = sparse(imageSize(1),imageSize(2));
%         template(row:row+patchSize(2)-1, col:col+patchSize(1)-1) = 1;
%         templates(:,globalCounter) = template(:);
%         templateInfo(1,:,globalCounter) = [rowCounter,colCounter];
%         templateInfo(2,:,globalCounter) = [row,col];
%         templateInfo(:,:,globalCounter);
%         imagesc(template); pause
%         rowCounter = rowCounter + 1;
%         globalCounter = globalCounter + 1;        
%     end
%     colCounter = colCounter + 1;
% end

% % Another specific way to create overlaps - it's best to visualize the templates to
% % understand what's happening here
% % Now let's generate mapping suitable for image patches. For the convention
% % of how patches are indexed and how the patches are vectorized see above.
% patchSize = [16 16];    % height and width
% repeatEvery = [16,16];  % Determines location of the patches which will compose the entire image. 
% offset = [8,8];
% imageSize = [64, 64];
% assert(all(mod(imageSize,patchSize)==0), 'Please choose imageSize(i) as a multiple of PatchSize(i).');
% imageSize = imageSize + offset;
% 
% %% way2
% %
% % For every patch, determining where the other patches overlap with it.
% % Deliberately using for loops to migrate the code more easily to C
% rowIdx = 1:repeatEvery(2):imageSize(2)-offset(2)-patchSize(2)+1;
% colIdx = 1:repeatEvery(1):imageSize(1)-offset(1)-patchSize(1)+1;
% templates = sparse(imageSize(1)*imageSize(2),numel(rowIdx)*numel(colIdx));
% templateInfo = zeros(2,2,numel(rowIdx)*numel(colIdx));
% %
% globalCounter = 1;
% colCounter = 1;
% for col=colIdx
%     rowCounter = 1;
%     for row=rowIdx
%         template = sparse(imageSize(1),imageSize(2));
%         template(row:row+patchSize(2)-1, col:col+patchSize(1)-1) = 1;
%         templates(:,globalCounter) = template(:);
%         templateInfo(1,:,globalCounter) = [rowCounter,colCounter];
%         templateInfo(2,:,globalCounter) = [row,col];
%         templateInfo(:,:,globalCounter);
% %         imagesc(template); pause
%         rowCounter = rowCounter + 1;
%         globalCounter = globalCounter + 1;        
%     end
%     colCounter = colCounter + 1;
% end
% theOtherImageBeginsAt = globalCounter;
% %
% colIdx = colIdx + offset(2);
% rowIdx = rowIdx + offset(1);
% for col=colIdx
%     rowCounter = 1;
%     for row=rowIdx
%         template = sparse(imageSize(1),imageSize(2));
%         template(row:row+patchSize(2)-1, col:col+patchSize(1)-1) = 1;
%         templates(:,globalCounter) = template(:);
%         templateInffo(1,:,globalCounter) = [rowCounter,colCounter];
%         templateInfo(2,:,globalCounter) = [row,col];
%         templateInfo(:,:,globalCounter);
%         %imagesc(template); pause
%         rowCounter = rowCounter + 1;
%         globalCounter = globalCounter + 1;        
%     end
%     colCounter = colCounter + 1;
% end
% % /way 2

%% way3
%  This pattern will rely more on inpainting and has less redundancy - less
%  overlaps
% patchSize = [16 16];    % height and width 
numPatches = [2,3];
repeatEvery = [8,8];  % Determines location of the patches which will compose the entire image. 
colIdx = 1:repeatEvery(2):numPatches(2)*repeatEvery(2);
rowIdx = 1:repeatEvery(1):numPatches(1)*repeatEvery(1);
imageSize = [rowIdx(end)+patchSize(1)-1, colIdx(end)+patchSize(2)-1]
templates = sparse(imageSize(1)*imageSize(2),numel(rowIdx)*numel(colIdx));
templateInfo = zeros(2,2,numel(rowIdx)*numel(colIdx));

% Specifying locations of seed patches, to be stored in 'toConsider'
locationsBool = true(numel(rowIdx)*numel(colIdx),1);
pattern = [1; 3; 5; 2; 4];
seedCols = repmat(pattern, [ceil(numel(rowIdx)/numel(pattern)),1]);
seedCols = seedCols(1:numel(rowIdx));
w = floor((numel(colIdx)-1)./5)
seedCols = [seedCols, ones(numel(rowIdx),w)*5];
seedCols = cumsum(seedCols,2);
seedRows = repmat((1:numel(rowIdx))',1,size(seedCols,2));
% locations = [seedRows(:), seedCols(:)];
locations = sub2ind([numel(rowIdx),numel(colIdx)],seedRows(:),seedCols(:));
locationsBool(locations) = false;


toConsider = true(numel(rowIdx)*numel(colIdx),1); % if true, the patch will get optimized. Otherwise, it won't
bla = zeros(imageSize);
globalCounter = 1;
colCounter = 1;
for col=colIdx
    rowCounter = 1;
    for row=rowIdx
        template = sparse(imageSize(1),imageSize(2));
        template(row:row+patchSize(2)-1, col:col+patchSize(1)-1) = 1;
        templates(:,globalCounter) = template(:);
        templateInfo(1,:,globalCounter) = [rowCounter,colCounter];
        templateInfo(2,:,globalCounter) = [row,col];
        bla = bla+template;
%         imagesc(template); pause
%
%         if(mod(colCounter,3)==1 && mod(rowCounter,3)==1 )
% %         if(false)
%             toConsider(globalCounter) = false;
% %             bla = bla+template;
%         end
        if(~locationsBool(globalCounter))
            toConsider(globalCounter) = false;
            bla = bla+template;
        end
        %
        rowCounter = rowCounter + 1;
        globalCounter = globalCounter + 1;
    end
    colCounter = colCounter + 1;
end

% % adding the "seed" patches
% offset = 6;
% repeatEvery = [24 24];
% rowIdx = offset+1:repeatEvery(2):imageSize(2)-patchSize(2)-1;
% colIdx = offset+1:repeatEvery(1):imageSize(1)-patchSize(2)-1;
% templates = [templates, sparse(imageSize(1)*imageSize(2),numel(rowIdx)*numel(colIdx))];
% colCounter = 1;
% for col=colIdx
%     rowCounter = 1;
%     for row=rowIdx
%         template = sparse(imageSize(1),imageSize(2));
%         template(row:row+patchSize(2)-1, col:col+patchSize(1)-1) = 1;
%         templates(:,globalCounter) = template(:);
%         templateInfo(1,:,globalCounter) = [rowCounter,colCounter];
%         templateInfo(2,:,globalCounter) = [row,col];
%         bla = bla+template;
%         toConsider(globalCounter) = false;
%         rowCounter = rowCounter + 1;
%         globalCounter = globalCounter + 1;
%     end
%     colCounter = colCounter + 1;
% end
% /way3


% Now we have all the templates (patch locations) vectorized. This way it's
% easy to find overlap of one template with the other, using matrix
% multiplication.
% Note that multipe channels at this point are still not taken into
% account.
% To reconstruct the original patch locations from their vectorized form,
% use it like this:
% i = 10; assert(i<numel(rowIdx)*numel(colIdx), 'Template index out of bounds. Select a smaller i.')
% reshape(templates(:,i),imageSize);
% templateInfo(1,:,i) holds the row,col index of patch i.
% templateInfo(2,:,i) holds the starting point of patch i in the global image space

% Finding where template 'i' overlaps with all the other templates, for all
% 'i's: overlaps = bsxfun(@times,templates(:,i),templates);
% Still, these indices are in the global image indexing space. What we are
% interested in is where the overlaps are in the patch indexing space
% (relative to the beginning of the patch). We'll use templateInfo(2,:,i)
% for that
S = cell(size(templates,2),size(templates,2),2);    % 2 matrices per {i,j} pair, describing one-to-one mapping in the overlaps. Matrix s in the problem formulation
for i=1:size(templates,2)
    overlaps = bsxfun(@times,templates(:,i),templates);
    ofInterest = any(overlaps);   % indices of patches which overlap with patch 'i'
    ofInterest(i) = false;
    ofInterest = find(ofInterest);
    % Removing overlaps with itself
    for j=ofInterest
        % Overlapping part in the global indexing space
        matricized = reshape(overlaps(:,j),[imageSize]);
        % Relative to the patch's 'i' position in the image
        matricized_i = matricized(...
            templateInfo(2,1,i):templateInfo(2,1,i)+patchSize(1)-1,...
            templateInfo(2,2,i):templateInfo(2,2,i)+patchSize(2)-1 );
        % Take channels into account - stack them one onto the other
        vectorized_i = repmat(matricized_i(:),[channels,1]);
        locations_i = find(vectorized_i);
        S1 = sparse(numel(vectorized_i),numel(vectorized_i));
        S1(sub2ind(size(S1),locations_i, locations_i)) = 1;
        % S1 is the transform matrix which will select only the vector
        % elements from 'locations' when used as S1*vectorized_1
        %        
        % Part relative to the patch's 'j' position in the image
        matricized_j = matricized(...
            templateInfo(2,1,j):templateInfo(2,1,j)+patchSize(1)-1,...
            templateInfo(2,2,j):templateInfo(2,2,j)+patchSize(2)-1 );

        % Take channels into account - stack them one onto the other
        vectorized_j = repmat(matricized_j(:),[channels,1]);
        locations_j = find(vectorized_j);
        % Now we have to determine one-to-one mapping from locations_j to
        % locations_i. For our specific case of patches being convex (square),
        % the order of appearance is sufficient to get the mapping.
        % Otherwise, you'll need to modify the code and keep note of which
        % pixel in j corresponds to which pixej in j.
        % In this case, one-to-one mapping is [locations_i, locations_j],
        % meaning that locations_j(i) maps to locations_i(i), for some
        % choice of i.
        S2 = sparse(numel(vectorized_j),numel(vectorized_j));
        S2(sub2ind(size(S2),locations_i, locations_j)) = 1;
        %
        % Finally, storing the S1 and S2
        S{i,j,1} = S1;
        S{i,j,2} = S2;
    end
%     %
%     % Showing how to properly place the template in global image
%     % indexing space from the local indexing space
%     inWholeImage = zeros(imageSize);
%     inWholeImage(templateInfo(2,1,i):templateInfo(2,1,i)+patchSize(1)-1,...
%         templateInfo(2,2,i):templateInfo(2,2,i)+patchSize(2)-1) = ...
%         reshape(vectorized_i,patchSize);
%     imagesc(inWholeImage); pause
end
%
% Another note of how S2 works: If 1 is placed on row(i,j), this means that
% location j will be moved to i (of course, sum of rows needs to be one,
% otherwise well get mixing/summation of the vector elements).

save overlaps

% % Trying ot the mapping
% v1 = [1,3,1,3,1,3,3]';
% v2 = [0,1.1,0,1.2,0,0,1.3]';
% % Let's say I want to measure divergence on locations 1,3,5
% locations = [1,3,5];
% S1 = zeros(7,7);
% S1(sub2ind(size(S1),locations, locations)) = 1;
% S1*v1
% onetoone = [1,2; 3,4; 5,7;]
% S2 = zeros(7,7);
% S2(sub2ind(size(S2),onetoone(:,1), onetoone(:,2))) = 1;
% % 'S2*v2' now maps v2 to match indexing of v1.
% % We can measure divergence between S1*v1 and S2*v2 - the required points
% % are aligned and otherwhere there are zeros, so divergence will not be
% % evaluated at these indeces.
% 
% 
% %% Another test
% T = [1, 0, 0;
%      1, 0, 0;
%      0, 1, 0;
%      0, 0, 1;
%      0, 0, 1;]
% v1 = [0 1 3]';
% v2 = [3 0 1]';
% % Suppose the two vectors overlap where the 3s are. This means that the
% % following hold
% S1 = zeros(5,5);
% locations = [4,5]
% S1(sub2ind(size(S1),locations, locations)) = 1;
% S1*T*v1
% onetoone = [4,1; 5,2];
% S2 = zeros(5,5);
% S2(sub2ind(size(S2),onetoone(:,1), onetoone(:,2))) = 1;
% S2*T*v2
% % Sweet. Now S1*T*v1 and S2*T*v2 can be matched.
