%% Playing aroung with a small test dictionary, consisting of horizontal,
vertical and diagonal elements
template = cell(4,1);
template{1} = zeros(16,16); % left
template{1}([8,9],[1:9]) = 1;
template{2} = fliplr(template{1}); % right
template{3} = template{1}'; % up
template{4} = flipud(template{3});  % down
%
% Some blurring
h = 1/3*ones(3,1);
H = h*h';
template{1} = filter2(H,template{1});
template{2} = filter2(H,template{2});
template{3} = filter2(H,template{3});
template{4} = filter2(H,template{4});

T = [];
idx = 1;
for i=1:4    
    combis = combnk(1:4,i);
    for j=1:size(combis,1)
        T(:,:,idx) = zeros(16,16);
        for k=1:size(combis,2)
            T(:,:,idx) = T(:,:,idx) | template{combis(j,k)};
        end
        idx = idx+1;
    end
end

cnt = 16;
for i=1:15
    T(:,:,cnt) = circshift(T(:,:,i),[0,8]);
    cnt = cnt+1;
    T(:,:,cnt) = circshift(T(:,:,i),[8,0]);
    cnt = cnt+1;
end

T = reshape(T,16*16,[]);
% Dictionary elements are now accross cols. To reshape element i back to
% image patch, invoke: 'reshape(T(:,i),16,16)'