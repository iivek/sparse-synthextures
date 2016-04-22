% Mali test min D(T*v|y) gdje je divergencija simetricna beta divergencija.
% Spustanjem po gradijenu po V optimiram gornji kriterij
%
% Napomena: Imam problema sa beta = 1 i beta = 0. To su special cases za
% koje imam problema sa gradijentom. Konkretno, za symmetric KL divergence.
% Cini se da normalizacja redaka pomaze - vidi http://www.uta.edu/faculty/rcli/Teaching/math5392/papers/bblpp2007.pdf
% Cini se cak da simetricna beta divergencija nije konveksna gdje i obicna
% Na svu srecu kl divergencija jeste.
%
% Imas izraze za l1 i l2 normalizaciju. ako koristis L1 normalizaciju, nema
% smisla penalizirati l1 normu kako bi dobio sparsity - jer neces. Tada
% treba penalizirati omjer l1 i l2 norme.
iterations = 1000;
epsilon = 1e-12;
beta = 2;
lambda = 0.1; % sparsity penalization

height = 20;
width = 1;
factors = 4;
Ttrue = rand(height,factors);
vtrue = rand(factors,width);
% vtrue = vtrue./repmat(sqrt(sum(vtrue.^2,1)),[factors,1]);    % unit l2 norm accross cols
% vtrue(1:2) = 0;
vtrue = vtrue./repmat(sum(vtrue,1),[factors,1]);    % unit l1 norm accross cols
y = Ttrue*vtrue;

T = rand(size(Ttrue));
%T = Ttrue;
v = rand(factors,width);

% switch beta
%     case 0,
%         % IS divergence
%     case 1,
%         % KL divergence 
%     otherwise
% end

obj = Inf;
for i=1:iterations    
%     % Unit l2 norm
%     % We can penalize L1 norm for promoting sparsity
%     negative = T'*((beta-1).*y.*(T*v).^(beta-2)+y.^(beta-1));
%     positive = beta.*T'*(T*v)+lambda;
%     numerator = negative + v.*repmat(sum(v.*positive,1),[factors,1]);
%     denominator = positive + v.*repmat(sum(v.*negative,1), [factors,1]);
%     v = v.*numerator./denominator;
%     v = v./repmat(sqrt(sum(v.^2,1)),[factors,1]);    % unit l2 norm accross cols
%     % objective:
%     curr = 1./(beta-1).*sum( (T*v-y).*((T*v).^(beta-1)-y.^(beta-1)) );
%     curr = curr + lambda*sum(v(:));
%     % /Unit l2 norm
    %
    % Unit l1 norm
    % In this case l1 penalization won't work. But we can penalize ratio of
    % l1 and l2 to promote sparsity (direct connections with Hoyer's sparsity)
    negative = T'*((beta-1).*y.*(T*v).^(beta-2)+y.^(beta-1))...
        + lambda.*v.*norm(v,1)./norm(v,2).^3;
    positive = beta.*T'*(T*v)...
        + lambda./norm(v,2);
    numerator = negative + repmat(sum(v.*positive,1),[factors,1]);
    denominator = positive + repmat(sum(v.*negative,1), [factors,1]); 
    v = v.*numerator./denominator;
    v = v./repmat(sum(v,1),[factors,1]);    % unit L1 norm accross cols
    % objective:
    curr = 1./(beta-1).*sum( (T*v-y).*((T*v).^(beta-1)-y.^(beta-1)) );
    curr = curr + lambda*norm(v,1)./norm(v,2);
    % /unit l1 norm
            
            
    %    
    obj = [obj curr];
    if obj(end)>obj(end-1)
        'broken'
        obj(end)-obj(end-1)
        %        break
    end
    subplot(2,1,1); plot([T*v,y]);
    subplot(2,1,2); plot(obj);
    drawnow
   
end

% Perfektno!







