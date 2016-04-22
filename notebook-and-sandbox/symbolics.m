%% l1/l2 measure of sparsity
M=5;
factors = 4;

v = sym('a', [factors,M]);
assume(v>0)
objective = sum(sum(v,1)./sqrt(sum(v.^2,1)));

for i=1:M
    for j=1:factors
        iWantThis(j,i) = diff(objective,v(j,i));
    end
end

negative = v.*repmat( sum(v,1)./sqrt(sum(v.^2,1)).^3, [factors,1]);
positive = 1./repmat(sqrt(sum(v.^2,1)), [factors,1]); 
iHaveThis = positive-negative;    
simplify(iHaveThis-iWantThis)