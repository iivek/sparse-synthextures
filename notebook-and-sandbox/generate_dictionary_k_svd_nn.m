%% Dictionary generation using K-SVD-NN
% Requires 'SelfSimSR-master/Lib/KSVD'
% The script assumes that X, the input matrix, exists in the workspace

extract_patches;

clear param
param.K = 256;
param.L = 10; 
param.numIteration = 200;
param.initialDictionary = rand(prod(patchSize),param.K);
param.InitializationMethod = 'GivenMatrix';
param.preserveDCAtom = 0;
param.displayProgress = 1;

warning('off','all')
[dictionary,output] = KSVD_NN(X,param);