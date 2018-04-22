% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% This file is part of the code available at
% https://github.com/iivek/sparse-synthextures
% which comes under GPL-3.0 license.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for selection = 1:6
    selection_ = selection; % TODO: 
    clear selection
    load(strcat('patches', num2str(selection_), '.mat'))
    selection = selection_;

    sparse_nmf
    selection_ = selection; % TODO: really ugly.
    clear selection
    save(strcat('dictionary', num2str(selection_), '.mat'))
    selection = selection_;
end