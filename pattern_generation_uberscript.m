% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% This file is part of the code available at
% https://github.com/iivek/sparse-synthextures
% which comes under GPL-3.0 license.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for selection = 1:6
    close all
    selection_ = selection; % TODO: really ugly.
    clear selection
    load(strcat('dictionary', num2str(selection_), '.mat'));
    selection = selection_;
    
    generate_overlaps;
    pattern_generation
    
    selection_ = selection; % TODO: really ugly.
    clear selection
    save(strcat('result', num2str(selection_), '.mat'))
    selection = selection_;
end