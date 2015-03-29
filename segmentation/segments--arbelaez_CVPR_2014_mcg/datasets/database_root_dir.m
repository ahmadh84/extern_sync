% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
%  University of California Berkeley (UCB) - USA
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  Pablo Arbelaez <arbelaez@berkeley.edu>
%  June 2014
% ------------------------------------------------------------------------ 
% This file is part of the MCG package presented in:
%    Arbelaez P, Pont-Tuset J, Barron J, Marques F, Malik J,
%    "Multiscale Combinatorial Grouping,"
%    Computer Vision and Pattern Recognition (CVPR) 2014.
% Please consider citing the paper if you use this code.
% ------------------------------------------------------------------------
% Adapt the paths to the place where you have downloaded PASCAL.
% ------------------------------------------------------------------------

function db_root_dir = database_root_dir( database )
if strcmp(database,'pascal2012')
    db_root_dir = '/data/images/everingham_IJCV_2010_pascalvoc';
elseif strcmp(database,'COCO')
    db_root_dir = '/data/images/lin_ECCV_2014_mscoco';
elseif strcmp(database,'SBD')
    db_root_dir = '/path/to/SBD/';
elseif strcmp(database,'bsds500')
    db_root_dir = '/path/to/BSDS500/';
else
    error(['Unknown database: ' database]);
end

end

