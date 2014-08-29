% Canidate Regions - version 1.0.  Code distributed under LGPL. 
% August 2014, University of California, Los Angeles
%
% Please, cite this work if you use the code: 
%   B. Bonev, A.L. Yuille, "A Fast and Simple Algorithm for Producing Candidate Regions", 
%   ECCV 2014, Zürich, Switzerland, 6th-12th September, 2014

classdef imageHandle < handle
    properties
        img; % Original in RGB
        sing;% RGB single
        labl;
        laba;
        labb;
        gray;% Grayscale single
        gy;  % Gradients
        gx;
        gyy;
        gxx;
        RC;  % [Rows, Columns]
        edg;
    end
    methods
        function obj = imageHandle(image, edges)
            if ischar(image)
                obj.img  = imread(image);
            else
                obj.img = image;
            end
            obj.sing = im2single(obj.img);
            %run /home/boyan/prog/ucla/vlfeat-0.9.16/toolbox/vl_setup
            if size(obj.img,3)==3
                %obj.lab = vl_xyz2lab(vl_rgb2xyz(obj.img));
                hsv = rgb2hsv(obj.img);
                obj.labl = hsv(:,:,1);
                obj.laba = hsv(:,:,2);
                obj.labb = hsv(:,:,3);
                obj.gray = single(rgb2gray(obj.img));
            else
                obj.gray = obj.sing;
                obj.labl = obj.gray;
                obj.laba = obj.gray;
                obj.labb = obj.gray;
            end

            %gradients:
            tmpgray = [obj.gray(1:end,1),   obj.gray(1:end,1), ...
                       obj.gray, ...
                       obj.gray(1:end,end), obj.gray(1:end,end)];
            tmpgray = [tmpgray(1,1:end);   tmpgray(1,1:end); ...
                       tmpgray; ...
                       tmpgray(end,1:end); tmpgray(end,1:end)];
            obj.gy  = conv2(tmpgray,[ 1 2 1; 0 0 0; -1 -2 -1]','same');
            obj.gx  = conv2(tmpgray,[ 1 2 1; 0 0 0; -1 -2 -1],'same');
            obj.gyy = conv2(obj.gy,[ 1 2 1; 0 0 0; -1 -2 -1]','same');
            obj.gxx = conv2(obj.gx,[ 1 2 1; 0 0 0; -1 -2 -1],'same');
            obj.gy  = obj.gy(3:end-2,3:end-2);
            obj.gx  = obj.gx(3:end-2,3:end-2);
            obj.gyy  = obj.gyy(3:end-2,3:end-2);
            obj.gxx  = obj.gxx(3:end-2,3:end-2);

            [ROWS,COLS,~]=size(obj.img); obj.RC=[ROWS,COLS];
            
            %edges:
            if exist('edges','var')
                obj.edg = edges;
            else
                obj.edg = 0.1*edge(obj.gray,'sobel',20) + ...
                    0.2*edge(obj.gray,'sobel',25) + ...
                    0.3*edge(obj.gray,'sobel',30) + ...
                    0.5*edge(obj.gray,'sobel',35);
                se = strel('square',2);
                %se = strel('disk',0);
                obj.edg = imdilate(obj.edg,se);
            end
        end
    end
end
