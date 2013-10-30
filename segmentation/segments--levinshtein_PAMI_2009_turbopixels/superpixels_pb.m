% Compute the superpixels given an image.
% img - the input image. Either RGB or Grayscale. Must be double in range
% 0..1
% numSuperpixels - number of superpixels
% display_int - Interval of frames at which the progress of the evolution
% will be displayed. 0 if not display is needed (default).
% contour_color - color of the superpixel boundaries (default is red)
%
% Returns:
%     sup_image - the final superpixel labeling.
%
function sup_image = superpixels_pb(img, numSuperpixels, display_int, contour_color, seeds, Pb)

    timeStep = 0.5;
    maxIterations = 500;
    
    if (nargin < 3 || isempty(display_int))
        display_int = 0;
    end
    
    if (nargin < 4 || isempty(contour_color))
        contour_color = [1,0,0];
    end
    
    if (nargin < 5 || isempty(seeds))
        seeds = [];
    end

    if (nargout < 4)
        phi = evolve_height_function_N(img, timeStep, maxIterations, 'superpixels', display_int, [], numSuperpixels, seeds, Pb);
    else
        [phi,frames] = evolve_height_function_N(img, timeStep, maxIterations, 'superpixels', display_int, [], numSuperpixels, seeds, Pb);
    end
    
    if (size(img,3) > 1)
        smooth_img = evolve_height_function_N(rgb2gray(img), 0.1, 10, 'curvature', 0, 0);
    else
        smooth_img = evolve_height_function_N(img, 0.1, 10, 'curvature', 0, 0);
    end

    speed2 = get_speed_based_on_gradient(Pb, [], 5, phi);
    boundary = get_superpixel_boundaries(phi,speed2);
    sup_image = get_segments_from_superpixel_boundaries(phi, speed2, boundary);
    %disp_img = display_logical(img, boundary, contour_color);
    