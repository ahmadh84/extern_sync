function [speed,speed_x,speed_y,doublet,curvature] = get_speed_based_on_gradient_pb(Pb, doDoublet, normSigma, phi, speed, speed_x, speed_y)
    
    % Flag controlling relative vs. absolute gradient magnitude
    if (nargin < 2 || isempty(doDoublet))
        doDoublet = 0;
    end
    
    if (nargin < 3 || isempty(normSigma))
        normSigma = 5;
    end

    
    if (nargin < 5 || isempty(speed))
        % normalized response s.t.
        %   - we control the gradient mag that is mapped to half height
        %   - the result is mapped to [0,127]
        d = bwdist(Pb>0.05);
        speed = 1-exp(-d / normSigma);
        [speed_x,speed_y] = height_function_der(speed);
    end
    
    if (doDoublet)
        [dx,dy,dxx,dyy,dxy] = height_function_der(phi);
        
        speed = get_full_speed(dx,dy,dxx,dyy,dxy,speed,speed_x,speed_y,1);
    end
