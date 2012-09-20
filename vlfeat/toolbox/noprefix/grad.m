function varargout = grad(varargin)
% VL_GRAD Compute the gradient of an image
%   [IX,IY] = VL_GRAD(I) returns the finite gradient components IX,IY
%   of the grayscale image I. I must be a two-dimensional matrix. The
%   function uses central differences for all but the boundary pixels,
%   for which it uses forward/backward differences as appropriate.
%
%   The function accepts the following options:
%
%   Type:: Central
%     Use either 'central', 'forward', or 'backward' differences for
%     all but the boundary pixels.
%
%   See also: GRADIENT(), VL_HELP().
[varargout{1:nargout}] = vl_grad(varargin{:});
