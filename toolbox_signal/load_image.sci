function M = load_image(name, n0, options)

// read an image from a file
//
//  M = load_image(name, n0, options);
//
//  Set options.nbdims=3 for color/3D data.
//
// Copyright (c) 2008 Gabriel Peyre
 
if argn(2)<2
	n0 = [];
end

options.null = 0;

M = read_bin(name, options);

if not(isempty(n0)) & (n0~=size(M, 1) | n0~=size(M, 2)) & argn(2)>=2
	M = image_resize(M,n0);
end

M = squeeze(M);

endfunction