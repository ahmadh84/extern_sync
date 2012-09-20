function M = read_bmp(name)

// read_bmp - read a 8 bit bw binary file
//
//   M = read_bmp(name);
//
//	Code of Kevin Obrejan

  if strindex(name, '.') == []
      name = strcat([name '.bmp']);
  end

  fid = mopen(name, 'rb');
  if fid<0
      error(['File ' name ' does not exists.']);
  end

  // sizes
  tmp = mget(2, 'ucl', fid); // poubelle
  tmp = mget(6, 'uil', fid);
  fileSize = tmp(1);
  offset = tmp(3);
  width = tmp(5);
  height = tmp(6);
  tmp = mget(2, 'usl', fid); // poubelle
  if tmp(2) ~= 24
    error(strcat(['File ' name ' must be a 24-bits bitmap and not a ' string(tmp(2)) '-bits']));
  end

  mseek(offset, fid);
  data = mget(fileSize-offset, 'ucl', fid);
  mclose(fid);

  lineWidth = 3*width+pmodulo(-3*width, 4);
  data = matrix(data, [lineWidth, height]);
  data = data(1:3*width, :);
  data = matrix(data(:), [3, width*height]);
  M = [];
  M = flipdim(matrix(floor((data(1, :) + data(2, :) + data(3, :))/3), [width, height])', 1);
endfunction