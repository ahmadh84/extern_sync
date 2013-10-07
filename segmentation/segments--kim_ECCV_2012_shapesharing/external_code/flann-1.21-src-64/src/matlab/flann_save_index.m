%Copyright 2008-2009  Marius Muja (mariusm@cs.ubc.ca). All rights reserved.
%Copyright 2008-2009  David G. Lowe (lowe@cs.ubc.ca). All rights reserved.

%THE BSD LICENSE
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions
%are met:
%
%1. Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
%2. Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in the
%   documentation and/or other materials provided with the distribution.
%
%THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
%IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
%OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
%IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
%INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
%NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
%THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%This function is added by Jaechul Kim (jaechul@cs.utexas.edu) as a matlab
%wrapper for c++ function, flann_save_index.

function rc = flann_save_index(index, filename)
%FLANN_SAVE_INDEX  Save the nearest-neighbors index
%Jaechul Kim, August 2009
%input: index - FLANN type(uint32 or uint64): nearest-neighbors index
%input: filename - String type: filename for saving index
%output: rc - return flag: 0-success, otherwise-fail

rc = nearest_neighbors('save_index',index, filename);