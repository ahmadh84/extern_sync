function [start_frames, end_frames] = generate_sequences(window_sz, ...
    num_frames, overlap, min_num_frames)
% Generates a set of windows given the size of the number of frames in a
% video sequence. You need to atleast specify the window size and the
% number of frames. The [] in the following list of arguments indicates the
% default value if not specified by the user
%
% @args:
%   window_sz: The number of frames in a window (the last window might have
%       less number of frames.
%
%   num_frames: The number of frames in the video sequence.
%
%   overlap [0]: The number of frames that need to overlap from one window to
%       the next. This is useful when your algorithm needs to have some
%       information overlap from the previous window.
%
%   min_num_frames [2]: Plus you can specify the number of frames that
%       should definitely present in each window. Otherwise that window is
%       dropped (usually happens with the last window). This is useful when
%       the information computed over each window needs to have a minimum
%       number of frames.

    if exist('overlap','var') ~= 1
        overlap = 0;
    end
    if exist('min_num_frames','var') ~= 1
        min_num_frames = 2;
    end
    
    assert(overlap < window_sz, ['The overlap size should be less than' ...
        ' the size of the window']);
    
    start_frames = 1:window_sz-overlap:num_frames;
    end_frames = window_sz:window_sz-overlap:num_frames;
    
    if isempty(end_frames) || end_frames(end) ~= num_frames
        end_frames(end+1) = num_frames;
    end
    if length(start_frames) > length(end_frames)
        start_frames = start_frames(1:length(end_frames));
    end
    
    assert(length(start_frames) == length(end_frames), ...
        'Number of sequences should be equivalent');
    
    % remove invalid size windows
    invalid_sz = end_frames - start_frames + 1 < min_num_frames;
    start_frames(invalid_sz) = [];
    end_frames(invalid_sz) = [];
end