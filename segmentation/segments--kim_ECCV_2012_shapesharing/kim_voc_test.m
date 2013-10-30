function kim_voc_test(SECTIONS, section_no, output_path)
    SetupPath;

    if ~exist('output_path', 'var'), output_path = '.'; end
    
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end
    
    if ~exist(fullfile(output_path, 'Scores'), 'dir')
        mkdir(fullfile(output_path, 'Scores'));
    end
    
    data_root = '/home/ahumayun/data/images/everingham_IJCV_2010_pascalvoc'; 
    
    data_dir = fullfile(data_root, 'JPEGImages');
    gt_dir = fullfile(data_root, 'SegmentationObject');

    voc_files = importdata(fullfile(data_root, 'ImageSets/Segmentation/val.txt'));
    
%     voc_files = {'2007_000032', '2007_000804', '2007_001311', ...
%                  '2007_001678', '2007_002105', '2007_002611', ...
%                  '2007_003189', '2007_003872', '2007_004537', ...
%                  '2007_004663', '2007_005173', '2007_006400', ...
%                  '2007_007165', '2007_007890', '2007_008140', ...
%                  '2007_008778', '2007_009348', '2008_000725', ...
%                  '2008_001610', '2008_002240', '2008_002425', ...
%                  '2008_002929', '2008_003329', '2008_004575', ...
%                  '2008_005197', '2008_005680', '2008_006254', ...
%                  '2008_006874', '2008_007313', '2008_008711', ...
%                  '2009_000603', '2009_002346', '2009_002382', ...
%                  '2009_002928', '2009_004656', '2009_004886', ...
%                  '2010_001851', '2010_003088', '2010_004951', ...
%                  '2010_005758', '2007_009084', '2010_000238', ...
%                  '2010_002868', '2010_003781'};
    
    cpmc_root = '/home/ahumayun/videovolumes/cpmc_src';
    addpath(fullfile(cpmc_root, 'code'));
    addpath('/home/ahumayun/videovolumes/extern_src');
    
    % find what files to work on
    if ~exist('SECTIONS', 'var'), SECTIONS = 1; end
    if ~exist('section_no', 'var'), section_no = 1; end
    [start_idx, end_idx] = give_section(length(voc_files), ...
        SECTIONS, section_no);
    
    fprintf('Computing from %d to %d\n', start_idx, end_idx);
    
    for voc_idx = start_idx:end_idx
        img_name = voc_files{voc_idx};
        img_filepath = fullfile(data_dir, [img_name, '.jpg']);
        
        [masks, timing] = ComputeSegment(img_filepath);
        additional_info.timing = timing;
        
        save(fullfile(output_path, [img_name '.mat']), ...
             'masks', 'additional_info');
         
        % compute scores; print results; save
        [Q, collated_scores] = SvmSegm_segment_quality(img_name, gt_dir, ...
                                                       masks, 'overlap');
        fprintf('\n');
        
        num_segs = size(masks,3);
        seg_time = timing.t_all;
        
        save(fullfile(output_path, 'Scores', [img_name '_scores.mat']), ...
             'Q', 'collated_scores', 'num_segs', 'seg_time');
        
        fprintf('>>>>> %s: Overlap %.4f\n\n', img_name, ...
                collated_scores.avg_best_overlap);
    end
end


function [start_idx, end_idx, section_sz] = give_section(total_len, ...
        SECTIONS, section_no)
    % get which params this program run will compute
    section_sz = ceil(total_len / SECTIONS);
    start_idx = section_sz * (section_no-1) + 1;
    end_idx = section_sz * section_no;
    if end_idx > total_len
        end_idx = total_len;
    end
    section_sz = end_idx - start_idx + 1;
end
