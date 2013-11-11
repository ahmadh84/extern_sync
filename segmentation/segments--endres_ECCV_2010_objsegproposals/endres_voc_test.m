function endres_voc_test(SECTIONS, section_no, output_path)
    proposals_startup;
    
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
        
        [proposal_data, superpixels, seg_time] = generate_only_segments(img_filepath);
        additional_info.timing.t_all = seg_time;
        
        % convert superpixel segments to masks
        masks = false([size(superpixels), length(proposal_data.final_regions)]);
        for seg_idx = 1:length(proposal_data.final_regions)
            masks(:,:,seg_idx) = ismember(superpixels, proposal_data.final_regions{seg_idx});
        end
        
        save(fullfile(output_path, [img_name '.mat']), ...
             'masks', 'additional_info', '-v7.3');
         
        % compute scores; print results; save
        [Q, collated_scores] = SvmSegm_segment_quality(img_name, gt_dir, ...
                                                       masks, 'overlap');
        fprintf('\n');
        
        num_segs = size(masks,3);
        
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