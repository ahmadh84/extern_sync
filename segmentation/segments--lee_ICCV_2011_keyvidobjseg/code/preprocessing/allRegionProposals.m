function allRegionProposals
num_tries = 10;

data_dir = '~/data/Images/SegTrack_201102';

videoName = 'birdfall2';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'cheetah';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'girl';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'monkeydog';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'parachute';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'penguin';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'penguin_rszd';
computeRegionProposals(videoName, data_dir, num_tries);


data_dir = '~/data/Images/SegTrack_2013';

videoName = 'bird_of_paradise_rszd';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'bmx_rszd';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'drift_rszd';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'frog';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'hummingbird_rszd';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'monkey';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'soldier';
computeRegionProposals(videoName, data_dir, num_tries);

videoName = 'worm';
computeRegionProposals(videoName, data_dir, num_tries);
end

