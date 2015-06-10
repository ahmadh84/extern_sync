#!/usr/bin/env python

import os
import os.path
import pwd
import sys
import time
import errno

import tempfile
import shutil

#dm_dirpath = "/gfsai/ai-group/users/ahumayun"

def attempt_create_dir(out_dir):
    time.sleep(5)
    try:
        # create the folder if necessary
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    except IOError, e:
        pass

    # print(pwd.getpwuid(os.getuid()).pw_name)
    # output_test_fn = os.path.join(out_dir, 'test.out')
    # try:
    #     fp = open(output_test_fn, 'w')
    #     fp.write("10\n")
    #     fp.close()
    # except IOError, e:
    #     print(e.errno)
    #     print(e)


def get_ochs_trajectories(vid_filepath, skip_frames):
    # create temporary directory
    dirpath = tempfile.mkdtemp()

    # use ffmpeg to convert all frames to ppm images
    ffmpeg_cmd = "ffmpeg -i %s %s/%06d.ppm" % (vid_filepath, tmp_vid_fname)
    proc = subprocess.Popen(ffmpeg_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    ret = proc.wait()

    # write all frame filenames to a bmf text file

    # create executable command

    # read trajectories from file

    # convert the trajectories file into binary format

    # delete the temporary files/folders created
    shutil.rmtree(dirpath)



if __name__ == "__main__":
    #print sys.argv
    params = sys.argv[1:]
    if len(params) > 1:
        params = params[:1] + params[1].split()
    print(params)
    # find the output file parameter
    out_param_idx = [i for i,p in enumerate(params) if p[-4:] == '-out']
    # attempt to create the output directory
    if len(out_param_idx) > 0 and out_param_idx[0]+1 < len(params):
        out_dir = params[out_param_idx[0] + 1]
        out_dir = os.path.dirname(out_dir)
        attempt_create_dir(out_dir)

    params = [p[:2].replace("--", "-") + p[2:] for p in params]
    params = ' '.join(params)
    #print params
    ret = os.system("export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH; %s %s" % (dm_dirpath, os.path.join(dm_dirpath, "deepmatching"), params))
    if ret == 0:
        print("success")
        sys.exit(0)
    else:
        sys.exit(-1)
