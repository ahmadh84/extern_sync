#!/usr/bin/env python

import os
import os.path
import pwd
import sys
import errno

dm_dirpath = "/gfsai/ai-group/users/ahumayun"

def attempt_create_dir(out_dir):
    # create the folder if necessary
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # print(pwd.getpwuid(os.getuid()).pw_name)
    # output_test_fn = os.path.join(out_dir, 'test.out')
    # try:
    #     fp = open(output_test_fn, 'w')
    #     fp.write("10\n")
    #     fp.close()
    # except IOError, e:
    #     print(e.errno)
    #     print(e)


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
    ret = os.system("export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH; %s %s" % 
                    (dm_dirpath, os.path.join(dm_dirpath, "deepmatching"), params)
