#!/usr/bin/env python

import os
import sys
import tempfile
import struct
import unittest
import math
from datetime import datetime

def read_output_matches_binary(tmp_fname):
    matches = []
    # read output
    with open(tmp_fname, "rb") as fd:
        byts = fd.read(4)
        num_matches = struct.unpack('<i', byts)[0]
        for x in range(num_matches):
            byts = fd.read(2*4 + 4*2)
            matches.append(list(struct.unpack('<HHHHff', byts)))
    return matches

def read_output_vid_matches_binary(tmp_fname):
    all_matches = []
    # read output
    with open(tmp_fname, "rb") as fd:
        byts = fd.read(8)
        tmp = struct.unpack('<ii', byts)
        num_frames = tmp[0]
        skip_frames = tmp[1]

        for fi in xrange(1, num_frames-(skip_frames), skip_frames+1):
            matches = []
            byts = fd.read(4)
            frame_no = struct.unpack('<i', byts)[0]
            byts = fd.read(4)
            num_matches = struct.unpack('<i', byts)[0]
            for x in range(num_matches):
                byts = fd.read(2*4 + 4*2)
                matches.append(list(struct.unpack('<HHHHff', byts)))
            all_matches.append(matches)

    return all_matches

def read_output_matches(tmp_fname, has_num_matches=True):
    matches = []
    # read output
    with open(tmp_fname) as fd:
        lines = fd.readlines()
        if has_num_matches:
            num_matches = int(lines[0])
            del lines[0]

        for l in lines:
            matches.append([float(x) for x in l.split()])

        if has_num_matches:
            assert(len(matches) == num_matches)
    return matches

def deep_matches_vids(exec_name, vid_frames_path, extra_params=""):
    all_matches = []
    # get a temporary filename
    tmp_vid_fname = next(tempfile._get_candidate_names()) + ".avi"
    tmp_fname = next(tempfile._get_candidate_names())
    ffmpeg_cmd = "ffmpeg -i %s -c:v huffyuv -pix_fmt rgb24 %s" % (vid_frames_path, tmp_vid_fname)
    cmd = '%s -i "%s" %s -b -out %s' % (exec_name, tmp_vid_fname, extra_params, tmp_fname)
    try:
        ret = os.system(ffmpeg_cmd)
        t1 = datetime.now()
        ret = os.system(cmd)
        t2 = datetime.now()
        delta = t2 - t1
        #sys.stdout.write("%s - %.2fs | " % (ret, delta.seconds + delta.microseconds/1E6))
        all_matches = read_output_vid_matches_binary(tmp_fname)
    finally:
        # delete the output file if it exists
        try:
            os.remove(tmp_vid_fname)
            os.remove(tmp_fname)
        except OSError:
            pass

    return all_matches

def deep_matches_ims(exec_name, im1_fp, im2_fp, binary_out=False, has_num_matches=True, extra_params=""):
    matches = []
    # get a temporary filename
    tmp_fname = next(tempfile._get_candidate_names())
    if binary_out:
        extra_params += " -b"
    cmd = '%s "%s" "%s" %s -out %s' % (exec_name, im1_fp, im2_fp, extra_params, tmp_fname)
    try:
        t1 = datetime.now()
        ret = os.system(cmd)
        t2 = datetime.now()
        delta = t2 - t1
        #sys.stdout.write("%s - %.2fs | " % (ret, delta.seconds + delta.microseconds/1E6))
        if binary_out:
            matches = read_output_matches_binary(tmp_fname)
        else:
            matches = read_output_matches(tmp_fname, has_num_matches)
    finally:
        # delete the output file if it exists
        try:
            os.remove(tmp_fname)
        except OSError:
            pass

    return matches

class ParserUnitTest(unittest.TestCase):
    """Unit test class to check the functionality of parse_file_for_jobs()
    """
    cp = "-nt 1"

    def cmp_matches(self, matches1, matches2):
        self.assertEqual(len(matches1), len(matches2), "The number of matches is not equal: %d != %d" % (len(matches1), len(matches2)))
        for m1, m2 in zip(matches1, matches2):
            self.assertEqual(m1[0], m2[0], "x1 coordinate not equal %d != %d" % (m1[0], m2[0]))
            self.assertEqual(m1[1], m2[1], "y1 coordinate not equal %d != %d" % (m1[1], m2[1]))
            self.assertEqual(m1[2], m2[2], "x2 coordinate not equal %d != %d" % (m1[2], m2[2]))
            self.assertEqual(m1[3], m2[3], "y2 coordinate not equal %d != %d" % (m1[3], m2[3]))
            self.assertTrue(math.fabs(m1[4] - m2[4]) < 1e-1, "maxima not equal %f != %f" % (m1[4], m2[4]))
            self.assertTrue(math.fabs(m1[5] - m2[5]) < 4, "score not equal %d != %d" % (m1[5], m2[5]))
    
    def run_video_test(self, seq_path, num_frames, skip_frames):
        all_matches = deep_matches_vids("./deepmatching", seq_path, extra_params=self.cp + " -k %d" % skip_frames)
        all_matches2 = []
        for fi in xrange(1, num_frames-skip_frames, skip_frames+1):
            frame1_num = fi
            frame2_num = fi + skip_frames + 1
            im1_fp = seq_path % frame1_num
            im2_fp = seq_path % frame2_num
            matches = deep_matches_ims("./deepmatching", im1_fp, im2_fp, True, extra_params=self.cp)
            all_matches2.append(matches)

        self.assertEqual(len(all_matches), len(all_matches2), "The number of match sets is unequal")
        for m1, m2 in zip(all_matches, all_matches2):
            self.cmp_matches(m1, m2)

    def test01_binary(self):
        """tests binary against to ascii file output
        """
        # compute deep matches for two images
        matches1 = deep_matches_ims("./deepmatching", "liberty1.png", "liberty2.png", True, extra_params=self.cp)
        matches2 = deep_matches_ims("./deepmatching", "liberty1.png", "liberty2.png", False, extra_params=self.cp)
        self.cmp_matches(matches1, matches2)

    def test02_binary(self):
        """tests binary against to ascii file output on dino
        """
        # compute deep matches for two images
        matches1 = deep_matches_ims("./deepmatching", "dino1.jpg", "dino2.jpg", True, extra_params=self.cp)
        matches2 = deep_matches_ims("./deepmatching", "dino1.jpg", "dino2.jpg", False, extra_params=self.cp)
        self.cmp_matches(matches1, matches2)

    def test03_binary(self):
        """tests binary against to ascii file output on climb
        """
        # compute deep matches for two images
        matches1 = deep_matches_ims("./deepmatching", "climb1.png", "climb2.png", True, extra_params=self.cp)
        matches2 = deep_matches_ims("./deepmatching", "climb1.png", "climb2.png", False, extra_params=self.cp)
        self.cmp_matches(matches1, matches2)

    def test04_static(self):
        """tests output against statically compiled binary
        """
        # compute deep matches for two images
        matches1 = deep_matches_ims("./deepmatching", "climb1.png", "climb2.png", True, extra_params=self.cp)
        matches2 = deep_matches_ims("./deepmatching-static", "climb1.png", "climb2.png", False, False, extra_params=self.cp)
        self.cmp_matches(matches1, matches2)

    def test05_video(self):
        """tests matches produced from a video
        """
        seq_path = "test_ims/soldier/soldier_%03d.png"
        self.run_video_test(seq_path, 11, 0)

    def test06_video(self):
        """tests matches produced from a video with 1 skipped frame
        """
        seq_path = "test_ims/soldier/soldier_%03d.png"
        self.run_video_test(seq_path, 11, 1)

    def test07_video(self):
        """tests matches produced from a video with 2 skipped frames
        """
        seq_path = "test_ims/soldier/soldier_%03d.png"
        self.run_video_test(seq_path, 11, 2)


if __name__ == "__main__":
    """Executes the unit tests in ParserUnitTest"""
    suite = unittest.TestLoader().loadTestsFromTestCase(ParserUnitTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
