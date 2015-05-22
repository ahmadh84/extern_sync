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
            byts = fd.read(4*6)
            matches.append(list(struct.unpack('<ffffff', byts)))
    return matches

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
        print("%s - %.2fs" % (ret, delta.seconds + delta.microseconds/1E6))
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
    def cmp_matches(self, matches1, matches2):
        self.assertEqual(len(matches1), len(matches2), "The number of matches is not equal: %d != %d" % (len(matches1), len(matches2)))
        for m1, m2 in zip(matches1, matches2):
            self.assertEqual(m1[0], m2[0], "x1 coordinate not equal")
            self.assertEqual(m1[1], m2[1], "y1 coordinate not equal")
            self.assertEqual(m1[2], m2[2], "x2 coordinate not equal")
            self.assertEqual(m1[3], m2[3], "y2 coordinate not equal")
            self.assertTrue(math.fabs(m1[4] - m2[4]) < 1e-5, "maxima not equal")
            self.assertTrue(math.fabs(m1[5] - m2[5]) < 1e-5, "score not equal")
    
    def test01_binary(self):
        """tests binary against to ascii file output
        """
        # compute deep matches for two images
        matches1 = deep_matches_ims("./deepmatching", "liberty1.png", "liberty2.png", True)
        matches2 = deep_matches_ims("./deepmatching", "liberty1.png", "liberty2.png", False)
        self.cmp_matches(matches1, matches2)

    def test02_static(self):
        """tests output against statically compiled binary
        """
        # compute deep matches for two images
        matches1 = deep_matches_ims("./deepmatching", "liberty1.png", "liberty2.png", True)
        matches2 = deep_matches_ims("./deepmatching-static", "liberty1.png", "liberty2.png", False, False)
        self.cmp_matches(matches1, matches2)

if __name__ == "__main__":
    """Executes the unit tests in ParserUnitTest"""
    suite = unittest.TestLoader().loadTestsFromTestCase(ParserUnitTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
