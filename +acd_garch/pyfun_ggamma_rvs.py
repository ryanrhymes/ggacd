#!/usr/bin/env python

import os
import sys
from scipy.stats import gengamma

if __name__=='__main__':
    shape = float(sys.argv[1])
    power = float(sys.argv[2])
    numrv = int(sys.argv[3])

    for i in xrange(numrv):
        rvs = gengamma.rvs(shape, power)
        print "%.15f" % rvs

    sys.exit(0)
