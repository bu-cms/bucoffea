#!/usr/bin/env python

from bucoffea.helpers.cutflow import print_cutflow
import sys
from coffea.util import load
print_cutflow(load(sys.argv[1]))