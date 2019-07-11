#!/usr/bin/env python

from bucoffea.monojet.monojetProcessor import debug_print_cutflows
import sys
from coffea.util import load
debug_print_cutflows(load(sys.argv[1]))