#!/usr/bin/env python

from bucoffea.helpers.cutflow import print_cutflow
import sys
from coffea.util import load

acc = None
for argument in sys.argv:
    if argument.endswith('.coffea'):
        tmp = load(argument)
        acc = tmp if not acc else acc + tmp
print_cutflow(acc)