#!/usr/bin/env python

import sys
import htcondor
from htcondor import JobEventType
import os
from collections import defaultdict
from tabulate import tabulate
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    WHITE = '\033[97m'
    GRAY = '\033[90m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

directories = sys.argv[1:]

logs = []
for directory in directories:
    logs.extend([os.path.join(directory, x) for x in filter(lambda x: x.startswith("log_"), os.listdir(directory))])

jels = {}
for log in logs:
    name = os.path.basename(log).replace("log_","").replace(".txt","")
    jels[name] = htcondor.JobEventLog(log)

while True:
    try:
        table = []
        for log, j in jels.items():
            for event in j.events(stop_after=0):
                latest = event
            try:
                ret = latest["ReturnValue"]
            except KeyError:
                ret = "-"
            col = ""

            if ret == 0:
                log = bcolors.OKGREEN + log
                ret = f"{ret}" + bcolors.ENDC
            elif ret != "-":
                log = bcolors.FAIL + log
                ret = f"{ret}" + bcolors.ENDC

            table.append([log, latest.cluster, str(JobEventType.values[latest.type]), ret ])
        print(tabulate(table, headers=["Name", "Cluster", "Status", "Return"]))
    except KeyboardInterrupt:
        break
    break