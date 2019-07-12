#!/usr/bin/env python

import sys
import htcondor
from htcondor import JobEventType
import os
from collections import defaultdict
from tabulate import tabulate

directories = sys.argv[1:]

logs = []
for directory in directories:
    logs.extend([os.path.join(directory, x) for x in filter(lambda x: x.startswith("log_"), os.listdir(directory))])

jels = {}
for log in logs:
    name = os.path.basename(log).replace("log_","").replace(".txt","")
    jels = { name : htcondor.JobEventLog(log) for log in logs}

while True:
    try:
        table = []
        for log, j in jels.items():
            for event in j.events(stop_after=0):
                latest = event
            try:
                ret = event["ReturnValue"]
            except KeyError:
                ret = "-"
            table.append([log, latest.cluster, str(JobEventType.values[latest.type]) ])
        print(tabulate(table, headers=["Name", "Cluster", "Status", "Return"]))
    except KeyboardInterrupt:
        break
    break