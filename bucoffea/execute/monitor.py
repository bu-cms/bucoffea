#!/usr/bin/env python

import sys
import htcondor
from htcondor import JobEventType
import os
from collections import defaultdict

from tabulate import tabulate
import curses

import time

def main(stdscr):
    curses.start_color()
    curses.use_default_colors()
    curses.init_pair(1, curses.COLOR_GREEN, -1)
    curses.init_pair(2, curses.COLOR_RED, -1)
    curses.init_pair(3, curses.COLOR_WHITE, -1)
    directories = sys.argv[1:]

    logs = []
    for directory in directories:
        for path, _, files in os.walk(directory):
            logs.extend([os.path.join(path, x) for x in files if x.startswith("log_")])

    jels = {}
    for log in logs:
        name = os.path.basename(log).replace("log_","").replace(".txt","")
        jels[name] = htcondor.JobEventLog(log)

    while True:
        stdscr.clear()
        colors = [3,3]
        try:
            table = []
            for log, j in jels.items():
                for event in j.events(stop_after=0):
                    latest = event
                try:
                    ret = latest["ReturnValue"]
                except KeyError:
                    ret = "-"
                if ret == 0:
                    colors.append(1)
                elif ret != "-":
                    colors.append(2)
                else:
                    colors.append(3)
                table.append([log, latest.cluster, str(JobEventType.values[latest.type]), ret ])
            tab = tabulate(table, headers=["Name", "Cluster", "Status", "Return"])
            for i,l in enumerate(tab.split('\n')):
                stdscr.addstr(i,0,l, curses.color_pair(colors[i]))
            stdscr.refresh()
            stdscr.nodelay(True)
        except KeyboardInterrupt:
            break
        time.sleep(5)

if __name__ == "__main__":
    curses.wrapper(main)