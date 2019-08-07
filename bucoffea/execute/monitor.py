#!/usr/bin/env python

import sys
import htcondor
from htcondor import JobEventType
import os
from collections import defaultdict

from tabulate import tabulate
import curses

import time
import logging

logger = logging.getLogger( 'monitor' )
format = '%(levelname)s (%(name)s) [%(asctime)s]: %(message)s'
date = '%F %H:%M:%S'
logging.basicConfig( level='DEBUG', format=format, datefmt=date, filename='monitor_log.txt')

INTERVAL=5

def main(stdscr):
    curses.start_color()
    curses.use_default_colors()
    curses.init_pair(1, curses.COLOR_GREEN, -1)
    curses.init_pair(2, curses.COLOR_RED, -1)
    curses.init_pair(3, curses.COLOR_WHITE, -1)

    directories = sys.argv[1:]

    def read_logs():
        logs = []
        for directory in directories:
            for path, _, files in os.walk(directory):
                logs.extend([os.path.join(path, x) for x in files if x.startswith("log_")])

        jels = {}
        for log in logs:
            name = os.path.basename(log).replace("log_","").replace(".txt","")
            try:
                jels[name] = htcondor.JobEventLog(log)
            except:
                logger.warn(f'Could not build JobEventLog for log: {log}')
                # raise
                continue
        return jels

    # Initiate new pad
    jels = read_logs()
    logger.info(f'Found {len(jels)} logs.')
    padlen = len(jels)+5
    pad = curses.newpad(padlen,curses.COLS)
    stdscr.nodelay(True)
    stdscr.timeout(0)
    pad.nodelay(True)
    pad_pos = 0


    while True:
        pad.refresh( pad_pos, 0, 1, 1, curses.LINES-1, curses.COLS-1)
        jels = read_logs()

        colors = [3,3]
        try:
            table = []
            counter = 1
            for log, j in sorted(jels.items()):
                first = None
                for event in j.events(stop_after=0):
                    if not first:
                        first = event
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
                table.append([counter, log, latest.cluster, str(JobEventType.values[latest.type])])
                j.close()
                counter += 1
            tab = tabulate(sorted(table), headers=["","Name", "Cluster", "Status", "Return"])
            for i,l in enumerate(tab.split('\n')):
                pad.addstr(i,0,l, curses.color_pair(colors[i]))


            # stdscr.refresh()
            # Scrolling
            start = time.time()
            stop = start
            while stop-start < INTERVAL:
                stop = time.time()

                pad.refresh( pad_pos, 0, 1, 1, curses.LINES-1, curses.COLS-1)
                cmd = stdscr.getch()
                if cmd == 'q':
                    sys.exit(0)
                elif  cmd == curses.KEY_HOME:
                    pad_pos = 0
                elif  cmd == curses.KEY_END:
                    pad_pos = len(jels) - curses.LINES + 1
                elif  cmd == curses.KEY_NPAGE:
                    pad_pos += 25
                elif  cmd == curses.KEY_PPAGE:
                    pad_pos -= 25
                elif  cmd == curses.KEY_DOWN:
                    pad_pos += 1
                elif cmd == curses.KEY_UP:
                    pad_pos -= 1

                # Post process
                if pad_pos < 0:
                    pad_pos = 0
                if pad_pos >= padlen:
                    pad_pos = padlen-1

                time.sleep(0.01)
        except KeyboardInterrupt:
            break

if __name__ == "__main__":
    curses.wrapper(main)