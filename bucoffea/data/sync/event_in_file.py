#!/usr/bin/env python
import numpy as np
import sys
import argparse
import uproot
import tqdm
def commandline():
    parser = argparse.ArgumentParser(prog='Check which file a given event is in.')
    parser.add_argument('files', type=str, help='Input files to use.',nargs='+')
    parser.add_argument('--event', type=int, help='Event to look for.')
    args = parser.parse_args()
    return args

def main():
    args =commandline()

    for path in tqdm.tqdm(args.files):
        f = uproot.open(path)

        events = f['Events'].array('event')
        if np.any(events==args.event):
            print(f'Found {args.event} {path}')
            break
        # if any(args.event == ):
    
    

if __name__ == "__main__":
    main()