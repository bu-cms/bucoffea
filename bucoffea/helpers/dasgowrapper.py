#!/usr/bin/env python

import subprocess


def das_go_query(query, json=False):
    cmd = ["dasgoclient", "--query", query]
    if json:
        cmd.append("-json")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    stdout, _ = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("Could not run DAS Go client query: {}.".format(query))

    return stdout




def das_go_query_json(query):
    proc = subprocess.Popen(
        ["dasgoclient", "--query", query, "-json"],
        stdout=subprocess.PIPE)
    stdout, _ = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("Could not run DAS Go client query: {}.".format(query))
    return json.loads(stdout)
    # rawlines =  stdout.splitlines()
    # lines = list(map(lambda x: x.decode("utf-8").strip(), rawlines))
    # lines = [l for l in lines if l]
    # print lines
    
