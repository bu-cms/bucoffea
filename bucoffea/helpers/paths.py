
import os
import bucoffea
import subprocess
import re

pjoin = os.path.join

def bucoffea_path(path_in_repo):
    return pjoin(bucoffea.__path__[0], path_in_repo)



def vo_proxy_path():
    """Finds the path where the VO proxy file is stored."""
    cmd = ["voms-proxy-info"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    stdout, _ = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("Could not run voms-proxy-info command.")

    # Parse path from output
    pathline =  [x for x in stdout.decode("utf-8").split("\n") if x.startswith("path")][0]
    vo_path = re.sub(".* /","/", pathline)

    return vo_path