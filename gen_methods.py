import subprocess
import sys

def run_bash_cmd(cmd_list):

    try:
        subprocess.check_output(" ".join(cmd_list), shell = True, stderr = subprocess.STDOUT, executable = "/bin/bash")

    except subprocess.CalledProcessError as e:

        print e.output
        sys.exit("bash cmd failed in run_bash_cmd. Exiting . . . ")
