import os
import sys
import datetime
import pathlib

def init_result_dir():
    is_already_exist = True
    if not os.path.isdir("./data"):
        os.mkdir("./data")
        is_already_exist = False
    
    if not os.path.isdir("./data/latest"):
        os.mkdir("./data/latest")
        is_already_exist = False

    if is_already_exist:
        oldname = "./data/latest"
        ts = pathlib.Path(oldname).stat().st_mtime
        dt = datetime.datetime.fromtimestamp(ts)
        newname = "./data/" + dt.strftime("%Y-%m-%d_%H-%M-%S")
        os.rename(oldname, newname)
        os.mkdir("./data/latest")
        
    os.mkdir("./data/latest/wf_real")
    os.mkdir("./data/latest/flux_real")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "init":
            init_result_dir()
    
    exit()