
import sys,os

# os.environ["TQDM_ASCII"] = "True"

from tqdm import tqdm

def progress(desc,iterable):
    return tqdm(iterable,desc=desc,ascii=True,ncols=80,smoothing=0.1)

import datetime
import time
class elapsed(object):
    def __init__(self,description):
        self.starttime = time.time()
        print(description+"...",end="",file=sys.stderr)
        sys.stderr.flush()

    def start(self):
        self.starttime = time.time()
    
    def finish(self):
        elapsed = time.time()-self.starttime
        asstr = str(datetime.timedelta(seconds=elapsed)).split('.')[0] #no milliseconds
        if asstr.startswith("0:"):
            asstr = asstr[2:]
        if asstr == "00:00":
            print(" done.",file=sys.stderr)
        else:
            print(" done. (%s elapsed)"%(asstr,),file=sys.stderr)
        sys.stderr.flush()

    def __enter__(self):
        self.start()
    
    def __exit__(self,*args):
        self.finish()

    def done(self):
        self.finish()
