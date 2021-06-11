import json  
import numpy as np 
import os 

adr = os.path.dirname(os.path.abspath(__file__))
    
def ilt(fun, T, maxFnEvals, method="cme"): 
    ''' Inverse Lapalce transform'''
    if method=="cme":
        if "cmeParams" not in globals():
            with open(os.path.join(adr, "iltcme.json")) as f:
                globals()["cmeParams"] = json.load(f)
        # find the most steep CME satisfying maxFnEvals
        params = cmeParams[0]
        for p in cmeParams:
            if p["cv2"] < params["cv2"] and p["n"]+1 <= maxFnEvals:
                params = p
        eta = np.concatenate(([params["c"]], np.array(params["a"]) + 1j*np.array(params["b"])))*params["mu1"]
        beta = np.concatenate(([1], 1 + 1j*np.arange(1,params["n"]+1)*params["omega"]))*params["mu1"]
        x = T
        res = eta.dot([fun(b/x) for b in beta]).real/x
    return res  
