import argparse  
import numpy as np  
from geondt import *
import matplotlib.pyplot as plt  

 
if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description='Test for GeoNDT')
    # Mandatory argument 
    parser.add_argument('file', metavar='f', type=str, help='The config file to retrieve parameters from')
    args = parser.parse_args()
    print(f"Config file should be located to: {args.file}")
    with open(args.file, "r") as f:
        para = json.load(f)  

    #BE one phase
    BE = one_phase_dynamic(**para["glossary"]["GlossDiv"]["GlossList"]["GlossEntry"])   
    yt = BE.run_f()   
    plt.plot(yt) 
    plt.xlabel('Time (ms)') 
    plt.ylabel('Displacement)') 
    plt.savefig('test.png')  