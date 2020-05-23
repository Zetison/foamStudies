import pandas as pd
import itertools
import os

dir = os.getcwd()
name = "postProcessing/forceCoeffs/0/coefficient.dat"

with open(name) as handle:
    *_comments, names = itertools.takewhile(lambda line: line.startswith('#'), handle)
    # This is not the most robust way, adjust for your needs :)
    names = names[1:].split()

headers = pd.read_table(name, header=0, names=names, sep='\t', comment='#')
headers.to_csv(name[:-4]+'.csv', index=False)

