#############################################################
# Script for single cell analysis workflow
# P. RIVAUD
# 2020/04
#############################################################

import pandas as pd

def sle(filename):
    df = pd.read_csv(filename, header=0)
    return df.library_id.values.tolist()

