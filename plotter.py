#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = 'data.out'

df = pd.read_csv(filename,
                 index_col=0,
                 delim_whitespace=True)
                #  names=['x', 'y1', 'y2', 'y3', 'y4', 'y5'],

df.plot()
plt.show()
