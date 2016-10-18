#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = 'data.out'

df = pd.read_csv(filename,
                 index_col=0,
                 delim_whitespace=True)

df.plot()
plt.show()
