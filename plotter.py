#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = 'data.out'

df = pd.read_csv(filename,
                 names=['x', 'FEM', 'analytical'],
                 delim_whitespace=True)

df.sort_values(by='x', inplace=True)

df.plot('x', ['FEM', 'analytical'], marker='o')
# df.plot('x', 'analytical', marker='s')
plt.show()
