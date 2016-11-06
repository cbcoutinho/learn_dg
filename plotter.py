#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = 'data.out'

df = pd.read_csv(filename,
                 names=['x', 'FEM'],
                 delim_whitespace=True)

df.sort_values(by='x', inplace=True)

#df.plot('x', ['FEM', 'analytical'], marker='o')
# df.plot('x', 'analytical', marker='s')

r = -5.0/0.1
x = np.linspace(0, 1, 100)

def analytical(x):
    return ( 1.0 - np.exp( r * x )) / ( 1.0 - np.exp( r ))

error = analytical(df.x)-df.FEM

plt.figure(1)
plt.subplot(211)
plt.plot(df.x, df.FEM, 'o-', x, analytical(x), '-')

plt.subplot(212)
plt.bar(df.x, error, width=0.01, align='center')

print(np.linalg.norm(error))



plt.show()
