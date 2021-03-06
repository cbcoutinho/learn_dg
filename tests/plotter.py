#!/usr/bin/python3

import os, sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("seaborn")

filename = "data.out"

df = pd.read_csv(filename, names=["x", "FEM"], delim_whitespace=True)

df.sort_values(by="x", inplace=True)
df.reset_index(drop=True, inplace=True)

# df.plot('x', ['FEM', 'analytical'], marker='o')
# df.plot('x', 'analytical', marker='s')

vel = -5.0
diff = 0.1
r = vel / diff
x = np.linspace(0, 1, 100)


def analytical(x):
    return (1.0 - np.exp(r * x)) / (1.0 - np.exp(r))


error = df.FEM - analytical(df.x)
# print(df.x, df.FEM, analytical(df.x))

plt.figure(1)
plt.subplot(211)
plt.plot(df.x, df.FEM, "o-", x, analytical(x), "-")
plt.xlim([np.min(x), np.max(x)])

width = df.x.diff()[1] * 0.75

plt.subplot(212)
plt.bar(df.x, error, width=width, align="center")
plt.xlim([np.min(x), np.max(x)])
# print(np.linalg.norm(error))

# print(df.x.diff()[1])

plt.show()
