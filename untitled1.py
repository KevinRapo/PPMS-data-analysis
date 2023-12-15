# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 14:34:58 2023

@author: kevin
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Get a list of named colors
colors = list(mcolors.XKCD_COLORS.keys())

# Sample data
data = [i for i in range(0,30)]

# Plot the data with different colors
for i, value in enumerate(data):
    color = colors[i % len(colors)]  # Use modulo to cycle through the available colors
    plt.bar(i, value, color=color, label=f'Data {i+1}')

# Set legend
#plt.legend()

# Show the plot
plt.show()
#%%

import pandas as pd
import matplotlib.pyplot as plt

# Sample DataFrame with points and colors
data = pd.DataFrame({
    'X': [1, 2, 3, 4, 5],
    'Y': [10, 8, 6, 4, 2],
    'Color': ['red', 'green', 'blue', 'purple', 'orange']
})

# Create a scatter plot with colors
plt.scatter(data['X'], data['Y'], c=data['Color'], s=100)

# Customize the plot (optional)
plt.title('Scatter Plot with Custom Colors')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')

# Show the plot
plt.show()
#%%

x = 0

try:
    y = 5/x
except ZeroDivisionError:
    print("Nulliga jagamine")
    x = 10
else:
    print("Timm jagamine")
    
print(f"Edasi? {y/x}")

#%%
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

# Sample data
x = np.linspace(0, 1e6, 100)
y = x ** 2

# Create a plot
fig, ax = plt.subplots()

# Plot the data
ax.plot(x, y)

# Define a custom tick formatter function
def scientific_notation_formatter(value, pos):
    """
    Formatter function to format tick labels in scientific notation.
    """
    # Use LaTeX-style scientific notation
    return "${:0.1e}$".format(value)

# Apply the custom formatter to the x-axis
ax.xaxis.set_major_formatter(FuncFormatter(scientific_notation_formatter))

plt.show()
#%%

import merged_v1_KR as PPMS

test = PPMS.askNewDatafilePath()

