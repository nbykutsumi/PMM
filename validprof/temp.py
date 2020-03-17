# %%
%matplotlib inline
import numpy as np
from numpy import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = sns.load_dataset('tips')
sns.violinplot(x='day', y='total_bill', hue='smoker', data=df, palette='Pastel1')




# %%
