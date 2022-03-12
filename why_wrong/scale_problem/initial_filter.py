import numpy as np
import pandas as pd
import pdb

raw = pd.read_csv("raw_big_data.txt.gz", sep="\t", header=0)
colnames = raw.columns.to_numpy()
raw = raw.to_numpy()

count_na = []
for i in range(raw.shape[1]):
  count_na.append(sum(np.isnan(raw[:,i])))


count_low = []
for i in range(raw.shape[1]):
  vals = raw[np.logical_not(np.isnan(raw[:,i])),i]
  u_len = len(np.unique(vals))
  if u_len == 2:
    count_bin = []
    u_bin = np.unique(vals)
    for j in range(2):
      count_bin.append(sum(vals == u_bin[j]))
    count_low.append(min(count_bin))
  elif u_len < 2:
    count_low.append(0)
  else:
    count_low.append(100000)



pdb.set_trace()

print("done")
