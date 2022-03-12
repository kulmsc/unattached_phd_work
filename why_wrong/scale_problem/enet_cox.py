from sksurv.linear_model import CoxnetSurvivalAnalysis
from sklearn.model_selection import KFold
from sksurv.metrics import concordance_index_censored
import numpy as np
import pandas as pd
import pdb

y = pd.read_csv("y_data")
x = pd.read_csv("x_data")

x = x.to_numpy()
y = y.to_numpy()

kf = KFold(n_splits=5,random_state=1,shuffle=True).split(y[:,1])

all_conc = []
all_coefs = []
n_fold = 0

for train_index, test_index in kf:
  print(n_fold)

  pdb.set_trace()
  train_y = y[train_index,:]
  train_x = x[train_index,:]

  test_y = y[test_index,:]
  test_x = x[test_index,:]

  str_train_y = np.core.records.fromarrays(train_y[:,[1,0]].transpose(), names='time, status', formats = '?, i8')
  str_test_y = np.core.records.fromarrays(test_y[:,[1,0]].transpose(), names='time, status', formats = '?, i8')

  estimator = CoxnetSurvivalAnalysis( normalize=False, # Do not normalize, the data is already normalized in our case
            l1_ratio=0.5, # Ratio between L1 and L2 regression penalties
            verbose=True,
            fit_baseline_model=True, # Fit the baselines
            copy_X=False)


  pdb.set_trace()
  pd.DataFrame(train_x).to_csv("train_x.csv")
  pd.DataFrame(train_y).to_csv("train_y.csv")

  estimator.fit(train_x, str_train_y)

  curr_conc = np.zeros((len(estimator.alphas_),2))
  curr_conc[:,0] = estimator.alphas_

  for i in range(len(estimator.alphas_)):

    pred = estimator.predict(test_x, alpha = estimator.alphas_[i])

    curr_conc[i,1] = concordance_index_censored(test_y[:,1].astype("bool"), test_y[:,0], pred)[0]


  all_conc.append(curr_conc)
  all_coefs.append(estimator.coef_)

  n_fold += 1


all_conc = sum(all_conc)/len(all_conc)
best_alpha = all_conc[all_conc[:,1] == max(all_conc[:,1]),0][0]

best_ind = np.where(all_conc[:,1] == max(all_conc[:,1]))[0][0]
best_coefs = np.zeros((all_coefs[0].shape[0]))

for i in range(len(all_coefs)):
  best_coefs += all_coefs[i][:,best_ind]
best_coefs = best_coefs/len(all_coefs)

pd.DataFrame(best_coefs).to_csv("best_coef.csv")


#estimator = CoxnetSurvivalAnalysis( normalize=False, l1_ratio=0.5, verbose=True, fit_baseline_model=True, copy_X=False, alphas = [best_alpha])

#str_y = np.core.records.fromarrays(y[:,[1,0]].transpose(), names='time, status', formats = '?, i8')
#estimator.fit(x, str_y)



#should set up a script to do cv
#just use the prediction function since the baseline hazard or survival will be the same for everyone
#then calculate concordance for each alpha
print("done")
