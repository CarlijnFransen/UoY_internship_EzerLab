import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.metrics import mean_squared_error,r2_score
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn import model_selection
import matplotlib.pyplot as plt
from sklearn.linear_model import RidgeCV, LassoCV,ElasticNetCV
from sklearn.datasets import load_diabetes
from matplotlib import pyplot as plt 
import statistics

def z_score(mean, sd, value):
    return((value - mean)/sd)



tpm_df = pd.read_csv("elastic_net_test/data/tnz_ws2_ril_TPM_table_only_control_2.csv", sep= ',',index_col = 0)

zscore_tpm_df = pd.read_csv("elastic_net_test/data/zscore_tpm_table_only_control.csv", sep= ',',index_col = 0,)
#print(zscore_tpm_df)

#tpm_df = tpm_df.dropna()
#tpm_df = tpm_df.to_numpy()
#tpm_df = tpm_df.reshape(tpm_df.shape[1:])
tpm_df = tpm_df.transpose()
print(tpm_df)
#print(tpm_df.shape)
##dms = pd.get_dummies(df[['League', 'Division', 'NewLeague']])

y_df = pd.read_csv("elastic_net_test\data\phenotypes.csv", sep=';',index_col = 0,)
##y_df = y_df.transpose()
y_df = y_df.squeeze()
#print(y_df.shape)
##X = df.drop(['Salary', 'League', 'Division', 'NewLeague'], axis=1).astype('float64')
##X = pd.concat([X_, dms[['League_N', 'Division_W', 'NewLeague_N']]], axis=1)

X_train, X_test, y_train, y_test = train_test_split(tpm_df, 
                                                    y_df, 
                                                    test_size=0.20, 
                                                    random_state=42)


enet_model = ElasticNet().fit(tpm_df, y_df)
y_predicted = enet_model.predict(tpm_df)

enet_model_2 = ElasticNet().fit(X_train,y_train)
y_test_pred =enet_model_2 .predict(X_test)
MSE_40_60_split = np.sqrt(mean_squared_error(y_test, y_test_pred ))
r2_40_60_split = r2_score(y_test, y_test_pred)
print("MSE ---> " + str(MSE_40_60_split))
print("r^2 ---> " + str(r2_40_60_split))
#no split
plt.plot(y_df, y_predicted, 'ob')
plt.show()


X_train, X_test, y_train, y_test = train_test_split(zscore_tpm_df, 
                                                    y_df, 
                                                    test_size=0.20, 
                                                    random_state=42)

enet_model = ElasticNet().fit(zscore_tpm_df, y_df)
y_predicted = enet_model.predict(zscore_tpm_df)

enet_model_2 = ElasticNet().fit(X_train,y_train)
y_test_pred =enet_model_2 .predict(X_test)
MSE_40_60_split = np.sqrt(mean_squared_error(y_test, y_test_pred ))
r2_40_60_split = r2_score(y_test, y_test_pred)
print("MSE ---> " + str(MSE_40_60_split))
print("r^2 ---> " + str(r2_40_60_split))
#no split
plt.plot(y_df, y_predicted, 'ob')
plt.show()
