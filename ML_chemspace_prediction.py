# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:57:42 2020

@author: Dani
"""

import pandas as pd
import numpy as np

# to plot things if needed
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score, cross_val_predict, train_test_split

# load file
biolog = pd.read_csv('D:\MRC_Postdoc\Pangenomic\Chem_space\ML_training_set.csv')

# separate variables into labels and matrix

y = biolog['Direction']
X = biolog.iloc[:, 1:85].to_numpy()

# create training and test sets for both X and y
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.05, random_state=0)

from sklearn.preprocessing import scale

X_train_scale = scale(X_train)







# create classifier and predictor
from sklearn.neighbors import KNeighborsClassifier

# initiate the classifier
knn_clf = KNeighborsClassifier()
knn_clf.fit(X_train, y_train)

# is the classifier classifying correctly?
knn_clf.predict([X_train[0]])

y_train[0]

cross_val_score(knn_clf, X_train, y_train, cv=3, scoring="accuracy", n_jobs=-1, verbose=3)


# let's do a grid search to fine tune the parameters from KNeighbourClassifier
from sklearn.model_selection import GridSearchCV
param_grid = [
    {'weights': ['uniform', 'distance'], 
     'n_neighbors': np.arange(1,80,1)}
    ]

# this way we are preparing a new classifier that will calculate accuracy from 
# the grid of parameters
grid_search = GridSearchCV(knn_clf, param_grid, cv=3,
                           scoring='accuracy',
                           return_train_score=True,
                           n_jobs=-1,
                           verbose=3)
# fit the classifier with the parameter grid
grid_search.fit(X_train, y_train)


grid_search.best_params_
grid_search.best_estimator_

cvres = grid_search.cv_results_
for mean_score, params in zip(cvres["mean_test_score"], cvres["params"]):
    print(np.sqrt(mean_score), params)

grid_search.best_params_

# with scaled values
# fit the classifier with the parameter grid
grid_search.fit(X_train_scale, y_train)


grid_search.best_params_
grid_search.best_estimator_

cvres = grid_search.cv_results_
for mean_score, params in zip(cvres["mean_test_score"], cvres["params"]):
    print(np.sqrt(mean_score), params)












#################
### RANDOM FOREST

from sklearn.ensemble import RandomForestClassifier

# initiate the classifier
rf_clf = RandomForestClassifier()
rf_clf.fit(X_train, y_train)

cross_val_score(rf_clf, X_train, y_train, cv=3, scoring="accuracy", n_jobs=-1, verbose=3)


param_grid = [
    {'criterion': ['gini', 'entropy'], 
     'n_estimators': np.arange(50,1000,50),
    'max_features': ['auto', 'sqrt', 'log2']}
    ]

# this way we are preparing a new classifier that will calculate accuracy from 
# the grid of parameters
grid_search = GridSearchCV(rf_clf, param_grid, cv=3,
                           scoring='accuracy',
                           return_train_score=True,
                           n_jobs=-1,
                           verbose=3)
# fit the classifier with the parameter grid
grid_search.fit(X_train, y_train)


grid_search.best_params_
grid_search.best_estimator_

cvres = grid_search.cv_results_
for mean_score, params in zip(cvres["mean_test_score"], cvres["params"]):
    print(np.sqrt(mean_score), params)

# roc curve
from sklearn.metrics import roc_auc_score, roc_curve

y_scores = rf_clf.fit(X_train, y_train).decision_function(X_test)


y_prob = rf_clf.predict_proba(X_test)

roc_auc_score(y_test, y_prob, multi_class="ovo",
                                  average="macro")

### RANDOM FOREST WITH SCALED VALUES



# initiate the classifier
rf_clf = RandomForestClassifier()
rf_clf.fit(X_train_scale, y_train)

cross_val_score(rf_clf, X_train_scale, y_train, cv=3, scoring="accuracy", n_jobs=-1, verbose=3)


param_grid = [
    {'criterion': ['gini', 'entropy'], 'n_estimators': [50, 100, 150, 200, 250, 300, 400, 500]}
    ]

# this way we are preparing a new classifier that will calculate accuracy from 
# the grid of parameters
grid_search = GridSearchCV(rf_clf, param_grid, cv=3,
                           scoring='accuracy',
                           return_train_score=True,
                           n_jobs=-1,
                           verbose=3)
# fit the classifier with the parameter grid
grid_search.fit(X_train_scale, y_train)


grid_search.best_params_
grid_search.best_estimator_

cvres = grid_search.cv_results_
for mean_score, params in zip(cvres["mean_test_score"], cvres["params"]):
    print(np.sqrt(mean_score), params)


##################
### MLPClassifier
from sklearn.neural_network import MLPClassifier

# initiate the classifier
MLP_clf = MLPClassifier()
MLP_clf.fit(X_train, y_train)

cross_val_score(MLP_clf, X_train, y_train, cv=3, scoring="accuracy", n_jobs=-1, verbose=3)



param_grid = [
    {'activation': ['identity', 'logistic', 'tanh', 'relu'], 
     'solver': ['lbfgs', 'sgd', 'adam'],
     'alpha': [0.0001, 0.001, 0.01, 0.1],
     'learning_rate':['constant', 'invscaling', 'adaptive']}
    ]

# this way we are preparing a new classifier that will calculate accuracy from 
# the grid of parameters
grid_search = GridSearchCV(MLP_clf, param_grid, cv=3,
                           scoring='accuracy',
                           return_train_score=True,
                           n_jobs=-1,
                           verbose=3)
# fit the classifier with the parameter grid
grid_search.fit(X_train, y_train)


grid_search.best_params_
grid_search.best_estimator_

cvres = grid_search.cv_results_
for mean_score, params in zip(cvres["mean_test_score"], cvres["params"]):
    print(np.sqrt(mean_score), params)



### Ridge Classifier
from sklearn.linear_model import RidgeClassifier

# initiate the classifier
RC_clf = RidgeClassifier()
RC_clf.fit(X_train, y_train)

cross_val_score(RC_clf, X_train, y_train, cv=3, scoring="accuracy", n_jobs=-1, verbose=3)



param_grid = [
    {'alpha': [0.5, 1, 2, 3, 4, 5], 
     'normalize': [True, False],
     'max_iter': [100, 500, 1000],
     'solver':['auto', 'svd', 'cholesky', 'lsqr', 'sag', 'saga']}
    ]

# this way we are preparing a new classifier that will calculate accuracy from 
# the grid of parameters
grid_search = GridSearchCV(RC_clf, param_grid, cv=3,
                           scoring='accuracy',
                           return_train_score=True,
                           n_jobs=-1,
                           verbose=3)
# fit the classifier with the parameter grid
grid_search.fit(X_train, y_train)


grid_search.best_params_
grid_search.best_estimator_

cvres = grid_search.cv_results_
for mean_score, params in zip(cvres["mean_test_score"], cvres["params"]):
    print(np.sqrt(mean_score), params)





#################################################################
#################################################################
#################################################################



# load file
drugbank = pd.read_csv('D:\MRC_Postdoc\Pangenomic\Chem_space\ML_predict.csv')

# separate variables into labels and matrix

drugs = drugbank.iloc[:, 1:85].to_numpy()



knn_clf = grid_search.best_estimator_
knn_clf.fit(X_train, y_train)

knn_clf.predict([drugs[0]])

knn_clf.predict_proba([drugs[0]])

labels = []

for i in range(drugs.shape[0]):
    temp = knn_clf.predict([drugs[i]])
    labels.append(temp[0])
    probs.append(probtemp)


map(knn_clf.predict, drugs)

probs = knn_clf.predict_proba(drugs)
df = pd.DataFrame(probs, columns =['Antag_prob', 'Neutral_prob', 'Synerg_prob'])  

drugbank['Direction'] = labels
drugbank['Antag_prob'] = df['Antag_prob']
drugbank['Neutral_prob'] = df['Neutral_prob']
drugbank['Synerg_prob'] = df['Synerg_prob']

finaldf = drugbank.iloc[:,[0,85,86, 87, 88]]

finaldf.to_excel('D:\MRC_Postdoc\Pangenomic\Chem_space\ML_predicted_KNN.xlsx')

