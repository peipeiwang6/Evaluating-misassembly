'''
input 1: your input file used for machine learning models
input 2: your classes, comma-delimited
'''
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.datasets import load_iris
from sklearn.feature_selection import SelectFromModel
import pandas as pd
import numpy as np
import random
import sys
from sklearn.ensemble import RandomForestClassifier
from math import sqrt

file = sys.argv[1]
Classes = sys.argv[2]

df = pd.read_csv(file, sep='\t', index_col = 0, header = 0)
df = df[df['Class'].isin(Classes.split(','))]
X = df.iloc[:,1:]
y = df.Class
clf = RandomForestClassifier(criterion='entropy', max_features= round(sqrt(len(df.columns)-1)), n_estimators=500, n_jobs=1)
clf = clf.fit(X, y)
importances = clf.feature_importances_ 
title = df.columns.tolist()[1:]
indices = np.argsort(importances)[::-1]
temp_imp = pd.DataFrame(importances, columns = ["imp"], index=title) 
imp = pd.DataFrame(temp_imp)
imp.to_csv(file + '_Tree-based_imp', index=True, header=True,sep="\t")
