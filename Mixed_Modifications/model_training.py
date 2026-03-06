import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, auc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.model_selection import GridSearchCV
from sklearn import svm
from sklearn.inspection import permutation_importance
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import time
from sklearn.decomposition import PCA
from matplotlib.colors import LinearSegmentedColormap
from sklearn.model_selection import StratifiedKFold
from sklearn.utils import resample
import shap
from xgboost import XGBClassifier

sns.set_style("white")

plt.rcParams.update({
    'font.size': 18,
    'font.family': 'Nimbus Roman',
    'axes.titlesize': 18,
    'axes.titleweight': 'bold',
    'axes.labelsize': 20,
    'axes.labelcolor': 'black',
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'xtick.color': 'black',
    'ytick.color': 'black',
    'legend.fontsize': 14,
    'legend.title_fontsize': 16,
    'legend.loc': 'upper center',  
    'axes.edgecolor': 'black',
    'ytick.left': False  
})

data = pd.read_csv('signal_analysis/training_input.tsv', delimiter = '\t')
data_withRinfo = data

conditions = [
    data_withRinfo['contig'] == "dA",
    data_withRinfo['contig'] == "dC",
    data_withRinfo['contig'] == "dG",
    data_withRinfo['contig'] == "dT",
    data_withRinfo['contig'] == "N6mA",
    data_withRinfo['contig'] == "5hmC",
    data_withRinfo['contig'] == "8oxoG",
    data_withRinfo['contig'] == "dI",
    data_withRinfo['contig'] == "BrdU",
    data_withRinfo['contig'] == "dU",
    data_withRinfo['contig'] == "didribose",
]

choices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
choices_general = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]
# Assign to new column
data_withRinfo['Type_specific'] = np.select(conditions, choices, default=np.nan)
data_withRinfo['Type_general'] = np.select(conditions, choices_general, default=np.nan)

data_N6mA = data_withRinfo[(data_withRinfo['contig'] =='N6mA') | (data_withRinfo['contig'] =="dA")]
data_5hmC = data_withRinfo[(data_withRinfo['contig'] =='5hmC') | (data_withRinfo['contig'] =="dC")]
data_8oxoG = data_withRinfo[(data_withRinfo['contig'] =='8oxoG') | (data_withRinfo['contig'] =="dG")]
data_dI = data_withRinfo[(data_withRinfo['contig'] =='dI') | (data_withRinfo['contig'] =="dG")]
data_BrdU = data_withRinfo[(data_withRinfo['contig'] =='BrdU') | (data_withRinfo['contig'] =="dT")]
data_dU = data_withRinfo[(data_withRinfo['contig'] =='dU') | (data_withRinfo['contig'] =="dT")]
data_abasic = data_withRinfo[(data_withRinfo['contig'] =='didribose') | (data_withRinfo['contig'] =="dT")]
data_N6mA = data_withRinfo[(data_withRinfo['contig'] =='N6mA') | (data_withRinfo['Type_general'] ==0)]
data_5hmC = data_withRinfo[(data_withRinfo['contig'] =='5hmC') | (data_withRinfo['Type_general'] ==0)]
data_8oxoG = data_withRinfo[(data_withRinfo['contig'] =='8oxoG') | (data_withRinfo['Type_general'] ==0)]
data_dI = data_withRinfo[(data_withRinfo['contig'] =='dI') | (data_withRinfo['Type_general'] ==0)]
data_BrdU = data_withRinfo[(data_withRinfo['contig'] =='BrdU') | (data_withRinfo['Type_general'] ==0)]
data_dU = data_withRinfo[(data_withRinfo['contig'] =='dU') | (data_withRinfo['Type_general'] ==0)]
data_abasic = data_withRinfo[(data_withRinfo['contig'] =='didribose') | (data_withRinfo['Type_general'] ==0)]

columns_to_drop = ['Type_specific', 'Type_general', 'read_name', 'contig']
X_N6mA = data_N6mA.drop(columns = columns_to_drop)
X_N6mA.columns = [str(col) for col in X_N6mA.columns]

## Putting response variable to y
y_N6mA = data_N6mA['Type_general']

## Train-test-split
X_train_N6mA, X_test_N6mA, y_train_N6mA, y_test_N6mA = train_test_split(X_N6mA, y_N6mA, train_size = 0.8, random_state = 42, stratify=y_N6mA)
X_train_N6mA.shape, X_test_N6mA.shape
# Combine X and y for resampling
df = pd.concat([X_train_N6mA, y_train_N6mA], axis=1)

# Find smallest class size
min_size = df['Type_general'].value_counts().min()

best_xgb_N6mA = XGBClassifier(
    random_state=30,        
    n_jobs=-1,              
    max_depth=20,           
    n_estimators=1000,     
    min_child_weight=10,    
    subsample=0.8,         
    colsample_bytree=0.8,   
    use_label_encoder=False, 
    eval_metric="logloss"    
)
best_xgb_5hmC = XGBClassifier(
    random_state=30,        
    n_jobs=-1,              
    max_depth=20,           
    n_estimators=1000,     
    min_child_weight=10,    
    subsample=0.8,         
    colsample_bytree=0.8,   
    use_label_encoder=False, 
    eval_metric="logloss" 
)   
best_xgb_8oxoG = XGBClassifier(
    random_state=30,        
    n_jobs=-1,              
    max_depth=20,           
    n_estimators=1000,     
    min_child_weight=10,    
    subsample=0.8,         
    colsample_bytree=0.8,   
    use_label_encoder=False, 
    eval_metric="logloss"    
)
best_xgb_dI = XGBClassifier(
    random_state=30,        
    n_jobs=-1,              
    max_depth=20,           
    n_estimators=1000,     
    min_child_weight=10,    
    subsample=0.8,         
    colsample_bytree=0.8,   
    use_label_encoder=False, 
    eval_metric="logloss"    
)
best_xgb_BrdU = XGBClassifier(
    random_state=30,        
    n_jobs=-1,              
    max_depth=20,           
    n_estimators=1000,     
    min_child_weight=10,    
    subsample=0.8,         
    colsample_bytree=0.8,   
    use_label_encoder=False, 
    eval_metric="logloss"    
)
best_xgb_dU = XGBClassifier(
    random_state=30,        
    n_jobs=-1,              
    max_depth=20,           
    n_estimators=1000,     
    min_child_weight=10,    
    subsample=0.8,         
    colsample_bytree=0.8,   
    use_label_encoder=False, 
    eval_metric="logloss"    
)
best_xgb_abasic = XGBClassifier(
    random_state=30,        
    n_jobs=-1,              
    max_depth=20,           
    n_estimators=1000,     
    min_child_weight=10,    
    subsample=0.8,         
    colsample_bytree=0.8,   
    use_label_encoder=False, 
    eval_metric="logloss"    
)
best_xgb_N6mA.fit(X_train_N6mA, y_train_N6mA)
y_score_N6mA = best_xgb_N6mA.predict_proba(X_test_N6mA)[:, 1]

# Compute ROC curve and ROC AUC
fpr_N6mA, tpr_N6mA, _ = roc_curve(y_test_N6mA, y_score_N6mA)
roc_auc_N6mA = auc(fpr_N6mA, tpr_N6mA)

columns_to_drop = ['Type_specific', 'Type_general', 'read_name', 'contig']
X_5hmC = data_5hmC.drop(columns = columns_to_drop)
X_5hmC.columns = [str(col) for col in X_5hmC.columns]
y_5hmC = data_5hmC['Type_general']
X_train_5hmC, X_test_5hmC, y_train_5hmC, y_test_5hmC = train_test_split(X_5hmC, y_5hmC, train_size = 0.8, random_state = 42, stratify=y_5hmC)
best_xgb_5hmC.fit(X_train_5hmC, y_train_5hmC)
y_score_5hmC = best_xgb_5hmC.predict_proba(X_test_5hmC)[:, 1]
fpr_5hmC, tpr_5hmC, _ = roc_curve(y_test_5hmC, y_score_5hmC)
roc_auc_5hmC = auc(fpr_5hmC, tpr_5hmC)
columns_to_drop = ['Type_specific', 'Type_general', 'read_name', 'contig']
X_8oxoG = data_8oxoG.drop(columns = columns_to_drop)
X_8oxoG.columns = [str(col) for col in X_8oxoG.columns]
y_8oxoG = data_8oxoG['Type_general']
X_train_8oxoG, X_test_8oxoG, y_train_8oxoG, y_test_8oxoG = train_test_split(X_8oxoG, y_8oxoG, train_size = 0.8, random_state = 42, stratify=y_8oxoG)
best_xgb_8oxoG.fit(X_train_8oxoG, y_train_8oxoG)
y_score_8oxoG = best_xgb_8oxoG.predict_proba(X_test_8oxoG)[:, 1]
fpr_8oxoG, tpr_8oxoG, _ = roc_curve(y_test_8oxoG, y_score_8oxoG)
roc_auc_8oxoG = auc(fpr_8oxoG, tpr_8oxoG)
columns_to_drop = ['Type_specific', 'Type_general', 'read_name', 'contig']
X_dI = data_dI.drop(columns = columns_to_drop)
X_dI.columns = [str(col) for col in X_dI.columns]
y_dI = data_dI['Type_general']
X_train_dI, X_test_dI, y_train_dI, y_test_dI = train_test_split(X_dI, y_dI, train_size = 0.8, random_state = 42, stratify=y_dI)
best_xgb_dI.fit(X_train_dI, y_train_dI)
y_score_dI = best_xgb_dI.predict_proba(X_test_dI)[:, 1]
fpr_dI, tpr_dI, _ = roc_curve(y_test_dI, y_score_dI)
roc_auc_dI = auc(fpr_dI, tpr_dI)
columns_to_drop = ['Type_specific', 'Type_general', 'read_name', 'contig']
X_BrdU = data_BrdU.drop(columns = columns_to_drop)
X_BrdU.columns = [str(col) for col in X_BrdU.columns]
y_BrdU = data_BrdU['Type_general']
X_train_BrdU, X_test_BrdU, y_train_BrdU, y_test_BrdU = train_test_split(X_BrdU, y_BrdU, train_size = 0.8, random_state = 42, stratify=y_BrdU)
best_xgb_BrdU.fit(X_train_BrdU, y_train_BrdU)
y_score_BrdU = best_xgb_BrdU.predict_proba(X_test_BrdU)[:, 1]
fpr_BrdU, tpr_BrdU, _ = roc_curve(y_test_BrdU, y_score_BrdU)
roc_auc_BrdU = auc(fpr_BrdU, tpr_BrdU)
columns_to_drop = ['Type_specific', 'Type_general', 'read_name', 'contig']
X_dU = data_dU.drop(columns = columns_to_drop)
X_dU.columns = [str(col) for col in X_dU.columns]
y_dU = data_dU['Type_general']
X_train_dU, X_test_dU, y_train_dU, y_test_dU = train_test_split(X_dU, y_dU, train_size = 0.8, random_state = 42, stratify=y_dU)
best_xgb_dU.fit(X_train_dU, y_train_dU)
y_score_dU = best_xgb_dU.predict_proba(X_test_dU)[:, 1]
fpr_dU, tpr_dU, _ = roc_curve(y_test_dU, y_score_dU)
roc_auc_dU = auc(fpr_dU, tpr_dU)
columns_to_drop = ['Type_specific', 'Type_general', 'read_name', 'contig']
X_abasic = data_abasic.drop(columns = columns_to_drop)
X_abasic.columns = [str(col) for col in X_abasic.columns]
y_abasic = data_abasic['Type_general']
X_train_abasic, X_test_abasic, y_train_abasic, y_test_abasic = train_test_split(X_abasic, y_abasic, train_size = 0.8, random_state = 42, stratify=y_abasic)
best_xgb_abasic.fit(X_train_abasic, y_train_abasic)
y_score_abasic = best_xgb_abasic.predict_proba(X_test_abasic)[:, 1]
fpr_abasic, tpr_abasic, _ = roc_curve(y_test_abasic, y_score_abasic)
roc_auc_abasic = auc(fpr_abasic, tpr_abasic)

plt.figure(figsize=(5, 5))
plt.plot(fpr_N6mA, tpr_N6mA, linewidth=2, label=f"N6mA (AUC = {roc_auc_N6mA:.2f})", color="#d47277")
plt.plot(fpr_5hmC, tpr_5hmC, linewidth=2, label=f"5hmC (AUC = {roc_auc_5hmC:.2f})", color="#698a9c")
plt.plot(fpr_8oxoG, tpr_8oxoG, linewidth=2, label=f"8oxoG (AUC = {roc_auc_8oxoG:.2f})", color="#ede5dc")
plt.plot(fpr_dI, tpr_dI, linewidth=2, label=f"dI (AUC = {roc_auc_dI:.2f})", color="#edd8af")
plt.plot(fpr_BrdU, tpr_BrdU, linewidth=2, label=f"BrdU (AUC = {roc_auc_BrdU:.2f})", color="#a6a6a6")
plt.plot(fpr_dU, tpr_dU, linewidth=2, label=f"dU (AUC = {roc_auc_dU:.2f})", color="#8a9cbe")
plt.plot(fpr_abasic, tpr_abasic, linewidth=2, label=f"abasic (AUC = {roc_auc_abasic:.2f})", color="#9abe8a")

plt.plot([0, 1], [0, 1], 'k--')  # Diagonal line
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="lower right", frameon=False)
plt.savefig("./ROC_allmods.pdf")

