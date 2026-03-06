import numpy as np
import pandas as pd
import sys
import dask.dataframe as dd
import joblib
from cuml.linear_model import LogisticRegression as cuLogisticRegression
from cuml.ensemble import RandomForestClassifier as cuRandomForestClassifier
from cuml.svm import SVC as cuSVC
from cuml.svm import LinearSVC
from cuml.model_selection import train_test_split as gpu_train_test_split
from cuml.metrics import accuracy_score as gpu_accuracy_score
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from xgboost import XGBClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import label_binarize
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeClassifier, LogisticRegression
from sklearn.model_selection import ParameterGrid
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, roc_curve, auc, confusion_matrix, recall_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.utils import resample
from sklearn.utils import shuffle
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from scipy.stats import sem, t
from sklearn.metrics import f1_score

import seaborn as sns
import matplotlib.pyplot as plt
import cupy as cp
import xgboost as xgb
import os
import shap
from bayes_opt import BayesianOptimization

input_data = sys.argv[1]
centralR_des = sys.argv[2]
DRD60_data = sys.argv[3]
UMI_data = sys.argv[4]
TelN_ctl = sys.argv[5]
pro_r1 = sys.argv[6]
pro_r2 = sys.argv[7]
wild_ygDNA = sys.argv[8]

for d in ["numba_cache", "cupy_cache", "xdg_cache", "tmp"]:
    os.makedirs(d, exist_ok=True)
os.environ['NUMBA_CACHE_DIR'] = './numba_cache'
os.environ['CUPY_CACHE_DIR'] = './cupy_cache'
os.environ['XDG_CACHE_HOME'] = './xdg_cache'
os.environ['TMPDIR'] = './tmp'

sns.set_style("white")

# Set font family and base size
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
    'axes.edgecolor': 'black'
})

data = pd.read_csv(input_data, delimiter = '\t')
centralR = pd.read_csv(centralR_des, delimiter = '\t')

data.columns = ['read_name', 'contig','ampl_mean_1st', 'ampl_mean_2nd', 'ampl_mean_3rd', 'ampl_mean_4th','ampl_mean_5th','ampl_mean_6th', 'ampl_mean_7th', 'ampl_mean_8th', 'ampl_mean_9th', 'model_deviation_1st', 'model_deviation_2nd', 'model_deviation_3rd', 'model_deviation_4th','model_deviation_5th','model_deviation_6th','model_deviation_7th','model_deviation_8th','model_deviation_9th','Type', 'dwell_time_R-18', 'dwell_time_R-15']
data_withRinfo = pd.merge(data, centralR, on = 'contig', how = 'left')
data_withRinfo = data_withRinfo[~data_withRinfo['read_name'].str.contains(r'(_upa|_downa|_downb|_downc)', regex=True)]

conditions = [
    data_withRinfo['Type'] == 0,
    data_withRinfo['cen_R'] == 'A',
    data_withRinfo['cen_R'] == 'C',
    data_withRinfo['cen_R'] == 'G',
    data_withRinfo['cen_R'] == 'T'
]

choices = [0, 1, 2, 3, 4]

# Assign ribo info to new column
data_withRinfo['rType'] = np.select(conditions, choices, default=np.nan)

# Features, currently excluding dwell time at all positions
features = ["ampl_mean_1st","ampl_mean_2nd","ampl_mean_3rd","ampl_mean_4th", "ampl_mean_5th", "ampl_mean_6th", "ampl_mean_7th", "ampl_mean_8th", "ampl_mean_9th",   
            "model_deviation_1st","model_deviation_2nd","model_deviation_3rd","model_deviation_4th", "model_deviation_5th", "model_deviation_6th", "model_deviation_7th", "model_deviation_8th", "model_deviation_9th",
            ]
X = data_withRinfo[features].to_numpy().astype(np.float32)
y = data_withRinfo['rType'].to_numpy().astype(np.int32)

X_train, X_test, y_train, y_test = gpu_train_test_split(X, y, test_size=0.2, random_state=42)

# PCA clustering:
imputer = SimpleImputer(strategy='median')
X_imputed = imputer.fit_transform(X)

# Run PCA on all data
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_imputed)
var_ratio = pca.explained_variance_ratio_ * 100

color_dict = {
    0:'#DB7176',  
    1:'#698A9C',
    2:'#A7C7E7',
    3:'#878580',
    4:'#C0B6A4'  
}

label_dict = {
    0:'dN',
    1:'rA',
    2:'rC',
    3:'rG',
    4:'rU'
}

# Map each label in y to its corresponding color
colors = [color_dict[label] for label in y]

# Create JointGrid
g = sns.JointGrid(height=5, ratio=2, space=0.1)

scatter = g.ax_joint.scatter(
    X_pca[:, 0], X_pca[:, 1], 
    c=colors, alpha=0.2, edgecolor=colors, rasterized = True
)
for val in sorted(set(y)):
    plt.scatter([], [], c=color_dict[val], label=label_dict[val])
plt.legend(title='Type', loc='best', fontsize=9, markerscale=0.7, title_fontsize=10, frameon=False)

# PC1 density distribution
for group in [0, 1, 2, 3, 4]:
    mask = (y == group)
    sns.kdeplot(
        x=X_pca[mask, 0],
        ax=g.ax_marg_x,
        color=color_dict[group],
        fill=True,
        alpha=0.5,
        label=label_dict[group]
    )
    

# PC2 density distribution
for group in [0, 1, 2, 3, 4]:
    mask = (y == group)
    sns.kdeplot(
        y=X_pca[mask, 1],
        ax=g.ax_marg_y,
        color=color_dict[group],
        fill=True,
        alpha=0.5,
        label=label_dict[group]
    )

g.ax_joint.set_xlabel(f'PC1 ({var_ratio[0]:.2f}%)')
g.ax_joint.set_ylabel(f'PC2 ({var_ratio[1]:.2f}%)')
plt.tight_layout()
plt.savefig("PCA_plot.pdf")

X_train, X_test, y_train, y_test = gpu_train_test_split(X_imputed, y, test_size=0.2, random_state=42)
# Try different models
model_configs = {
    "cuML_logistic": (cuLogisticRegression(), {"C": [0.01]}),
    "cuML_rf": (cuRandomForestClassifier(), {"n_estimators": [300,200], "max_depth": [20]}),
    "cuML_svc": (cuSVC(probability=True), {"C": [1]}),
    "cuML_linear_svc": (LinearSVC(), {"C": [0.1]}),
    "xgb_gpu": (XGBClassifier(tree_method='hist', device="cuda", eval_metric='mlogloss'),
                {"max_depth": [6], "n_estimators": [300, 200]}),
    "sklearn_ridge": (RidgeClassifier(), {"alpha": [0.1]}),
    "sklearn_lasso": (LogisticRegression(penalty='l1', solver='saga', max_iter=1000), {"C": [1.0]}),
    "sklearn_elastic": (LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, max_iter=1000),
                    {"C": [1.0], "l1_ratio": [0.3]})
}

# 10-fold Cross-validation setup
kf = KFold(n_splits=10, shuffle=True, random_state=42)
scores = {}
best_models = {}
cv_scores_per_model = {}

# Evaluate each model with hyperparameter tuning
for name, (model_proto, param_grid) in model_configs.items():
    best_cv_score = -np.inf
    best_params = None
    best_fold_scores = None

    for params in ParameterGrid(param_grid):
        fold_scores = []

        for train_idx, val_idx in kf.split(X_train):
            X_tr, X_val = X_train[train_idx], X_train[val_idx]
            y_tr, y_val = y_train[train_idx], y_train[val_idx]

            model = type(model_proto)(**params)  # fresh instance with current hyperparams
            model.fit(X_tr, y_tr)
            y_pred = model.predict(X_val)
            acc = gpu_accuracy_score(y_val, y_pred)
            fold_scores.append(acc)

        mean_score = np.mean(fold_scores)
        if mean_score > best_cv_score:
            best_cv_score = mean_score
            best_params = params
            best_fold_scores = fold_scores

    scores[name] = best_cv_score
    best_models[name] = type(model_proto)(**best_params)
    cv_scores_per_model[name] = best_fold_scores

desired_order = [
    "cuML_logistic",
    "sklearn_ridge", "sklearn_lasso", "sklearn_elastic", "cuML_rf","xgb_gpu", "cuML_svc",
    "cuML_linear_svc"
]
custom_labels = [
    "Logistic Regression",
    "Ridge Regression", "LASSO Regression", "Elastic Net", "Random Forest", "XGBoost", "Support Vector Classifier",
    "Linear Support Vector Classifier"
]


lower_errors = [np.mean(cv_scores_per_model[name]) - np.min(cv_scores_per_model[name]) for name in desired_order]
upper_errors = [np.max(cv_scores_per_model[name]) - np.mean(cv_scores_per_model[name]) for name in desired_order]
asymmetric_errors = [lower_errors, upper_errors]

yticks = np.arange(0.5, 1.05, 0.05)
model_names = list(cv_scores_per_model.keys())
means = [np.mean(cv_scores_per_model[name]) for name in desired_order]
stds = [np.std(cv_scores_per_model[name]) for name in desired_order]

plt.figure(figsize=(12, 6))
plt.bar(
    range(len(desired_order)), means, 
    yerr=asymmetric_errors, capsize=10, error_kw=dict(lw=2), color="#698A9C"
)
plt.yticks(yticks)
plt.ylim(0.5, 1.0)
plt.ylabel("Accuracy")
plt.title("10-fold Cross-validation Accuracy (mean ± range) per Model")
plt.xticks(ticks=range(len(desired_order)), labels=custom_labels, rotation=45, ha='right')
plt.tight_layout()
plt.savefig("model_selection.pdf", dpi=300)

# Train XGBoost with tuned parameters
X_train, X_test, y_train, y_test = gpu_train_test_split(X, y, test_size=0.2, random_state=42)

custom_color = "#DB7176"
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_blue", ["#ffffff", custom_color])

# Setup
custom_labels = {
    0: "dN",
    1: "rA",
    2: "rC",
    3: "rG",
    4: "rU"
}
label_names = [custom_labels[i] for i in range(len(custom_labels))]

classes = np.unique(y_train)
n_classes = len(classes)
best_params = {'colsample_bytree': np.float64(0.7721), 'gamma': np.float64(0.3028), 'learning_rate': np.float64(0.1768), 'max_depth': np.float64(8.81), 'min_child_weight': np.float64(5.81), 'n_estimators': np.float64(1.182e+03), 'reg_alpha': np.float64(8.457), 'reg_lambda': np.float64(6.806), 'subsample': np.float64(0.7801)}
best_model = XGBClassifier(
    max_depth=int(best_params['max_depth']),
    learning_rate=best_params['learning_rate'],
    n_estimators=int(best_params['n_estimators']),
    gamma=best_params['gamma'],
    min_child_weight=best_params['min_child_weight'],
    subsample=best_params['subsample'],
    colsample_bytree=best_params['colsample_bytree'],
    reg_alpha=best_params['reg_alpha'],
    reg_lambda=best_params['reg_lambda'],
    tree_method='hist',
    device='cuda',
    eval_metric='mlogloss',
    random_state=42
)

# Cross-validation setup
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
fpr_grid = np.linspace(0, 1, 100)
tpr_all = []  
auc_folds = []
conf_matrices = []

for train_idx, val_idx in cv.split(X_train, y_train):
    X_tr, X_val = X_train[train_idx], X_train[val_idx]
    y_tr, y_val = y_train[train_idx], y_train[val_idx]

    y_val_bin = label_binarize(y_val, classes=classes)

    model = best_model
    model.fit(X_tr, y_tr)
    y_score = model.predict_proba(X_val)

    # Loop over classes
    for i in range(n_classes):
        fpr, tpr, _ = roc_curve(y_val_bin[:, i], y_score[:, i])
        interp_tpr = np.interp(fpr_grid, fpr, tpr)
        interp_tpr[0] = 0.0
        tpr_all.append(interp_tpr)

    # Macro-average AUC for reference
    auc_macro = np.mean([
        auc(roc_curve(y_val_bin[:, i], y_score[:, i])[0],
            roc_curve(y_val_bin[:, i], y_score[:, i])[1])
        for i in range(n_classes)
    ])
    auc_folds.append(auc_macro)

    y_pred = model.predict(X_val)

    # Confusion Matrix
    cm = confusion_matrix(y_val, y_pred, labels=np.arange(n_classes))
    cm_normalized = cm.astype(float) / cm.sum(axis=1, keepdims=True)
    conf_matrices.append(cm_normalized)

# Convert and aggregate
tpr_all = np.array(tpr_all)  # shape = (n_classes × n_folds, len(fpr_grid))
mean_tpr = np.mean(tpr_all, axis=0)
std_tpr = np.std(tpr_all, axis=0, ddof=1)
se_tpr = std_tpr / np.sqrt(tpr_all.shape[0])
ci_upper = np.clip(mean_tpr + 1.96 * se_tpr, 0, 1)
ci_lower = np.clip(mean_tpr - 1.96 * se_tpr, 0, 1)
mean_auc = np.mean(auc_folds)
mean_cm = np.mean(conf_matrices, axis=0)

final_model = best_model
final_model.fit(X_train, y_train)

y_test_bin = label_binarize(y_test, classes=[0, 1, 2, 3, 4])
y_test_proba = final_model.predict_proba(X_test)
y_test_pred = final_model.predict(X_test)

# ROC per class
fpr_test = np.linspace(0, 1, 100)
tpr_test_list = []
auc_test_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_test_proba[:, i])
    interp_tpr = np.interp(fpr_test, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_test_list.append(interp_tpr)
    auc_test_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_test = np.mean(tpr_test_list, axis=0)
auc_macro_test = np.mean(auc_test_list)
cm_test = confusion_matrix(y_test, y_test_pred, labels=np.arange(n_classes))
cm_test_normalized = cm_test.astype(float) / cm_test.sum(axis=1, keepdims=True)

# DRD60 validation
data_DRD60 = pd.read_csv(DRD60_data, delimiter = '\t')
data_DRD60.columns = ['read_name', 'contig','ampl_mean_1st', 'ampl_mean_2nd', 'ampl_mean_3rd', 'ampl_mean_4th','ampl_mean_5th','ampl_mean_6th', 'ampl_mean_7th', 'ampl_mean_8th', 'ampl_mean_9th', 'model_deviation_1st', 'model_deviation_2nd', 'model_deviation_3rd', 'model_deviation_4th','model_deviation_5th','model_deviation_6th','model_deviation_7th','model_deviation_8th','model_deviation_9th','Type', 'dwell_time_R-18', 'dwell_time_R-15']
data_DRD60_withRinfo = pd.merge(data_DRD60, centralR, on = 'contig', how = 'left')

conditions = [
    data_DRD60_withRinfo['Type'] == 0,
    data_DRD60_withRinfo['cen_R'] == 'A',
    data_DRD60_withRinfo['cen_R'] == 'C',
    data_DRD60_withRinfo['cen_R'] == 'G',
    data_DRD60_withRinfo['cen_R'] == 'T'
]

data_DRD60_withRinfo['rType'] = np.select(conditions, choices, default=np.nan)
X_DRD60 = data_DRD60_withRinfo[features].to_numpy().astype(np.float32)
y_DRD60 = data_DRD60_withRinfo['rType'].to_numpy().astype(np.int32)
y_DRD60_bin = label_binarize(y_DRD60, classes=[0, 1, 2, 3, 4])
y_DRD60_proba = final_model.predict_proba(X_DRD60)
y_DRD60_pred = final_model.predict(X_DRD60)
# ROC per class
fpr_DRD60 = np.linspace(0, 1, 100)
tpr_DRD60_list = []
auc_DRD60_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_DRD60_bin[:, i], y_DRD60_proba[:, i])
    interp_tpr = np.interp(fpr_DRD60, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_DRD60_list.append(interp_tpr)
    auc_DRD60_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_DRD60 = np.mean(tpr_DRD60_list, axis=0)
auc_macro_DRD60 = np.mean(auc_DRD60_list)
cm_DRD60 = confusion_matrix(y_DRD60, y_DRD60_pred, labels=np.arange(n_classes))
cm_DRD60_normalized = cm_DRD60.astype(float) / cm_DRD60.sum(axis=1, keepdims=True)

#Mixed-SRICO validation
data_UMI = pd.read_csv(UMI_data, delimiter = '\t')
data_UMI.columns = ['read_name', 'contig','ampl_mean_1st', 'ampl_mean_2nd', 'ampl_mean_3rd', 'ampl_mean_4th','ampl_mean_5th','ampl_mean_6th', 'ampl_mean_7th', 'ampl_mean_8th', 'ampl_mean_9th', 'model_deviation_1st', 'model_deviation_2nd', 'model_deviation_3rd', 'model_deviation_4th','model_deviation_5th','model_deviation_6th','model_deviation_7th','model_deviation_8th','model_deviation_9th','Type', 'dwell_time_R-18', 'dwell_time_R-15']
data_UMI_withRinfo = pd.merge(data_UMI, centralR, on = 'contig', how = 'left')
data_teln_realctl = pd.read_csv(TelN_ctl, delimiter = '\t')
data_teln_realctl.columns = ['read_name', 'contig','ampl_mean_1st', 'ampl_mean_2nd', 'ampl_mean_3rd', 'ampl_mean_4th','ampl_mean_5th','ampl_mean_6th', 'ampl_mean_7th', 'ampl_mean_8th', 'ampl_mean_9th', 'model_deviation_1st', 'model_deviation_2nd', 'model_deviation_3rd', 'model_deviation_4th','model_deviation_5th','model_deviation_6th','model_deviation_7th','model_deviation_8th','model_deviation_9th','Type', 'dwell_time_R-18', 'dwell_time_R-15']
data_teln_realctl_withRinfo = pd.merge(data_teln_realctl, centralR, on = 'contig', how = 'left')
combined_UMI = pd.concat([data_UMI_withRinfo, data_teln_realctl_withRinfo], axis=0, ignore_index=True)

conditions = [
    combined_UMI['Type'] == 0,
    combined_UMI['cen_R'] == 'A',
    combined_UMI['cen_R'] == 'C',
    combined_UMI['cen_R'] == 'G',
    combined_UMI['cen_R'] == 'T'
]

combined_UMI['rType'] = np.select(conditions, choices, default=np.nan)
X_UMI = combined_UMI[features].to_numpy().astype(np.float32)
y_UMI = combined_UMI['rType'].to_numpy().astype(np.int32)
y_UMI_bin = label_binarize(y_UMI, classes=[0, 1, 2, 3, 4])
y_UMI_proba = final_model.predict_proba(X_UMI)
y_UMI_pred = final_model.predict(X_UMI)
# ROC per class
fpr_UMI = np.linspace(0, 1, 100)
tpr_UMI_list = []
auc_UMI_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_UMI_bin[:, i], y_UMI_proba[:, i])
    interp_tpr = np.interp(fpr_UMI, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_UMI_list.append(interp_tpr)
    auc_UMI_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_UMI = np.mean(tpr_UMI_list, axis=0)
auc_macro_UMI = np.mean(auc_UMI_list)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)

# Plot ROC
plt.figure(figsize=(7, 7))
plt.plot(fpr_grid, mean_tpr,
         label=f"10-fold CV with 95% CI (Mean AUC = {mean_auc:.2f})",
         color="#698A9C", lw=2)

plt.plot(fpr_test, mean_tpr_test,
         label=f"Hold-out Test-set (AUC = {auc_macro_test:.2f})",
         color="#C0B6A4", lw=2)

plt.plot(fpr_DRD60, mean_tpr_DRD60,
         label=f"DRD60 Test-set (AUC = {auc_macro_DRD60:.2f})",
         color="#DB7176", lw=2)

plt.plot(fpr_UMI, mean_tpr_UMI,
         label=f"Mixed SRICO Test-set (AUC = {auc_macro_UMI:.2f})",
         color="#878580", lw=2)

plt.fill_between(fpr_grid, ci_lower, ci_upper, color="#698A9C", alpha=0.4)
plt.plot([0, 1], [0, 1], 'k--', lw=1)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Macro-average ROC for Multiclass XGboost')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("ROC_XGBoost.pdf", dpi=300)

# Plot Confusion Matrix
custom_color = "#DB7176"
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_blue", ["#ffffff", custom_color])
threshold = 0.5
y_HO_pred = []

for prob in y_test_proba:
    # If no prob exceeds threshold, fallback to argmax
    if (prob >= threshold).any():
        pred_class = np.argmax(prob * (prob >= threshold))
    else:
        pred_class = np.argmax(prob)
    y_HO_pred.append(pred_class)

y_HO_pred = np.array(y_HO_pred)
cm_HO = confusion_matrix(y_test, y_HO_pred, labels=np.arange(n_classes))
cm_HO_normalized = cm_HO.astype(float) / cm_HO.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_HO_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in Hold-out Test-set")
plt.tight_layout()
plt.savefig("Holdout_Confusion.pdf", dpi=300)

y_DRD60_pred = []

for prob in y_DRD60_proba:
    # If no prob exceeds threshold, fallback to argmax
    if (prob >= threshold).any():
        pred_class = np.argmax(prob * (prob >= threshold))
    else:
        pred_class = np.argmax(prob)
    y_DRD60_pred.append(pred_class)

y_DRD60_pred = np.array(y_DRD60_pred)
cm_DRD60 = confusion_matrix(y_DRD60, y_DRD60_pred, labels=np.arange(n_classes))
cm_DRD60_normalized = cm_DRD60.astype(float) / cm_DRD60.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_DRD60_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in DRD60 Test-set")
plt.tight_layout()
plt.savefig("DRD60_Confusion.pdf", dpi=300)

y_UMI_pred = []

for prob in y_UMI_proba:
    # Classes meeting the threshold
    over_thresh = np.where(prob >= threshold)[0]
    
    if over_thresh.size > 0:
        # Pick the one with the highest probability among those over the threshold
        pred_class = over_thresh[np.argmax(prob[over_thresh])]
    else:
        # Fallback to the usual best guess
        pred_class = 0
    
    y_UMI_pred.append(pred_class)

y_UMI_pred = np.array(y_UMI_pred)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_UMI_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in Mixed SRICO Test-set")
plt.tight_layout()
plt.savefig("Mixed-SRICO_Confusion.pdf", dpi=300)

# Feature Correlation
X_train = pd.DataFrame(X_train, columns=features)
corr_matrix = X_train.corr()
custom_cmap = LinearSegmentedColormap.from_list(
    "custom_corr", ["#698A9C", "white", "#DB7176"]
)

plt.figure(figsize=(8,8))
sns.heatmap(corr_matrix, annot=False, fmt=".2f", cmap=custom_cmap, 
            square=True, linewidths=0.5, cbar_kws={"shrink": 0.8})

plt.title("Feature Correlation", fontsize=16, fontweight='bold')
plt.xticks(
    np.arange(len(features)) + 0.5,  
    features,
    rotation=45,
    ha='right',
    size=10,  
    fontweight='normal'
)

plt.yticks(
    np.arange(len(features)) + 0.5,  
    features,
    rotation=0,
    size=10,  
    fontweight='normal'
)

plt.tight_layout()
plt.savefig("Feature_Correlation.pdf", dpi=300)

# Feature Importance
custom_labels = {
    0: "dN",
    1: "rA",
    2: "rC",
    3: "rG",
    4: "rU"
}
class_names_list = [custom_labels[i] for i in range(len(custom_labels))]
custom_colors = ['#A7C7E7', '#C0B6A4', '#878580', '#698A9C', '#DB7176']

explainer = shap.TreeExplainer(final_model)
X_df = data_withRinfo[features]
X_sample = X_df.sample(n=5000, random_state=42)
#X_sample_np = X_sample.to_numpy()
shap_values = explainer.shap_values(X_sample)

shap.summary_plot(
    shap_values, X_sample,
    plot_type="bar",
    class_names=class_names_list, 
    show=False
)
# Get the axis and bars
ax = plt.gca()
bars = ax.patches

# Number of features (assumes all classes have the same number of bars)
num_features = len(bars) // len(custom_colors)

# Assign colors per class (not per bar)
for class_idx in range(len(custom_colors)):
    # Get all bars for the current class
    start_idx = class_idx * num_features
    end_idx = (class_idx + 1) * num_features
    for bar in bars[start_idx:end_idx]:
        bar.set_color(custom_colors[class_idx])

# Update legend colors to match
#legend = ax.legend_
#for legend_handle, color in zip(legend.legendHandles, custom_colors):
    #legend_handle.set_color(color)

# Optional: move legend to bottom right
plt.legend(loc='lower right', frameon=True)
plt.xlabel("Mean SHAP Value")
plt.ylabel("Feature")
plt.tight_layout()
plt.savefig("Feature_Importance.pdf", dpi=300)

# Try to also include adjacent true negatives in training

data_withRinfo = pd.merge(data, centralR, on = 'contig', how = 'left')
data_withRinfo = data_withRinfo[~data_withRinfo['read_name'].str.contains(r'(_downc)', regex=True)]

conditions = [
    data_withRinfo['Type'] == 0,
    data_withRinfo['cen_R'] == 'A',
    data_withRinfo['cen_R'] == 'C',
    data_withRinfo['cen_R'] == 'G',
    data_withRinfo['cen_R'] == 'T'
]

data_withRinfo['rType'] = np.select(conditions, choices, default=np.nan)

# Features and target
features = ["ampl_mean_2nd","ampl_mean_3rd","ampl_mean_4th", "ampl_mean_5th", "ampl_mean_6th", "ampl_mean_7th",    
            "model_deviation_2nd","model_deviation_3rd","model_deviation_4th", "model_deviation_5th", "model_deviation_6th", "model_deviation_7th", 
            ]
X = data_withRinfo[features].to_numpy().astype(np.float32)
y = data_withRinfo['rType'].to_numpy().astype(np.int32)

X_train, X_test, y_train, y_test = gpu_train_test_split(X, y, test_size=0.2, random_state=42)
X_train_df = pd.DataFrame(X_train, columns=features)
y_train_df = pd.Series(y_train, name='rType')

df = pd.concat([X_train_df, y_train_df], axis=1)

# Balance rNMP classes by down-sampling to the smallest class
min_size = df['rType'].value_counts().min()
df_balanced = (
    df.groupby('rType')
      .apply(lambda x: x.sample(min_size, random_state=42))
      .reset_index(drop=True)
)

# Split back into X and y
X_train_bal = df_balanced.drop(columns=['rType'])
y_train_bal = df_balanced['rType']

custom_color = "#DB7176"
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_blue", ["#ffffff", custom_color])

# Setup
custom_labels = {
    0: "dN",
    1: "rA",
    2: "rC",
    3: "rG",
    4: "rU"
}
label_names = [custom_labels[i] for i in range(len(custom_labels))]

classes = np.unique(y_train_bal)
n_classes = len(classes)
best_params = {'colsample_bytree': np.float64(0.7721), 'gamma': np.float64(0.3028), 'learning_rate': np.float64(0.1768), 'max_depth': np.float64(8.81), 'min_child_weight': np.float64(5.81), 'n_estimators': np.float64(1.182e+03), 'reg_alpha': np.float64(8.457), 'reg_lambda': np.float64(6.806), 'subsample': np.float64(0.7801)}
best_model = XGBClassifier(
    max_depth=int(best_params['max_depth']),
    learning_rate=best_params['learning_rate'],
    n_estimators=int(best_params['n_estimators']),
    gamma=best_params['gamma'],
    min_child_weight=best_params['min_child_weight'],
    subsample=best_params['subsample'],
    colsample_bytree=best_params['colsample_bytree'],
    reg_alpha=best_params['reg_alpha'],
    reg_lambda=best_params['reg_lambda'],
    tree_method='hist',
    device='cuda',
    eval_metric='mlogloss',
    random_state=42
)

# Cross-validation setup
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
fpr_grid = np.linspace(0, 1, 100)
tpr_all = [] 
auc_folds = []
conf_matrices = []


X_train_bal_cp = cp.asarray(X_train_bal.values)
y_train_bal_cp = cp.asarray(y_train_bal.values)
y_train_bal_np = y_train_bal.values  

for train_idx, val_idx in cv.split(X_train_bal, y_train_bal_np):
    # Use indices to slice CuPy arrays
    X_tr = cp.take(X_train_bal_cp, train_idx, axis=0)
    X_val = cp.take(X_train_bal_cp, val_idx, axis=0)
    y_tr = cp.take(y_train_bal_cp, train_idx)
    y_val = cp.take(y_train_bal_cp, val_idx)

    y_val_cpu = cp.asnumpy(y_val)
    y_val_bin = label_binarize(y_val_cpu, classes=classes)

    model = best_model
    model.fit(X_tr, y_tr)

    y_score = cp.asnumpy(model.predict_proba(X_val))

    for i in range(n_classes):
        fpr, tpr, _ = roc_curve(y_val_bin[:, i], y_score[:, i])
        interp_tpr = np.interp(fpr_grid, fpr, tpr)
        interp_tpr[0] = 0.0
        tpr_all.append(interp_tpr)

    auc_macro = np.mean([
        auc(roc_curve(y_val_bin[:, i], y_score[:, i])[0],
            roc_curve(y_val_bin[:, i], y_score[:, i])[1])
        for i in range(n_classes)
    ])
    auc_folds.append(auc_macro)

    y_pred = model.predict(X_val)

    # Confusion Matrix
    cm = confusion_matrix(y_val_cpu, y_pred, labels=np.arange(n_classes))
    cm_normalized = cm.astype(float) / cm.sum(axis=1, keepdims=True)
    conf_matrices.append(cm_normalized)


# Convert and aggregate
tpr_all = np.array(tpr_all)  
mean_tpr = np.mean(tpr_all, axis=0)
std_tpr = np.std(tpr_all, axis=0, ddof=1)
se_tpr = std_tpr / np.sqrt(tpr_all.shape[0])
ci_upper = np.clip(mean_tpr + 1.96 * se_tpr, 0, 1)
ci_lower = np.clip(mean_tpr - 1.96 * se_tpr, 0, 1)
mean_auc = np.mean(auc_folds)
mean_cm = np.mean(conf_matrices, axis=0)

final_model = best_model
final_model.fit(X_train_bal, y_train_bal)

y_test_bin = label_binarize(y_test, classes=[0, 1, 2, 3, 4])
y_test_proba = final_model.predict_proba(X_test)
y_test_pred = final_model.predict(X_test)

# ROC per class
fpr_test = np.linspace(0, 1, 100)
tpr_test_list = []
auc_test_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_test_proba[:, i])
    interp_tpr = np.interp(fpr_test, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_test_list.append(interp_tpr)
    auc_test_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_test = np.mean(tpr_test_list, axis=0)
auc_macro_test = np.mean(auc_test_list)
cm_test = confusion_matrix(y_test, y_test_pred, labels=np.arange(n_classes))
cm_test_normalized = cm_test.astype(float) / cm_test.sum(axis=1, keepdims=True)

data_DRD60_withRinfo = pd.merge(data_DRD60, centralR, on = 'contig', how = 'left')
conditions = [
    data_DRD60_withRinfo['Type'] == 0,
    data_DRD60_withRinfo['cen_R'] == 'A',
    data_DRD60_withRinfo['cen_R'] == 'C',
    data_DRD60_withRinfo['cen_R'] == 'G',
    data_DRD60_withRinfo['cen_R'] == 'T'
]

data_DRD60_withRinfo['rType'] = np.select(conditions, choices, default=np.nan)

X_DRD60 = data_DRD60_withRinfo[features].to_numpy().astype(np.float32)
y_DRD60 = data_DRD60_withRinfo['rType'].to_numpy().astype(np.int32)
y_DRD60_bin = label_binarize(y_DRD60, classes=[0, 1, 2, 3, 4])
y_DRD60_proba = final_model.predict_proba(X_DRD60)
y_DRD60_pred = final_model.predict(X_DRD60)
fpr_DRD60 = np.linspace(0, 1, 100)
tpr_DRD60_list = []
auc_DRD60_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_DRD60_bin[:, i], y_DRD60_proba[:, i])
    interp_tpr = np.interp(fpr_DRD60, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_DRD60_list.append(interp_tpr)
    auc_DRD60_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_DRD60 = np.mean(tpr_DRD60_list, axis=0)
auc_macro_DRD60 = np.mean(auc_DRD60_list)
cm_DRD60 = confusion_matrix(y_DRD60, y_DRD60_pred, labels=np.arange(n_classes))
cm_DRD60_normalized = cm_DRD60.astype(float) / cm_DRD60.sum(axis=1, keepdims=True)

data_UMI_withRinfo = pd.merge(data_UMI, centralR, on = 'contig', how = 'left')
data_teln_realctl_withRinfo = pd.merge(data_teln_realctl, centralR, on = 'contig', how = 'left')
combined_UMI = pd.concat([data_UMI_withRinfo, data_teln_realctl_withRinfo], axis=0, ignore_index=True)

conditions = [
    combined_UMI['Type'] == 0,
    combined_UMI['cen_R'] == 'A',
    combined_UMI['cen_R'] == 'C',
    combined_UMI['cen_R'] == 'G',
    combined_UMI['cen_R'] == 'T'
]

combined_UMI['rType'] = np.select(conditions, choices, default=np.nan)

X_UMI = combined_UMI[features].to_numpy().astype(np.float32)
y_UMI = combined_UMI['rType'].to_numpy().astype(np.int32)
y_UMI_bin = label_binarize(y_UMI, classes=[0, 1, 2, 3, 4])
y_UMI_proba = final_model.predict_proba(X_UMI)
y_UMI_pred = final_model.predict(X_UMI)

# ROC per class
fpr_UMI = np.linspace(0, 1, 100)
tpr_UMI_list = []
auc_UMI_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_UMI_bin[:, i], y_UMI_proba[:, i])
    interp_tpr = np.interp(fpr_UMI, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_UMI_list.append(interp_tpr)
    auc_UMI_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_UMI = np.mean(tpr_UMI_list, axis=0)
auc_macro_UMI = np.mean(auc_UMI_list)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)

# Plot
plt.figure(figsize=(7, 7))
plt.plot(fpr_grid, mean_tpr,
         label=f"10-fold CV with 95% CI (Mean AUC = {mean_auc:.2f})",
         color="#698A9C", lw=2)

# Test ROC
plt.plot(fpr_test, mean_tpr_test,
         label=f"Hold-out Test-set (AUC = {auc_macro_test:.2f})",
         color="#C0B6A4", lw=2)

plt.plot(fpr_DRD60, mean_tpr_DRD60,
         label=f"DRD60 Test-set (AUC = {auc_macro_DRD60:.2f})",
         color="#DB7176", lw=2)

plt.plot(fpr_UMI, mean_tpr_UMI,
         label=f"Mixed SRICO Test-set (AUC = {auc_macro_UMI:.2f})",
         color="#878580", lw=2)

plt.fill_between(fpr_grid, ci_lower, ci_upper, color="#698A9C", alpha=0.4)
plt.plot([0, 1], [0, 1], 'k--', lw=1)

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Macro-average ROC for Multiclass XGboost')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("ROC_XGBoost_final.pdf", dpi=300)

threshold = 0.5
y_HO_pred = []

for prob in y_test_proba:
    # If no prob exceeds threshold, fallback to argmax
    if (prob >= threshold).any():
        pred_class = np.argmax(prob * (prob >= threshold))
    else:
        pred_class = np.argmax(prob)
    y_HO_pred.append(pred_class)

y_HO_pred = np.array(y_HO_pred)
cm_HO = confusion_matrix(y_test, y_HO_pred, labels=np.arange(n_classes))
cm_HO_normalized = cm_HO.astype(float) / cm_HO.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_HO_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in Hold-out Test-set")
plt.tight_layout()
plt.savefig("Holdout_Confusion_final.pdf", dpi=300)

y_DRD60_pred = []

for prob in y_DRD60_proba:
    # If no prob exceeds threshold, fallback to argmax
    if (prob >= threshold).any():
        pred_class = np.argmax(prob * (prob >= threshold))
    else:
        pred_class = np.argmax(prob)
    y_DRD60_pred.append(pred_class)

y_DRD60_pred = np.array(y_DRD60_pred)
cm_DRD60 = confusion_matrix(y_DRD60, y_DRD60_pred, labels=np.arange(n_classes))
cm_DRD60_normalized = cm_DRD60.astype(float) / cm_DRD60.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_DRD60_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in DRD60 Test-set")
plt.tight_layout()
plt.savefig("DRD60_Confusion_final.pdf", dpi=300)

y_UMI_pred = []

for prob in y_UMI_proba:
    over_thresh = np.where(prob >= threshold)[0]
    
    if over_thresh.size > 0:
        # Pick the one with the highest probability among those over the threshold
        pred_class = over_thresh[np.argmax(prob[over_thresh])]
    else:
        # Fallback to normal dN
        pred_class = 0
    
    y_UMI_pred.append(pred_class)

y_UMI_pred = np.array(y_UMI_pred)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_UMI_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in Mixed UMI Test-set")
plt.tight_layout()
plt.savefig("Mixed-SRICO_Confusion_final.pdf", dpi=300)

# Feature correlation
X_train_bal = pd.DataFrame(X_train_bal, columns=features)
corr_matrix = X_train_bal.corr()
custom_cmap = LinearSegmentedColormap.from_list(
    "custom_corr", ["#698A9C", "white", "#DB7176"]
)

plt.figure(figsize=(8,8))
sns.heatmap(corr_matrix, annot=False, fmt=".2f", cmap=custom_cmap, 
            square=True, linewidths=0.5, cbar_kws={"shrink": 0.8})

plt.title("Feature Correlation", fontsize=16, fontweight='bold')
plt.xticks(
    np.arange(len(features)) + 0.5,  
    features,
    rotation=45,
    ha='right',
    size=12,  
    fontweight='normal'
)

plt.yticks(
    np.arange(len(features)) + 0.5,  
    features,
    rotation=0,
    size=12,  
    fontweight='normal'
)

plt.tight_layout()
plt.savefig("Feature_Correlation_final.pdf", dpi=300)

custom_labels = {
    0: "dN",
    1: "rA",
    2: "rC",
    3: "rG",
    4: "rU"
}
class_names_list = [custom_labels[i] for i in range(len(custom_labels))]
custom_colors = ['#A7C7E7', '#C0B6A4', '#878580', '#698A9C', '#DB7176']

explainer = shap.TreeExplainer(final_model)
X_df = data_withRinfo[features]
X_sample = X_df.sample(n=5000, random_state=42)
#X_sample_np = X_sample.to_numpy()
shap_values = explainer.shap_values(X_sample)

shap.summary_plot(
    shap_values, X_sample,
    plot_type="bar",
    class_names=class_names_list, 
    show=False
)
ax = plt.gca()
bars = ax.patches

num_features = len(bars) // len(custom_colors)

for class_idx in range(len(custom_colors)):
    start_idx = class_idx * num_features
    end_idx = (class_idx + 1) * num_features
    for bar in bars[start_idx:end_idx]:
        bar.set_color(custom_colors[class_idx])

plt.legend(loc='lower right', frameon=True)
plt.xlabel("Mean SHAP Value")
plt.ylabel("Feature")
plt.tight_layout()
plt.savefig("Feature_Importance_final.pdf", dpi=300)

# Try lower rN ratio
data_UMI_withRinfo = pd.merge(data_UMI, centralR, on = 'contig', how = 'left')
data_teln_realctl_withRinfo = pd.merge(data_teln_realctl, centralR, on = 'contig', how = 'left')
combined_UMI = pd.concat([data_UMI_withRinfo, data_teln_realctl_withRinfo], axis=0, ignore_index=True)

conditions = [
    combined_UMI['Type'] == 0,
    combined_UMI['cen_R'] == 'A',
    combined_UMI['cen_R'] == 'C',
    combined_UMI['cen_R'] == 'G',
    combined_UMI['cen_R'] == 'T'
]

combined_UMI['rType'] = np.select(conditions, choices, default=np.nan)

X_UMI = combined_UMI[features].to_numpy().astype(np.float32)
y_UMI = combined_UMI['rType'].to_numpy().astype(np.int32)
# Original data
X = combined_UMI[features]
y = combined_UMI['rType']
# Separate class 1 and others
X_class1 = X[y == 0]
X_others = X[y != 0]
y_class1 = y[y == 0]
y_others = y[y != 0]

# Higher dN ratio
desired_ratio = 200
n_class1 = len(y_class1)
n_others_required = int(n_class1 / desired_ratio)

if len(y_others) > n_others_required:
    X_others_sampled = X_others.sample(n=n_others_required, random_state=42)
    y_others_sampled = y_others.loc[X_others_sampled.index]
else:
    X_others_sampled = X_others
    y_others_sampled = y_others

X_balanced = pd.concat([X_class1, X_others_sampled])
y_balanced = pd.concat([y_class1, y_others_sampled])

from sklearn.utils import shuffle
X_balanced, y_balanced = shuffle(X_balanced, y_balanced, random_state=42)

X_UMI = X_balanced.to_numpy().astype(np.float32)
y_UMI = y_balanced.to_numpy().astype(np.int32)

from sklearn.preprocessing import label_binarize
y_UMI_bin = label_binarize(y_UMI, classes=[0, 1, 2, 3, 4])
y_UMI_proba = final_model.predict_proba(X_UMI)
y_UMI_pred = final_model.predict(X_UMI)
# ROC per class
fpr_UMI = np.linspace(0, 1, 100)
tpr_UMI_list = []
auc_UMI_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_UMI_bin[:, i], y_UMI_proba[:, i])
    interp_tpr = np.interp(fpr_UMI, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_UMI_list.append(interp_tpr)
    auc_UMI_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_UMI = np.mean(tpr_UMI_list, axis=0)
auc_macro_UMI = np.mean(auc_UMI_list)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)

plt.figure(figsize=(7, 7))
plt.plot(fpr_DRD60, mean_tpr_DRD60,
         label=f"DRD60 Test-set; rN:dN = 1:200 (AUC = {auc_macro_DRD60:.2f})",
         color="#DB7176", lw=2)

plt.plot(fpr_UMI, mean_tpr_UMI,
         label=f"Mixed SRICO Test-set; rN:dN = 1:200 (AUC = {auc_macro_UMI:.2f})",
         color="#878580", lw=2)

plt.plot([0, 1], [0, 1], 'k--', lw=1)

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Macro-average ROC for Multiclass XGboost')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("ROC_lower_rN_ratio.pdf", dpi=300)

custom_color = "#DB7176"
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_blue", ["#ffffff", custom_color])

# Setup
custom_labels = {
    0: "dN",
    1: "rA",
    2: "rC",
    3: "rG",
    4: "rU"
}
label_names = [custom_labels[i] for i in range(len(custom_labels))]

threshold = 0.5
y_DRD60_pred = []

for prob in y_DRD60_proba:
    over_thresh = np.where(prob >= threshold)[0]
    
    if over_thresh.size > 0:
        pred_class = over_thresh[np.argmax(prob[over_thresh])]
    else:
        # Fallback to dN
        pred_class = 0
    
    y_DRD60_pred.append(pred_class)

y_DRD60_pred = np.array(y_DRD60_pred)
cm_DRD60 = confusion_matrix(y_DRD60, y_DRD60_pred, labels=np.arange(n_classes))
cm_DRD60_normalized = cm_DRD60.astype(float) / cm_DRD60.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_DRD60_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("DRD60 Test-set; rN:dN = 1:200")
plt.tight_layout()
plt.savefig("DRD60_Confusion_lower_rN_ratio.pdf")

y_UMI_pred = []

for prob in y_UMI_proba:
    over_thresh = np.where(prob >= threshold)[0]
    
    if over_thresh.size > 0:
        pred_class = over_thresh[np.argmax(prob[over_thresh])]
    else:
        # Fallback to dN
        pred_class = 0
    
    y_UMI_pred.append(pred_class)

y_UMI_pred = np.array(y_UMI_pred)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_UMI_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mixed UMI Test-set; rN:dN = 1:200")
plt.tight_layout()
plt.savefig("Mixed_SRICO_Confusion_lower_rN_ratio.pdf")

# Plot F1 score over ratio of rN
data_UMI_withRinfo = pd.merge(data_UMI, centralR, on = 'contig', how = 'left')
data_teln_realctl_withRinfo = pd.merge(data_teln_realctl, centralR, on = 'contig', how = 'left')
combined_UMI = pd.concat([data_UMI_withRinfo, data_teln_realctl_withRinfo], axis=0, ignore_index=True)

conditions = [
    combined_UMI['Type'] == 0,
    combined_UMI['cen_R'] == 'A',
    combined_UMI['cen_R'] == 'C',
    combined_UMI['cen_R'] == 'G',
    combined_UMI['cen_R'] == 'T'
]

choices = [0, 1, 2, 3, 4]
combined_UMI['rType'] = np.select(conditions, choices, default=np.nan)

combined_UMI.head(5)

X_UMI = combined_UMI[features].to_numpy().astype(np.float32)
y_UMI = combined_UMI['rType'].to_numpy().astype(np.int32)

X = combined_UMI[features]
y = combined_UMI['rType']

def safe_indexing(X, idx):
    if isinstance(X, (pd.DataFrame, pd.Series)):
        return X.iloc[idx]
    else:  
        return X[idx]
    
ratios = [1, 10, 20, 50, 100, 200]
thresholds = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
n_repeats = 1

results_overall = {thr: [] for thr in thresholds}
results_overall_ci = {thr: [] for thr in thresholds}

for threshold in thresholds:
    for desired_ratio in ratios:
        f1all_repeats = []

        for rep in range(n_repeats):
            X_class0 = X[y == 0]
            X_others = X[y != 0]
            y_class0 = y[y == 0]
            y_others = y[y != 0]

            n_others = len(y_others)
            n_class0_required = int(n_others * desired_ratio)

            rng = np.random.RandomState(rep)
            idx0 = rng.choice(len(y_class0), n_class0_required, replace=True)
            X_class0_up = safe_indexing(X_class0, idx0)
            y_class0_up = safe_indexing(y_class0, idx0)

            if isinstance(X, pd.DataFrame):
                X_balanced = pd.concat([X_class0_up, X_others])
                y_balanced = pd.concat([y_class0_up, y_others])
            else:
                X_balanced = np.vstack([X_class0_up, X_others])
                y_balanced = np.concatenate([y_class0_up, y_others])

            X_balanced, y_balanced = shuffle(X_balanced, y_balanced, random_state=rep)

            X_eval = X_balanced.to_numpy().astype(np.float32)
            y_eval = y_balanced.to_numpy().astype(np.int32)

            y_proba = final_model.predict_proba(X_eval)
            y_pred = np.argmax(y_proba, axis=1)
            max_proba = np.max(y_proba, axis=1)
            y_pred[max_proba < threshold] = 0

            f1all = f1_score(y_eval, y_pred, average="macro", zero_division=0)
            f1all_repeats.append(f1all)

        results_overall[threshold].append(np.mean(f1all_repeats))
        results_overall_ci[threshold].append(
            (np.min(f1all_repeats), np.max(f1all_repeats))
        )

colors = {
    0.6: "#878580",  
    0.95: "#DB7176",  
    0.5: "#698a9c",  
    0.7: "#edd8af",  
    0.8: "#A7C7E7",   
    0.9: "#9abe8a"   
}

plt.figure(figsize=(5, 5))

for i, threshold in enumerate(thresholds):
    means = results_overall[threshold]
    ci = results_overall_ci[threshold]

    plt.plot(ratios, means, lw=2, label=f"{threshold}", color=colors[threshold])
    plt.fill_between(
        ratios,
        [low for (low, high) in ci],
        [high for (low, high) in ci],
        alpha=0.15,
        color=colors[threshold]
    )

plt.ylim(0, 1)
plt.xlabel("Ratio dN:rN")
plt.ylabel("Macro F1 (all classes)")
plt.grid(True, alpha=0.25)
plt.legend(title="Threshold", loc="lower left")
plt.tight_layout()
plt.savefig("F1_score_over_ratio_Mixed_SRICO_range.pdf")

# Include 20x more training data 
folder_path_r1 = pro_r1
folder_path_r2 = pro_r2
folder_path_minion = "./"

patterns = [f"{folder_path_r1}*.tsv", f"{folder_path_r2}*.tsv", f"{folder_path_minion}*.tsv"]

ddf = dd.read_csv(patterns, sep='\t')
ddf.columns = [
    'read_name', 'contig',
    'ampl_mean_1st', 'ampl_mean_2nd', 'ampl_mean_3rd', 'ampl_mean_4th', 'ampl_mean_5th',
    'ampl_mean_6th', 'ampl_mean_7th', 'ampl_mean_8th', 'ampl_mean_9th',
    'model_deviation_1st', 'model_deviation_2nd', 'model_deviation_3rd', 'model_deviation_4th',
    'model_deviation_5th', 'model_deviation_6th', 'model_deviation_7th', 'model_deviation_8th',
    'model_deviation_9th',
    'Type', 'dwell_time_R-18', 'dwell_time_R-15'
]
centralR_design = centralR
ddf_merged = ddf.merge(centralR_design, on='contig', how='left')

def assign_rType(df):
    df['cen_R'] = df['cen_R'].fillna('')

    mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    df['rType'] = df['cen_R'].map(mapping)

    df.loc[df['Type'] == 0, 'rType'] = 0

    return df
meta = ddf_merged._meta.copy()
meta['rType'] = pd.Series(dtype='float64')

ddf_new = ddf_merged.map_partitions(assign_rType, meta=meta)

ddf_new = ddf_new[~ddf_new['read_name'].str.contains(r'_downc', regex=True)]

final_df = ddf_new.compute()
data_withRinfo = final_df

# Features and target
features = ["ampl_mean_2nd","ampl_mean_3rd","ampl_mean_4th", "ampl_mean_5th", "ampl_mean_6th", "ampl_mean_7th",    
            "model_deviation_2nd","model_deviation_3rd","model_deviation_4th", "model_deviation_5th", "model_deviation_6th", "model_deviation_7th", 
            ]
X = data_withRinfo[features].to_numpy().astype(np.float32)
y = data_withRinfo['rType'].to_numpy().astype(np.int32)

X_train, X_test, y_train, y_test = gpu_train_test_split(X, y, test_size=0.2, random_state=42)
X_train_df = pd.DataFrame(X_train, columns=features)
y_train_df = pd.Series(y_train, name='rType')

# Combine features + target
df = pd.concat([X_train_df, y_train_df], axis=1)

# Balance class size
min_size = df['rType'].value_counts().min()

df_balanced = (
    df.groupby('rType')
      .apply(lambda x: x.sample(min_size, random_state=42))
      .reset_index(drop=True)
)

X_train_bal = df_balanced.drop(columns=['rType'])
y_train_bal = df_balanced['rType']

# Train
classes = np.unique(y_train_bal)
n_classes = len(classes)
best_params = {'colsample_bytree': np.float64(0.7721), 'gamma': np.float64(0.3028), 'learning_rate': np.float64(0.1768), 'max_depth': np.float64(8.81), 'min_child_weight': np.float64(5.81), 'n_estimators': np.float64(1.182e+03), 'reg_alpha': np.float64(8.457), 'reg_lambda': np.float64(6.806), 'subsample': np.float64(0.7801)}
best_model = XGBClassifier(
    max_depth=int(best_params['max_depth']),
    learning_rate=best_params['learning_rate'],
    n_estimators=int(best_params['n_estimators']),
    gamma=best_params['gamma'],
    min_child_weight=best_params['min_child_weight'],
    subsample=best_params['subsample'],
    colsample_bytree=best_params['colsample_bytree'],
    reg_alpha=best_params['reg_alpha'],
    reg_lambda=best_params['reg_lambda'],
    tree_method='hist',
    device='cuda',
    eval_metric='mlogloss',
    random_state=42
)

final_model = best_model
final_model.fit(X_train_bal, y_train_bal)

y_test_bin = label_binarize(y_test, classes=[0, 1, 2, 3, 4])
y_test_proba = final_model.predict_proba(X_test)
y_test_pred = final_model.predict(X_test)

# ROC per class
fpr_test = np.linspace(0, 1, 100)
tpr_test_list = []
auc_test_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_test_proba[:, i])
    interp_tpr = np.interp(fpr_test, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_test_list.append(interp_tpr)
    auc_test_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_test = np.mean(tpr_test_list, axis=0)
auc_macro_test = np.mean(auc_test_list)
cm_test = confusion_matrix(y_test, y_test_pred, labels=np.arange(n_classes))
cm_test_normalized = cm_test.astype(float) / cm_test.sum(axis=1, keepdims=True)

# DRD60
data_DRD60_withRinfo = pd.merge(data_DRD60, centralR_design, on = 'contig', how = 'left')
conditions = [
    data_DRD60_withRinfo['Type'] == 0,
    data_DRD60_withRinfo['cen_R'] == 'A',
    data_DRD60_withRinfo['cen_R'] == 'C',
    data_DRD60_withRinfo['cen_R'] == 'G',
    data_DRD60_withRinfo['cen_R'] == 'T'
]

choices = [0, 1, 2, 3, 4]
data_DRD60_withRinfo['rType'] = np.select(conditions, choices, default=np.nan)

data_DRD60_withRinfo.head(5)
X_DRD60 = data_DRD60_withRinfo[features].to_numpy().astype(np.float32)
y_DRD60 = data_DRD60_withRinfo['rType'].to_numpy().astype(np.int32)
y_DRD60_bin = label_binarize(y_DRD60, classes=[0, 1, 2, 3, 4])
y_DRD60_proba = final_model.predict_proba(X_DRD60)
y_DRD60_pred = final_model.predict(X_DRD60)

fpr_DRD60 = np.linspace(0, 1, 100)
tpr_DRD60_list = []
auc_DRD60_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_DRD60_bin[:, i], y_DRD60_proba[:, i])
    interp_tpr = np.interp(fpr_DRD60, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_DRD60_list.append(interp_tpr)
    auc_DRD60_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_DRD60 = np.mean(tpr_DRD60_list, axis=0)
auc_macro_DRD60 = np.mean(auc_DRD60_list)
cm_DRD60 = confusion_matrix(y_DRD60, y_DRD60_pred, labels=np.arange(n_classes))
cm_DRD60_normalized = cm_DRD60.astype(float) / cm_DRD60.sum(axis=1, keepdims=True)

data_teln_realctl_withRinfo = pd.merge(data_teln_realctl, centralR_design, on = 'contig', how = 'left')
combined_UMI = pd.concat([data_UMI_withRinfo, data_teln_realctl_withRinfo], axis=0, ignore_index=True)

conditions = [
    combined_UMI['Type'] == 0,
    combined_UMI['cen_R'] == 'A',
    combined_UMI['cen_R'] == 'C',
    combined_UMI['cen_R'] == 'G',
    combined_UMI['cen_R'] == 'T'
]

choices = [0, 1, 2, 3, 4]
combined_UMI['rType'] = np.select(conditions, choices, default=np.nan)

X_UMI = combined_UMI[features].to_numpy().astype(np.float32)
y_UMI = combined_UMI['rType'].to_numpy().astype(np.int32)

y_UMI_bin = label_binarize(y_UMI, classes=[0, 1, 2, 3, 4])
y_UMI_proba = final_model.predict_proba(X_UMI)
y_UMI_pred = final_model.predict(X_UMI)

# ROC per class
fpr_UMI = np.linspace(0, 1, 100)
tpr_UMI_list = []
auc_UMI_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_UMI_bin[:, i], y_UMI_proba[:, i])
    interp_tpr = np.interp(fpr_UMI, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_UMI_list.append(interp_tpr)
    auc_UMI_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_UMI = np.mean(tpr_UMI_list, axis=0)
auc_macro_UMI = np.mean(auc_UMI_list)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)

# Plot
plt.figure(figsize=(7, 7))

# Test ROC
plt.plot(fpr_test, mean_tpr_test,
         label=f"Hold-out Test-set (AUC = {auc_macro_test:.2f})",
         color="#C0B6A4", lw=2)

plt.plot(fpr_DRD60, mean_tpr_DRD60,
         label=f"DRD60 Test-set (AUC = {auc_macro_DRD60:.2f})",
         color="#DB7176", lw=2)

plt.plot(fpr_UMI, mean_tpr_UMI,
         label=f"Mixed SRICO Test-set (AUC = {auc_macro_UMI:.2f})",
         color="#878580", lw=2)


# Diagonal line
plt.plot([0, 1], [0, 1], 'k--', lw=1)

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Macro-average ROC for Multiclass XGboost')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("./ROC_Prome.pdf")

# Confusion Matrix for Mixed-SRICO set
threshold = 0.5
y_UMI_pred = []

for prob in y_UMI_proba:
    over_thresh = np.where(prob >= threshold)[0]
    
    if over_thresh.size > 0:
        pred_class = over_thresh[np.argmax(prob[over_thresh])]
    else:
        pred_class = 0
    
    y_UMI_pred.append(pred_class)

y_UMI_pred = np.array(y_UMI_pred)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_UMI_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in Mixed UMI Test-set")
plt.tight_layout()
plt.savefig("./CM_prom_Mixed_SRICO.pdf")

# Plot recalls over the sizes of training set
ratios = [1, 0.1, 0.01, 0.001, 0.0001, 0]  
abs_num = [int(len(X_train_bal) * r) for r in ratios]

n_repeats = 10  

# Store per-class results
recalls_mean_UMI = {cls: [] for cls in range(5)}
recalls_range_UMI = {cls: [] for cls in range(5)}

for r in ratios:
    if r == 0:
        # baseline: random (no training)
        for cls in range(5):
            recalls_mean_UMI[cls].append(0.0)
            recalls_range_UMI[cls].append((0.0, 0.0))
        continue

    n_samples = max(int(len(X_train_bal) * r), len(np.unique(y_train_bal)))

    # store recall results across repeats
    recall_repeats_UMI = {cls: [] for cls in range(5)}

    for rep in range(n_repeats):
        X_train_sub, y_train_sub = resample(
            X_train_bal, y_train_bal,
            replace=False,
            n_samples=n_samples,
            stratify=y_train_bal if n_samples >= len(np.unique(y_train_bal)) else None,
            random_state=rep 
        )

        # Train model
        model = best_model
        model.fit(X_train_sub, y_train_sub)

        y_proba_UMI = model.predict_proba(X_UMI)
        y_pred_UMI = np.argmax(y_proba_UMI, axis=1)
        max_proba_UMI = np.max(y_proba_UMI, axis=1)
        y_pred_UMI[max_proba_UMI < 0.5] = 0

        # Confusion matrix
        cm_UMI = confusion_matrix(y_UMI, y_pred_UMI, labels=[0, 1, 2, 3, 4])

        # Compute recall per class
        for i in range(5):
            TP_UMI = cm_UMI[i, i]
            FN_UMI = cm_UMI[i, :].sum() - TP_UMI
            recall_UMI = TP_UMI / (TP_UMI + FN_UMI) if (TP_UMI + FN_UMI) > 0 else 0
            recall_repeats_UMI[i].append(recall_UMI)

    for i in range(5):
        recalls_mean_UMI[i].append(np.mean(recall_repeats_UMI[i]))
        recalls_range_UMI[i].append((np.min(recall_repeats_UMI[i]), np.max(recall_repeats_UMI[i])))

colors = {
    0: "#878580",
    1: "#DB7176",
    2: "#698a9c",
    3: "#edd8af",
    4: "#9abe8a"
}
labels = {0: "dN", 1: "rA", 2: "rC", 3: "rG", 4: "rU"}

plt.figure(figsize=(8, 6))
for cls, mean_vals in recalls_mean_UMI.items():
    plt.plot(abs_num, mean_vals, lw=2, color=colors[cls], label=labels[cls])
    plt.fill_between(
        abs_num,
        [low for (low, high) in recalls_range_UMI[cls]],
        [high for (low, high) in recalls_range_UMI[cls]],
        color=colors[cls],
        alpha=0.2
    )

plt.xscale("log")
plt.axvline(x=686743 * 0.8, color="black", linestyle="--", linewidth=1.5, alpha=0.8)
plt.ylim(0, 1)
plt.xlabel("Training-set Size (log scale)")
plt.ylabel("Recall")
plt.title("Recalls of rN (Prob.≥0.5) in Mixed SRICO test-set upon down-sampled training-set")
plt.grid(True)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("./Downsampling_Recalls_MixedSRICO.pdf")

joblib.dump(final_model, '../../Single_rNMP_predictor/XGB_SRICO_1041.pkl')

# Include 20x promethion training data with 40,000 Wild-type yeast genomic DNA data
wild= pd.read_csv(wild_ygDNA, delimiter = '\t')
wild_index = wild.rename(columns={'Type': 'rType'})
wild_index_sampled=wild_index.sample(n=40000, random_state=42)
merged_data = pd.concat([data_withRinfo, wild_index_sampled], ignore_index=True)

features = ["ampl_mean_2nd","ampl_mean_3rd","ampl_mean_4th", "ampl_mean_5th", "ampl_mean_6th", "ampl_mean_7th",    
            "model_deviation_2nd","model_deviation_3rd","model_deviation_4th", "model_deviation_5th", "model_deviation_6th", "model_deviation_7th", 
            ]
X = merged_data[features].to_numpy().astype(np.float32)
y = merged_data['rType'].to_numpy().astype(np.int32)
X_train, X_test, y_train, y_test = gpu_train_test_split(X, y, test_size=0.2, random_state=42)
X_train_df = pd.DataFrame(X_train, columns=features)
y_train_df = pd.Series(y_train, name='rType')

# Combine features + target
df = pd.concat([X_train_df, y_train_df], axis=1)

# Balance class size
min_size = df['rType'].value_counts().min()

df_balanced = (
    df.groupby('rType')
      .apply(lambda x: x.sample(min_size, random_state=42))
      .reset_index(drop=True)
)

X_train_bal = df_balanced.drop(columns=['rType'])
y_train_bal = df_balanced['rType']

# Cross-validation setup
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
fpr_grid = np.linspace(0, 1, 100)
tpr_all = [] 
auc_folds = []

X_train_bal_cp = cp.asarray(X_train_bal.values)
y_train_bal_cp = cp.asarray(y_train_bal.values)
y_train_bal_np = y_train_bal.values  

for train_idx, val_idx in cv.split(X_train_bal, y_train_bal_np):
    # Use indices to slice CuPy arrays
    X_tr = cp.take(X_train_bal_cp, train_idx, axis=0)
    X_val = cp.take(X_train_bal_cp, val_idx, axis=0)
    y_tr = cp.take(y_train_bal_cp, train_idx)
    y_val = cp.take(y_train_bal_cp, val_idx)

    y_val_cpu = cp.asnumpy(y_val)
    y_val_bin = label_binarize(y_val_cpu, classes=classes)

    model = best_model
    model.fit(X_tr, y_tr)

    y_score = cp.asnumpy(model.predict_proba(X_val))

    for i in range(n_classes):
        fpr, tpr, _ = roc_curve(y_val_bin[:, i], y_score[:, i])
        interp_tpr = np.interp(fpr_grid, fpr, tpr)
        interp_tpr[0] = 0.0
        tpr_all.append(interp_tpr)

    auc_macro = np.mean([
        auc(roc_curve(y_val_bin[:, i], y_score[:, i])[0],
            roc_curve(y_val_bin[:, i], y_score[:, i])[1])
        for i in range(n_classes)
    ])
    auc_folds.append(auc_macro)

    y_pred = model.predict(X_val)
  In [18]:  # Convert and aggregate
tpr_all = np.array(tpr_all)  
mean_tpr = np.mean(tpr_all, axis=0)
std_tpr = np.std(tpr_all, axis=0, ddof=1)
se_tpr = std_tpr / np.sqrt(tpr_all.shape[0])
ci_upper = np.clip(mean_tpr + 1.96 * se_tpr, 0, 1)
ci_lower = np.clip(mean_tpr - 1.96 * se_tpr, 0, 1)
mean_auc = np.mean(auc_folds)

# Train
classes = np.unique(y_train_bal)
n_classes = len(classes)
best_params = {'colsample_bytree': np.float64(0.7721), 'gamma': np.float64(0.3028), 'learning_rate': np.float64(0.1768), 'max_depth': np.float64(8.81), 'min_child_weight': np.float64(5.81), 'n_estimators': np.float64(1.182e+03), 'reg_alpha': np.float64(8.457), 'reg_lambda': np.float64(6.806), 'subsample': np.float64(0.7801)}
best_model = XGBClassifier(
    max_depth=int(best_params['max_depth']),
    learning_rate=best_params['learning_rate'],
    n_estimators=int(best_params['n_estimators']),
    gamma=best_params['gamma'],
    min_child_weight=best_params['min_child_weight'],
    subsample=best_params['subsample'],
    colsample_bytree=best_params['colsample_bytree'],
    reg_alpha=best_params['reg_alpha'],
    reg_lambda=best_params['reg_lambda'],
    tree_method='hist',
    device='cuda',
    eval_metric='mlogloss',
    random_state=42
)
  
final_model = best_model
final_model.fit(X_train_bal, y_train_bal)

y_test_bin = label_binarize(y_test, classes=[0, 1, 2, 3, 4])
y_test_proba = final_model.predict_proba(X_test)
y_test_pred = final_model.predict(X_test)

# ROC per class
mixed_20X_small_fpr_test = np.linspace(0, 1, 100)
mixed_20X_small_tpr_test_list = []
mixed_20X_small_auc_test_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_test_proba[:, i])
    interp_tpr = np.interp(mixed_20X_small_fpr_test, fpr, tpr)
    interp_tpr[0] = 0.0
    mixed_20X_small_tpr_test_list.append(interp_tpr)
    mixed_20X_small_auc_test_list.append(auc(fpr, tpr))

# Macro-average ROC
mixed_20X_small_mean_tpr_test = np.mean(mixed_20X_small_tpr_test_list, axis=0)
mixed_20X_small_auc_macro_test = np.mean(mixed_20X_small_auc_test_list)
cm_test = confusion_matrix(y_test, y_test_pred, labels=np.arange(n_classes))
cm_test_normalized = cm_test.astype(float) / cm_test.sum(axis=1, keepdims=True)

# Feature Correlation 
X_train_bal = pd.DataFrame(X_train_bal, columns=features)
corr_matrix = X_train_bal.corr()
custom_cmap = LinearSegmentedColormap.from_list(
    "custom_corr", ["#698A9C", "white", "#DB7176"]
)
plt.figure(figsize=(8,8))
sns.heatmap(corr_matrix, annot=False, fmt=".2f", cmap=custom_cmap, 
            square=True, linewidths=0.5, cbar_kws={"shrink": 0.8})
plt.title("Feature Correlation", fontsize=16, fontweight='bold')
plt.xticks(
    np.arange(len(features)) + 0.5,  
    features,
    rotation=45,
    ha='right',
    size=10,  
    fontweight='normal'
)
plt.yticks(
    np.arange(len(features)) + 0.5,  
    features,
    rotation=0,
    size=10,  
    fontweight='normal'
)
plt.tight_layout()
plt.savefig("Feature_correlation_final_ygDNA.pdf", dpi=300)

# Feature Importance
custom_labels = {
    0: "dN",
    1: "rA",
    2: "rC",
    3: "rG",
    4: "rU"
}
class_names_list = [custom_labels[i] for i in range(len(custom_labels))]
custom_colors = ['#A7C7E7', '#C0B6A4', '#878580', '#698A9C', '#DB7176']

explainer = shap.TreeExplainer(final_model)
X_df = merged_data[features]
X_sample = X_df.sample(n=5000, random_state=42)
#X_sample_np = X_sample.to_numpy()
shap_values = explainer.shap_values(X_sample)

shap.summary_plot(
    shap_values, X_sample,
    plot_type="bar",
    class_names=class_names_list, 
    show=False
)
# Get the axis and bars
ax = plt.gca()
bars = ax.patches

# Number of features (assumes all classes have the same number of bars)
num_features = len(bars) // len(custom_colors)

# Assign colors per class (not per bar)
for class_idx in range(len(custom_colors)):
    # Get all bars for the current class
    start_idx = class_idx * num_features
    end_idx = (class_idx + 1) * num_features
    for bar in bars[start_idx:end_idx]:
        bar.set_color(custom_colors[class_idx])

# Update legend colors to match
#legend = ax.legend_
#for legend_handle, color in zip(legend.legendHandles, custom_colors):
    #legend_handle.set_color(color)

# Optional: move legend to bottom right
plt.legend(loc='lower right', frameon=True)
plt.xlabel("Mean SHAP Value")
plt.ylabel("Feature")
plt.tight_layout()
plt.savefig("Feature_importance_final_ygDNA.pdf", dpi=300)

# DRD60
data_DRD60_withRinfo = pd.merge(data_DRD60, centralR_design, on = 'contig', how = 'left')
conditions = [
    data_DRD60_withRinfo['Type'] == 0,
    data_DRD60_withRinfo['cen_R'] == 'A',
    data_DRD60_withRinfo['cen_R'] == 'C',
    data_DRD60_withRinfo['cen_R'] == 'G',
    data_DRD60_withRinfo['cen_R'] == 'T'
]
choices = [0, 1, 2, 3, 4]
data_DRD60_withRinfo['rType'] = np.select(conditions, choices, default=np.nan)

X_DRD60 = data_DRD60_withRinfo[features].to_numpy().astype(np.float32)
y_DRD60 = data_DRD60_withRinfo['rType'].to_numpy().astype(np.int32)
y_DRD60_bin = label_binarize(y_DRD60, classes=[0, 1, 2, 3, 4])
y_DRD60_proba = final_model.predict_proba(X_DRD60)
y_DRD60_pred = final_model.predict(X_DRD60)

fpr_DRD60 = np.linspace(0, 1, 100)
tpr_DRD60_list = []
auc_DRD60_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_DRD60_bin[:, i], y_DRD60_proba[:, i])
    interp_tpr = np.interp(fpr_DRD60, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_DRD60_list.append(interp_tpr)
    auc_DRD60_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_DRD60 = np.mean(tpr_DRD60_list, axis=0)
auc_macro_DRD60 = np.mean(auc_DRD60_list)
cm_DRD60 = confusion_matrix(y_DRD60, y_DRD60_pred, labels=np.arange(n_classes))
cm_DRD60_normalized = cm_DRD60.astype(float) / cm_DRD60.sum(axis=1, keepdims=True)

# Mixed-SRICO 
data_teln_realctl_withRinfo = pd.merge(data_teln_realctl, centralR_design, on = 'contig', how = 'left')
combined_UMI = pd.concat([data_UMI_withRinfo, data_teln_realctl_withRinfo], axis=0, ignore_index=True)

conditions = [
    combined_UMI['Type'] == 0,
    combined_UMI['cen_R'] == 'A',
    combined_UMI['cen_R'] == 'C',
    combined_UMI['cen_R'] == 'G',
    combined_UMI['cen_R'] == 'T'
]
choices = [0, 1, 2, 3, 4]
combined_UMI['rType'] = np.select(conditions, choices, default=np.nan)

X_UMI = combined_UMI[features].to_numpy().astype(np.float32)
y_UMI = combined_UMI['rType'].to_numpy().astype(np.int32)

y_UMI_bin = label_binarize(y_UMI, classes=[0, 1, 2, 3, 4])
y_UMI_proba = final_model.predict_proba(X_UMI)
y_UMI_pred = final_model.predict(X_UMI)

# ROC per class
fpr_UMI = np.linspace(0, 1, 100)
tpr_UMI_list = []
auc_UMI_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_UMI_bin[:, i], y_UMI_proba[:, i])
    interp_tpr = np.interp(fpr_UMI, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_UMI_list.append(interp_tpr)
    auc_UMI_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_UMI = np.mean(tpr_UMI_list, axis=0)
auc_macro_UMI = np.mean(auc_UMI_list)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)

# Plot ROC
plt.figure(figsize=(7, 7))

# Test ROC
plt.plot(fpr_grid, mean_tpr,
         label=f"10-fold CV with 95% CI (Mean AUC = {mean_auc:.2f})",
         color="#698A9C", lw=2)

plt.plot(mixed_20X_small_fpr_test, mixed_20X_small_mean_tpr_test,
         label=f"Hold-out Test-set (AUC = {mixed_20X_small_auc_macro_test:.2f})",
         color="#C0B6A4", lw=2)

plt.plot(fpr_DRD60, mean_tpr_DRD60,
         label=f"DRD60 Test-set (AUC = {auc_macro_DRD60:.2f})",
         color="#DB7176", lw=2)

plt.plot(fpr_UMI, mean_tpr_UMI,
         label=f"Mixed SRICO Test-set (AUC = {auc_macro_UMI:.2f})",
         color="#878580", lw=2)

plt.fill_between(fpr_grid, ci_lower, ci_upper, color="#698A9C", alpha=0.4)

# Diagonal line
plt.plot([0, 1], [0, 1], 'k--', lw=1)

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Macro-average ROC for Multiclass XGboost')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("ROC_prome_final_ygDNA.pdf")

# Confusion Matrix for Hold-out Test
# Setup
custom_labels = {
    0: "dN",
    1: "rA",
    2: "rC",
    3: "rG",
    4: "rU"
}
label_names = [custom_labels[i] for i in range(len(custom_labels))]
custom_color = "#DB7176"
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_blue", ["#ffffff", custom_color])
threshold = 0.5
y_HO_pred = []

for prob in y_test_proba:
    # If no prob exceeds threshold, fallback to argmax
    if (prob >= threshold).any():
        pred_class = np.argmax(prob * (prob >= threshold))
    else:
        pred_class = np.argmax(prob)
    y_HO_pred.append(pred_class)

y_HO_pred = np.array(y_HO_pred)
cm_HO = confusion_matrix(y_test, y_HO_pred, labels=np.arange(n_classes))
cm_HO_normalized = cm_HO.astype(float) / cm_HO.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_HO_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in Hold-out Test-set")
plt.tight_layout()
plt.savefig("CM_Holdout_Confusion_final_ygDNA.pdf", dpi=300)

# Confusion matrix for DRD60
y_DRD60_pred = []

for prob in y_DRD60_proba:
    # If no prob exceeds threshold, fallback to argmax
    if (prob >= threshold).any():
        pred_class = np.argmax(prob * (prob >= threshold))
    else:
        pred_class = np.argmax(prob)
    y_DRD60_pred.append(pred_class)

y_DRD60_pred = np.array(y_DRD60_pred)
cm_DRD60 = confusion_matrix(y_DRD60, y_DRD60_pred, labels=np.arange(n_classes))
cm_DRD60_normalized = cm_DRD60.astype(float) / cm_DRD60.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_DRD60_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in DRD60 Test-set")
plt.tight_layout()
plt.savefig("CM_DRD60_final_ygDNA.pdf", dpi=300)

# Confusion Matrix for Mixed-SRICO
y_UMI_pred = []

for prob in y_UMI_proba:
    # Classes meeting the threshold
    over_thresh = np.where(prob >= threshold)[0]
    
    if over_thresh.size > 0:
        # Pick the one with the highest probability among those over the threshold
        pred_class = over_thresh[np.argmax(prob[over_thresh])]
    else:
        # Fallback to the usual best guess
        pred_class = 0
    
    y_UMI_pred.append(pred_class)

y_UMI_pred = np.array(y_UMI_pred)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_UMI_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mean Proportion in Mixed SRICO Test-set")
plt.tight_layout()
plt.savefig("CM_Mixed-SRICO_final_ygDNA.pdf", dpi=300)

# Try lower rN ratio
data_UMI_withRinfo = pd.merge(data_UMI, centralR, on = 'contig', how = 'left')
data_teln_realctl_withRinfo = pd.merge(data_teln_realctl, centralR, on = 'contig', how = 'left')
combined_UMI = pd.concat([data_UMI_withRinfo, data_teln_realctl_withRinfo], axis=0, ignore_index=True)

conditions = [
    combined_UMI['Type'] == 0,
    combined_UMI['cen_R'] == 'A',
    combined_UMI['cen_R'] == 'C',
    combined_UMI['cen_R'] == 'G',
    combined_UMI['cen_R'] == 'T'
]

combined_UMI['rType'] = np.select(conditions, choices, default=np.nan)

X_UMI = combined_UMI[features].to_numpy().astype(np.float32)
y_UMI = combined_UMI['rType'].to_numpy().astype(np.int32)
# Original data
X = combined_UMI[features]
y = combined_UMI['rType']
# Separate class 1 and others
X_class1 = X[y == 0]
X_others = X[y != 0]
y_class1 = y[y == 0]
y_others = y[y != 0]

# Higher dN ratio
desired_ratio = 200
n_class1 = len(y_class1)
n_others_required = int(n_class1 / desired_ratio)

if len(y_others) > n_others_required:
    X_others_sampled = X_others.sample(n=n_others_required, random_state=42)
    y_others_sampled = y_others.loc[X_others_sampled.index]
else:
    X_others_sampled = X_others
    y_others_sampled = y_others

X_balanced = pd.concat([X_class1, X_others_sampled])
y_balanced = pd.concat([y_class1, y_others_sampled])

from sklearn.utils import shuffle
X_balanced, y_balanced = shuffle(X_balanced, y_balanced, random_state=42)

X_UMI = X_balanced.to_numpy().astype(np.float32)
y_UMI = y_balanced.to_numpy().astype(np.int32)

from sklearn.preprocessing import label_binarize
y_UMI_bin = label_binarize(y_UMI, classes=[0, 1, 2, 3, 4])
y_UMI_proba = final_model.predict_proba(X_UMI)
y_UMI_pred = final_model.predict(X_UMI)

# ROC per class
fpr_UMI = np.linspace(0, 1, 100)
tpr_UMI_list = []
auc_UMI_list = []

for i in range(5):
    fpr, tpr, _ = roc_curve(y_UMI_bin[:, i], y_UMI_proba[:, i])
    interp_tpr = np.interp(fpr_UMI, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_UMI_list.append(interp_tpr)
    auc_UMI_list.append(auc(fpr, tpr))

# Macro-average ROC
mean_tpr_UMI = np.mean(tpr_UMI_list, axis=0)
auc_macro_UMI = np.mean(auc_UMI_list)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)

plt.figure(figsize=(7, 7))
plt.plot(fpr_DRD60, mean_tpr_DRD60,
         label=f"DRD60 Test-set (AUC = {auc_macro_DRD60:.2f})",
         color="#DB7176", lw=2)

plt.plot(fpr_UMI, mean_tpr_UMI,
         label=f"Mixed SRICO Test-set (AUC = {auc_macro_UMI:.2f})",
         color="#878580", lw=2)

plt.plot([0, 1], [0, 1], 'k--', lw=1)

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Macro-average ROC for Multiclass XGboost')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("ROC_lower_rN_ratio_final_ygDNA.pdf", dpi=300)

# Confusion matrix for Mixed-SCRICO in lower rN ratio
y_UMI_pred = []

for prob in y_UMI_proba:
    over_thresh = np.where(prob >= threshold)[0]
    
    if over_thresh.size > 0:
        pred_class = over_thresh[np.argmax(prob[over_thresh])]
    else:
        # Fallback to dN
        pred_class = 0
    
    y_UMI_pred.append(pred_class)

y_UMI_pred = np.array(y_UMI_pred)
cm_UMI = confusion_matrix(y_UMI, y_UMI_pred, labels=np.arange(n_classes))
cm_UMI_normalized = cm_UMI.astype(float) / cm_UMI.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_UMI_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Mixed UMI Test-set; rN:dN = 1:200")
plt.tight_layout()
plt.savefig("CM_Mixed_SRICO_lower_rN_ratio_final_ygDNA.pdf")

# Plot F1 score over ratio of rN
data_UMI_withRinfo = pd.merge(data_UMI, centralR, on = 'contig', how = 'left')
data_teln_realctl_withRinfo = pd.merge(data_teln_realctl, centralR, on = 'contig', how = 'left')
combined_UMI = pd.concat([data_UMI_withRinfo, data_teln_realctl_withRinfo], axis=0, ignore_index=True)

conditions = [
    combined_UMI['Type'] == 0,
    combined_UMI['cen_R'] == 'A',
    combined_UMI['cen_R'] == 'C',
    combined_UMI['cen_R'] == 'G',
    combined_UMI['cen_R'] == 'T'
]

choices = [0, 1, 2, 3, 4]
combined_UMI['rType'] = np.select(conditions, choices, default=np.nan)

combined_UMI.head(5)

X_UMI = combined_UMI[features].to_numpy().astype(np.float32)
y_UMI = combined_UMI['rType'].to_numpy().astype(np.int32)

X = combined_UMI[features]
y = combined_UMI['rType']

def safe_indexing(X, idx):
    if isinstance(X, (pd.DataFrame, pd.Series)):
        return X.iloc[idx]
    else:  
        return X[idx]
    
ratios = [1, 10, 20, 50, 100, 200]
thresholds = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
n_repeats = 1

results_overall = {thr: [] for thr in thresholds}
results_overall_ci = {thr: [] for thr in thresholds}

for threshold in thresholds:
    for desired_ratio in ratios:
        f1all_repeats = []

        for rep in range(n_repeats):
            X_class0 = X[y == 0]
            X_others = X[y != 0]
            y_class0 = y[y == 0]
            y_others = y[y != 0]

            n_others = len(y_others)
            n_class0_required = int(n_others * desired_ratio)

            rng = np.random.RandomState(rep)
            idx0 = rng.choice(len(y_class0), n_class0_required, replace=True)
            X_class0_up = safe_indexing(X_class0, idx0)
            y_class0_up = safe_indexing(y_class0, idx0)

            if isinstance(X, pd.DataFrame):
                X_balanced = pd.concat([X_class0_up, X_others])
                y_balanced = pd.concat([y_class0_up, y_others])
            else:
                X_balanced = np.vstack([X_class0_up, X_others])
                y_balanced = np.concatenate([y_class0_up, y_others])

            X_balanced, y_balanced = shuffle(X_balanced, y_balanced, random_state=rep)

            X_eval = X_balanced.to_numpy().astype(np.float32)
            y_eval = y_balanced.to_numpy().astype(np.int32)

            y_proba = final_model.predict_proba(X_eval)
            y_pred = np.argmax(y_proba, axis=1)
            max_proba = np.max(y_proba, axis=1)
            y_pred[max_proba < threshold] = 0

            f1all = f1_score(y_eval, y_pred, average="macro", zero_division=0)
            f1all_repeats.append(f1all)

        results_overall[threshold].append(np.mean(f1all_repeats))
        results_overall_ci[threshold].append(
            (np.min(f1all_repeats), np.max(f1all_repeats))
        )

colors = {
    0.6: "#878580",  
    0.95: "#DB7176",  
    0.5: "#698a9c",  
    0.7: "#edd8af",  
    0.8: "#A7C7E7",   
    0.9: "#9abe8a"   
}

plt.figure(figsize=(5, 5))

for i, threshold in enumerate(thresholds):
    means = results_overall[threshold]
    ci = results_overall_ci[threshold]

    plt.plot(ratios, means, lw=2, label=f"{threshold}", color=colors[threshold])
    plt.fill_between(
        ratios,
        [low for (low, high) in ci],
        [high for (low, high) in ci],
        alpha=0.15,
        color=colors[threshold]
    )

plt.ylim(0, 1)
plt.xlabel("Ratio dN:rN")
plt.ylabel("Macro F1 (all classes)")
plt.grid(True, alpha=0.25)
plt.legend(title="Threshold", loc="lower left")
plt.tight_layout()
plt.savefig("F1_score_over_ratio_Mixed_SRICO_range_final_ygDNA.pdf")

joblib.dump(final_model, '../../Single_rNMP_predictor/XGB_SCRICO_final_ygDNA.pkl')