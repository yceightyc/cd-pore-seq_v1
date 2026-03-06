import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, auc, confusion_matrix, accuracy_score, classification_report
import matplotlib.font_manager as fm
import matplotlib as mpl
import seaborn as sns
import time
from sklearn.decomposition import PCA
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from sklearn.utils import resample
from sklearn.multiclass import OneVsRestClassifier
import shap

# Set plot style
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
    'axes.edgecolor': 'black',
    'ytick.left': False  
})

colour_dict = {
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

# Import data
data = pd.read_csv('./data/model_input_final.tsv', delimiter = '\t')
data['Type'] = data['Type'].replace({'Ctl': 0, 'Ribo' : 1})
centralR_design = pd.read_csv('./input/contig_centralR.tsv', delimiter = '\t')
data_withRinfo = pd.merge(data, centralR_design, on = 'contig', how = 'left')

conditions = [
    data_withRinfo['Type'] == 0,
    data_withRinfo['cen_R'] == 'A',
    data_withRinfo['cen_R'] == 'C',
    data_withRinfo['cen_R'] == 'G',
    data_withRinfo['cen_R'] == 'T'
]

choices = [0, 1, 2, 3, 4]

# Annotate central rN type
data_withRinfo['rType'] = np.select(conditions, choices, default=np.nan)

# Get feature variables, remove offset positions R-18 and R-15 as base contexts are varied, leading to possible overfitting
columns_to_drop = ['Type', 'read_name', 'contig', 'dwell_time_R-18', 'dwell_time_R-15','ampl_mean_R-18', 'ampl_mean_R-15', 'model_deviation_R-18', 'model_deviation_R-15','cen_R', 'rType']
X = data_withRinfo.drop(columns = columns_to_drop)
X.columns = [str(col) for col in X.columns]

# Impute median value to the NA
imputer = SimpleImputer(strategy='median')
X = imputer.fit_transform(X)

# Get response variable to y
y = data_withRinfo['rType']

# Run PCA on all data (after imputation)
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)
var_ratio = pca.explained_variance_ratio_ * 100

# Create PCA plot
colours = [colour_dict[label] for label in y]
g = sns.JointGrid(height=5, ratio=2, space=0.1)

# --- Main Scatter Plot --- 
scatter = g.ax_joint.scatter(X_pca[:, 0], X_pca[:, 1], c=colours, alpha=0.2, edgecolor=colours, rasterized = True) 
for val in sorted(set(y)): 
    plt.scatter([], [], c=colour_dict[val], label=label_dict[val]) 
plt.legend(title='Type', loc='best', fontsize=9, markerscale=0.7, title_fontsize=10) 

# --- Marginal KDE Plots (Split by Group) --- 
for group in [0, 1, 2, 3, 4]: 
    mask = (y == group) 
    sns.kdeplot(
        x=X_pca[mask, 0], 
        ax=g.ax_marg_x, 
        color=colour_dict[group], 
        fill=True, 
        alpha=0.5, 
        label=label_dict[group] ) 

for group in [0, 1, 2, 3, 4]: 
    mask = (y == group) 
    sns.kdeplot( 
        y=X_pca[mask, 1], 
        ax=g.ax_marg_y, 
        color=colour_dict[group], 
        fill=True, alpha=0.5, 
        label=label_dict[group] ) 

# --- Labels--- 
    g.ax_joint.set_xlabel(f'PC1 ({var_ratio[0]:.2f}%)') 
    g.ax_joint.set_ylabel(f'PC2 ({var_ratio[1]:.2f}%)') 

plt.tight_layout()
plt.savefig("./data/PCA_plot.pdf")

# Train-test-split
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size = 0.8, random_state = 42, stratify=y)
X_train.shape, X_test.shape

# Train RF
best_rf = RandomForestClassifier(random_state=30, n_jobs=-1, max_depth=20, min_samples_leaf=5,
                                 n_estimators=1000, oob_score=True)
best_rf.fit(X_train, y_train)

# Get list of classes
classes = best_rf.classes_
n_classes = len(classes)

# Binarize y_test for ROC
y_test_bin = label_binarize(y_test, classes=classes)

# Predict probabilities
y_score = best_rf.predict_proba(X_test)
y_pred = best_rf.predict(X_test)

# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()

for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

plt.figure(figsize=(7, 7))
for i in range(n_classes):
    label_name = label_dict.get(classes[i], f"Class {classes[i]}")
    plt.plot(fpr[i], tpr[i], lw=2, color=colour_dict[i], label=f'{label_name} (AUC = {roc_auc[i]:.2f})')

# Generate ROC plot
plt.plot([0, 1], [0, 1], 'k--') 
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve for Multiclass Random Forest')
plt.legend(loc='lower right')
plt.savefig("./data/multi-class_ROC.pdf")

# Generate confusion matrix showing recalls
custom_color = "#DB7176"
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_blue", ["#ffffff", custom_color])
label_names = [label_dict[i] for i in range(len(label_dict))]
cm_test = confusion_matrix(y_test, y_pred, normalize='true')
cm_test_normalized = cm_test.astype(float) / cm_test.sum(axis=1, keepdims=True)
plt.figure(figsize=(6, 5))
sns.heatmap(cm_test_normalized, annot=True, fmt=".2f", cmap=custom_cmap,
            xticklabels=label_names, yticklabels=label_names,
            cbar_kws={'label': 'Proportion'})
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Recalls in Hold-out Test-set")
plt.tight_layout()
plt.savefig("./data/CM_plot.pdf")

# Feature Correlation
features = ["ampl_mean_4th", "ampl_mean_5th", "ampl_mean_6th", "ampl_mean_7th", "ampl_mean_8th",   
            "model_deviation_4th", "model_deviation_5th", "model_deviation_6th", "model_deviation_7th", "model_deviation_8th",
            "dwell_time_4th", "dwell_time_5th", "dwell_time_6th", "dwell_time_7th", "dwell_time_8th"
            ]
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
plt.savefig("./data/Feature_Correlation.pdf")

# Feature Importance
custom_labels = {
    0: "dN",
    1: "rA",
    2: "rC",
    3: "rG",
    4: "rU"
}
class_names_list = [custom_labels[i] for i in range(len(custom_labels))]

custom_colors = ['#DB7176', '#698A9C', '#A7C7E7', '#878580', '#C0B6A4']

explainer = shap.TreeExplainer(best_rf)
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
plt.savefig("./data/Feature_Importance.pdf")