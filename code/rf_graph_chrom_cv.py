import os
from optparse import OptionParser
from sklearn.model_selection import *
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.ensemble import RandomForestClassifier
import joblib
import networkx as nx

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option('-i', '--input_file',help='Samples for training the model')
parser.add_option('-o', '--output_path',help='Output path for storing the training results')
(opts, args) = parser.parse_args()
input_file = opts.input_file
output_path = opts.output_path
if not output_path.endswith('/'):
    output_path=output_path+'/'
if not os.path.exists(output_path):
    os.makedirs(output_path)

#################################Loading data#################################
print('1. Loading data')
sample = pd.read_csv(input_file)
f = []
f_suffix = ['motif_strand', 'motif_score', 'ctcf_signalValue', 'rad_signalValue', 'age']
f_suffix.extend([c + str(num) for num in range(19) for c in ['A', 'G', 'C', 'T']])
f_from = ['from_' + name for name in f_suffix]
f_to = ['to_' + name for name in f_suffix]
f.extend(f_from)
f.extend(f_to)
f.append('distance')
f.extend(['between_motif',
          'between_positive_motif',
          'between_negative_motif',
          'between_ctcf_score',
          'between_positive_ctcf_score',
          'between_negative_ctcf_score',
          'between_ctcf_signalValue',
          'between_positive_ctcf_signalValue',
          'between_negative_ctcf_signalValue'
          ])
feature = sample[f].values
label = sample.frequency.values >= 1
#################################训练模型#################################
print('2. Cross validation')

random_state = np.random.RandomState(0)
n_estimator = 500
rf_base = RandomForestClassifier(max_depth=16, max_features=45, n_estimators=n_estimator, random_state=random_state,
    min_samples_split=60,
    n_jobs=-1)
rf_graph = RandomForestClassifier(max_depth=16, max_features=45, n_estimators=n_estimator, random_state=random_state,
    min_samples_split=60,
    n_jobs=-1)

group_kfold = GroupKFold(n_splits=10)
group_kfold.get_n_splits(sample.index,groups=sample['chr'])
y_predict_prob = np.zeros(label.shape[0])


def get_graph_distance(graph, row):
    graph_restricted = nx.restricted_view(graph, [], [(row.from_motif_index, row.to_motif_index)])
    try:
        distance = nx.shortest_path_length(graph_restricted, source=row.from_motif_index, target=row.to_motif_index,
                                           weight='path_weight')
    except:
        distance = 100000
    return distance


for i, cv_index in enumerate(group_kfold.split(sample.index,groups=sample['chr'])):
    print('fold: %d/10' % (i))
    train_index, test_index = cv_index
    X_train, X_test = feature[train_index], feature[test_index]
    y_train, y_test = label[train_index], label[test_index]
    print('2.1 Training the basic model')
    rf_base.fit(X_train, y_train)
    sample_predict = rf_base.predict_proba(feature)[:, 1]
    sample_predict[sample_predict <= 0] = 0.000001
    sample_predict = -np.log(sample_predict)
    sample['path_weight'] = sample_predict
    print('2.2 Constructing the network and computing GCP')
    # conv_sample=sample[(sample.from_motif_strand==1)&(sample.to_motif_strand==0)]
    graph = nx.convert_matrix.from_pandas_edgelist(sample, source='from_motif_index', target='to_motif_index',
                                                   edge_attr=True)
    shorest = [get_graph_distance(graph, row) for index, row in sample.iterrows()]
    shorest = list(np.exp(-np.array(shorest)))
    new_f = f.copy()
    new_f.append('gcp')
    sample['gcp'] = shorest
    new_feature = sample[new_f].values
    X_train_new, X_test_new = new_feature[train_index], new_feature[test_index]
    print('2.3 Training the final model')
    rf_graph.fit(X_train_new, y_train)
    y_test_proba = rf_graph.predict_proba(X_test_new)[:, 1].ravel()
    y_predict_prob[test_index] = y_test_proba

np.save(output_path + "/cross_val_predict_prob.npy", y_predict_prob)

print('3. Printing the results')
precision, recall, thresholds = precision_recall_curve(label, y_predict_prob)
f1 = 2 * precision * recall / (precision + recall)
best_f1 = np.max(f1)
best_f1_index = np.argmax(f1)
au_pr = metrics.auc(recall, precision)

plt.figure(figsize=(7, 7))
plt.title('Precision Recall Curve')
plt.plot(recall, precision, 'b', label='AUPRC = %0.4f' % au_pr)
plt.plot(recall[best_f1_index], precision[best_f1_index], 'g^',
         label='Best F1 = (%0.4f,%0.4f)%0.4f' % (recall[best_f1_index], precision[best_f1_index], best_f1))
plt.vlines(recall[best_f1_index], 0, precision[best_f1_index], colors='g', linestyles='-.')
plt.hlines(precision[best_f1_index], 0, recall[best_f1_index], colors='g', linestyles='-.')
plt.legend(loc='lower left')
plt.plot([1, 0], [0, 1], 'r--')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('Precision score')
plt.xlabel('Recall score')
plt.savefig(output_path + 'Precision Recall Curve.svg', format='svg')

fpr, tpr, thresholds = roc_curve(label, y_predict_prob)
au_roc = auc(fpr, tpr)

plt.figure(figsize=(7, 7))
plt.title('Receiver Operating Characteristic')
plt.plot(fpr, tpr, 'b', label='AUROC = %0.4f' % au_roc)
print('AUROC = %0.4f' % au_roc)
plt.legend(loc='lower right')
plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('True Positive')
plt.xlabel('False Positive')
plt.savefig(output_path + 'Receiver Operating Characteristic.svg', format='svg')

#################################使用全部数据#################################
print('4. Using all the data to train the model')
rf_base.fit(feature, label)
sample_predict = rf_base.predict_proba(feature)[:, 1]
sample_predict[sample_predict <= 0] = 0.000001
sample_predict = -np.log(sample_predict)
sample['path_weight'] = sample_predict
graph = nx.convert_matrix.from_pandas_edgelist(sample, source='from_motif_index', target='to_motif_index',
                                               edge_attr=True)
shorest = [get_graph_distance(graph, row) for index, row in sample.iterrows()]
shorest = list(np.exp(-np.array(shorest)))
new_f = f.copy()
new_f.append('gcp')
sample['gcp'] = shorest
new_feature = sample[new_f].values
rf_graph.fit(new_feature, label)

joblib.dump(rf_base, output_path + 'rf_base.model')
joblib.dump(rf_graph, output_path + 'rf_graph.model')
