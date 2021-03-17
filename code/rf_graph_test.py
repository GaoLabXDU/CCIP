import os
import numpy as np
from sklearn import metrics
import pandas as pd
import joblib
from sklearn.metrics import roc_curve,precision_recall_curve,auc
import networkx as nx
from optparse import OptionParser



def get_graph_distance(graph, row):
    graph_restricted = nx.restricted_view(graph, [], [(row.from_motif_index, row.to_motif_index)])
    try:
        distance = nx.shortest_path_length(graph_restricted, source=row.from_motif_index, target=row.to_motif_index,
                                           weight='path_weight')
    except:
        distance = 100000
    return distance
def run(model_path,sample_file,output_path,model_cell,sample_cell):
    model_base_path = model_path+'/rf_base.model'
    model_graph_path = model_path+'/rf_graph.model'

    data_path = sample_file
    print('Reading in features...')
    sample = pd.read_csv(data_path)
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

    rf_base = joblib.load(model_base_path)
    rf_graph=joblib.load(model_graph_path)
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
    y_predict_prob=rf_graph.predict_proba(new_feature)


    np.save(output_path + "%s_%s_ccip_prob.npy" % (model_cell, sample_cell), y_predict_prob[:, 1])

    precision, recall, thresholds = precision_recall_curve(label, y_predict_prob[:, 1])
    au_pr = metrics.auc(recall, precision)
    fpr, tpr, thresholds = roc_curve(label, y_predict_prob[:, 1])
    au_roc = auc(fpr, tpr)
    print("auroc for the model: %0.4f"%(au_roc))

parser = OptionParser()
parser.add_option('-o', '--output_path',help='Output path for storing the testing results')
parser.add_option('-m', '--model_path',help='Model file for predicting')
parser.add_option('-s', '--sample_file',help='Sample file for predicting')
parser.add_option('-M', '--model_cell',help='Model cell for predicting')
parser.add_option('-S', '--sample_cell',help='Sample cell for predicting')
(opts, args) = parser.parse_args()
output_path = opts.output_path
model_path = opts.model_path
sample_file = opts.sample_file
model_cell = opts.model_cell
sample_cell = opts.sample_cell

if not output_path.endswith('/'):
    output_path=output_path+'/'
if not os.path.exists(output_path):
    os.makedirs(output_path)

run(model_path,sample_file,output_path,model_cell,sample_cell)
