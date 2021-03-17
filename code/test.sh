python generate_pairs.py -c ../data/HeLa-S3/CTCF_peak.bed\
                         -m ../data/HeLa-S3/fimo.csv\
                         -p ../data/HeLa-S3/hela_ctcf.interactions.intra.bedpe\
                         -o ../data/HeLa-S3/output
                         
python generate_features.py -c ../data/HeLa-S3/CTCF_peak.bed\
                            -m ../data/HeLa-S3/fimo.csv\
                            -p ../data/HeLa-S3/hela_ctcf.interactions.intra.bedpe\
                            -o ../data/HeLa-S3/output\
                            -r ../data/HeLa-S3/rad21.narrowPeak\
                            -a ../data/HeLa-S3/CTCF_age.bed
python rf_graph.py -i ../data/HeLa-S3/output/sample.csv -o ../data/HeLa-S3/output/model


python generate_pairs.py -c ../data/GM12878/CTCF_peak.bed\
                         -m ../data/GM12878/fimo.csv\
                         -p ../data/GM12878/gm12878_ctcf.interactions.intra.bedpe\
                         -o ../data/GM12878/output
                         
python generate_features.py -c ../data/GM12878/CTCF_peak.bed\
                            -m ../data/GM12878/fimo.csv\
                            -p ../data/GM12878/gm12878_ctcf.interactions.intra.bedpe\
                            -o ../data/GM12878/output\
                            -r ../data/GM12878/rad21.narrowPeak\
                            -a ../data/GM12878/CTCF_age.bed
python rf_graph.py -i ../data/GM12878/output/sample.csv -o ../data/GM12878/output/model



python generate_pairs.py -c ../data/K562/CTCF_peak.bed\
                         -m ../data/K562/fimo.csv\
                         -p ../data/K562/k562_ctcf.interactions.intra.bedpe\
                         -o ../data/K562/output
                         
python generate_features.py -c ../data/K562/CTCF_peak.bed\
                            -m ../data/K562/fimo.csv\
                            -p ../data/K562/k562_ctcf.interactions.intra.bedpe\
                            -o ../data/K562/output\
                            -r ../data/K562/rad21.narrowPeak\
                            -a ../data/K562/CTCF_age.bed
python rf_graph.py -i ../data/K562/output/sample.csv -o ../data/K562/output/model



python generate_pairs.py -c ../data/MCF-7/CTCF_peak.bed\
                         -m ../data/MCF-7/fimo.csv\
                         -p ../data/MCF-7/mcf_ctcf.interactions.intra.bedpe\
                         -o ../data/MCF-7/output
                         
python generate_features.py -c ../data/MCF-7/CTCF_peak.bed\
                            -m ../data/MCF-7/fimo.csv\
                            -p ../data/MCF-7/mcf_ctcf.interactions.intra.bedpe\
                            -o ../data/MCF-7/output\
                            -r ../data/MCF-7/rad21.narrowPeak\
                            -a ../data/MCF-7/CTCF_age.bed
python rf_graph.py -i ../data/MCF-7/output/sample.csv -o ../data/MCF-7/output/model




python rf_graph_test.py -o ../cross_cell_val\
                        -m ../data/GM12878/output/model\
                        -s ../data/MCF-7/output/sample.csv\
                        -M GM12878\
                        -S MCF-7