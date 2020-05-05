import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

p_prefix = sys.argv[1]
f_mtx = sys.argv[2]
f_doublet = sys.argv[3]
f_hist = sys.argv[4]
f_tSNE = sys.argv[5]

counts_matrix = scipy.io.mmread(f_mtx).T.tocsc()

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.07)


doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

print(doublet_scores)
print(np.count_nonzero(predicted_doublets))
print(np.count_nonzero(doublet_scores>0.2))

plt.hist(doublet_scores, normed=True, log=True)
plt.savefig(f_hist, dpi=300)

predicted_doublets = scrub.call_doublets(threshold=0.2)

np.savetxt(f_doublet, doublet_scores)

#scrub.plot_histogram()

#print('Running UMAP...')
#scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

# # Uncomment to run tSNE - slow
print('Running tSNE...')
scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))

# # Uncomment to run force layout - slow
# print('Running ForceAtlas2...')
# scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5. n_iter=1000))
    
print('Done.')

#scrub.plot_embedding('UMAP', order_points=True)
scrub.plot_embedding('tSNE', order_points=True)
# scrub.plot_embedding('FA', order_points=True)

plt.savefig(f_tSNE, dpi=300)
