import numpy as np

def calculate_fdr(score_array, label_array):
    # Arguments should be sorted such that best scores come first
    # This function is very similar to that provided by Percolator
    qvalues = []

    score_array = np.atleast_1d(score_array)
    label_array = np.atleast_1d(label_array)
    ntargets = np.cumsum(label_array == "target")
    ndecoys = np.arange(1, ntargets.shape[0] + 1) - ntargets
    qvalues = (ndecoys + 1)/ntargets
    for ind in np.arange(qvalues.shape[0])[::-1]:
         if ind - 1 >= 0 and score_array[ind] == score_array[ind-1]:
             qvalues[ind - 1] = qvalues[ind]

         if ind + 1 < qvalues.shape[0]:
             qvalues[ind] = min(qvalues[ind], qvalues[ind + 1])

    return qvalues
