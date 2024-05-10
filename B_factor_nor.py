from cif_read import cif_read
import numpy as np

def b_normal(filename_1, chain):

    [RNA_name_chain, N_nucleotide, nucleotide_index, nucleotide_name, nucleotide_B_factor, nucleotide_coordinate] = cif_read(filename_1, chain)

    B_median = np.median(nucleotide_B_factor)
    B_mean = np.mean(nucleotide_B_factor)
    B_std = np.std(nucleotide_B_factor)

    l = len(nucleotide_B_factor)
    outlier_index = []
    x = []
    for i in range(l):
        x.append(abs(nucleotide_B_factor[i] - B_median))
    x_median = np.median(x)

    for j in range(l):
        if (0.6745*(nucleotide_B_factor[j] - B_median)/x_median) >= 3.5:
            outlier_index.append(j)
        else:
            continue

    for k in range(len(outlier_index)):
        del_index = outlier_index[k] - k
        nucleotide_B_factor.pop(del_index)
        nucleotide_name.pop(del_index)

    nucleotide_name1 = nucleotide_name
    nucleotide_B_factor1 = nucleotide_B_factor

    B_mean1 = np.mean(nucleotide_B_factor1)
    B_std1 = np.std(nucleotide_B_factor1)

    for i in range(len(nucleotide_B_factor1)):
        nucleotide_B_factor1[i] = (nucleotide_B_factor1[i] - B_mean1)/B_std1
        nucleotide_B_factor1[i] = "%.2f" % nucleotide_B_factor1[i]
    return B_mean, B_std, B_mean1, B_std1, nucleotide_B_factor1, nucleotide_name1, outlier_index
