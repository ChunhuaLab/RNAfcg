from cif_read import cif_read
from B_factor_nor import b_normal
import math
import numpy as np
import networkx as nx


def letter(t):
    global a
    a = 0
    if t == "A":
        a = "0 0 0 1 "
    elif t == "G":
        a = "0 0 1 0 "
    elif t == "C":
        a = "0 1 0 0 "
    elif t == "U":
        a = "1 0 0 0 "
    return a

names = locals()
def FRI(RNAID_chain):
    [RNA_name_chain, N_nucleotide, nucleotide_index, nucleotide_name, nucleotide_B_factor, nucleotide_coordinate] = \
        cif_read(RNAID_chain[0:4].lower() + ".cif", RNAID_chain[5:7].strip())
    [B_mean, B_std, B_mean1, B_std1, nucleotide_B_factor1, nucleotide_name1, outlier_index] = \
        b_normal(RNAID_chain[0:4].lower() + ".cif", RNAID_chain[5:7].strip())

    n1, k1 = 6, 1       # Parameters of the kernel function
    for n in range(n1 - 3, n1 + 4):
        names["flexible_" + str(n)] = []
        names["flexible_guiyihua_" + str(n)] = []
        num = len(nucleotide_coordinate)
        C_M = np.zeros([num, num])
        for i1 in range(num):
            for j1 in range(num):
                if i1 == j1:
                    continue
                elif i1 != j1:
                    dis = math.sqrt((nucleotide_coordinate[i1][0] - nucleotide_coordinate[j1][0]) ** 2
                                    + (nucleotide_coordinate[i1][1] - nucleotide_coordinate[j1][1]) ** 2
                                    + (nucleotide_coordinate[i1][2] - nucleotide_coordinate[j1][2]) ** 2)
                    aa = math.exp(-(dis / n) ** k1)
                    C_M[i1, j1] = aa
            C_M[i1, i1] = sum(C_M[i1, :])

            if C_M[i1, i1] == 0:
                C_M[i1, i1] = math.exp(-710)

        ''' Flexibility index '''
        for i2 in range(num):
            if 1 / C_M[i2, i2] == float('inf'):
                names["flexible_" + str(n)].append(1 / math.exp(-400))
            else:
                names["flexible_" + str(n)].append(1 / (C_M[i2, i2]))
        for i3 in range(len(outlier_index)):
            del_index = outlier_index[i3] - i3
            names["flexible_" + str(n)].pop(del_index)

        RI_mean = np.mean(names["flexible_" + str(n)])
        RI_std = np.std(names["flexible_" + str(n)])
        for i4 in range(len(names["flexible_" + str(n)])):
            names["flexible_guiyihua_" + str(n)].append(str("%.4f" % ((names["flexible_" + str(n)][i4] - RI_mean) / RI_std)))

    return flexible_guiyihua_3, flexible_guiyihua_4, flexible_guiyihua_5, flexible_guiyihua_6, flexible_guiyihua_7, flexible_guiyihua_8, flexible_guiyihua_9

def complex_network(RNAID_chain, CN_cutoff):
    [RNA_name_chain, N_nucleotide, nucleotide_index, nucleotide_name, nucleotide_B_factor, nucleotide_coordinate] = \
        cif_read(RNAID_chain[0:4].lower() + ".cif", RNAID_chain[5:7].strip())
    [B_mean, B_std, B_mean1, B_std1, nucleotide_B_factor1, nucleotide_name1, outlier_index] = \
        b_normal(RNAID_chain[0:4].lower() + ".cif", RNAID_chain[5:7].strip())
    num = len(nucleotide_coordinate)
    com_M = np.zeros([num, num])
    for i1 in range(num):
        for j1 in range(num):
            dis = math.sqrt((nucleotide_coordinate[i1][0] - nucleotide_coordinate[j1][0]) ** 2
                            + (nucleotide_coordinate[i1][1] - nucleotide_coordinate[j1][1]) ** 2
                            + (nucleotide_coordinate[i1][2] - nucleotide_coordinate[j1][2]) ** 2)
            if i1 == j1:
                continue
            elif dis <= CN_cutoff:
                com_M[i1, j1] = 1
            elif dis > CN_cutoff:
                continue

    G = nx.Graph()
    G = nx.from_numpy_matrix(np.matrix(com_M))
    Graph_information = nx.info(G)

    cluster_coefficient = nx.clustering(G)
    degrees_cen = nx.degree_centrality(G)
    closenesses_cen = nx.closeness_centrality(G)
    betweennesses = nx.betweenness_centrality(G)

    du_centrality = []
    jujixishu = []
    jiejin_centrality = []
    zhongjie_centrality = []

    for i2 in range(num):
        du_centrality.append(str("%.4f" % (degrees_cen[i2])))
        jujixishu.append(str("%.4f" % (cluster_coefficient[i2])))
        jiejin_centrality.append(str("%.4f" % (closenesses_cen[i2])))
        zhongjie_centrality.append(str("%.4f" % (betweennesses[i2])))
    for i3 in range(len(outlier_index)):
        del_index = outlier_index[i3] - i3
        du_centrality.pop(del_index)
        jujixishu.pop(del_index)
        jiejin_centrality.pop(del_index)
        zhongjie_centrality.pop(del_index)
    return du_centrality, jujixishu, jiejin_centrality, zhongjie_centrality

def ASA(RNAID_chain):
    ASA_list = []
    rsa_file = open(RNAID_chain[0:7].strip() + ".rsa", "r")
    [B_mean, B_std, B_mean1, B_std1, nucleotide_B_factor_nor1, nucleotide_name1, outlier_index] = b_normal(
        RNAID_chain[0:4].lower() + ".cif", RNAID_chain[5:7].strip())
    rsa_file = rsa_file.readlines()
    n = 0
    m = 0
    cc = []
    for j in rsa_file:
        if j[0:3] == "RES":
            m += 1
            if n in outlier_index:
                n += 1
                continue
            elif j[5:7].strip() in ["A", "G", "C", "U"]:
                cc.append(float(j[15:22]))
            n += 1
    average = np.mean(cc)
    std_asa = np.std(cc)
    for j1 in cc:
        aa = "%.4f" % ((j1 - average) / std_asa)
        ASA_list.append(aa + "   ")
    return ASA_list

def ss(RNAID_chain):
    ss_list = []
    cif_file = open(RNAID_chain[0:4] + ".cif", "r")
    cif_file = cif_file.readlines()
    basepair_file = open("SS_" + RNAID_chain + ".csv", "r")
    baspa = basepair_file.readlines()
    sequence = []
    res_ID = []
    ss = []

    for ato in cif_file:
        ato = ato.split()
        if ato == "\n":
            continue
        elif len(ato) <= 3:
            continue
        elif ato[0] == 'END':
            break
        elif ((ato[0] == 'ATOM') and (ato[18] == RNAID_chain[5:7].strip()) and (ato[19] == "\"C1\'\"") and (ato[17] in ['A', 'U', 'C', 'G'])):
            # elif (ato[0:4] == 'ATOM') and (ato[20:22].strip() == RNAID_chain[5:7].strip()) and (ato[72:74].strip() == RNA_chain_CIF2) and (ato[13:16] == "C1\'") and (ato[17:20] in ['A', 'U', 'C', 'G']):
            sequence.append(ato[17])
            res_ID.append(ato[16])
    pair1_ID = []
    pair2_ID = []
    W_C = []
    for bpa in baspa:
        pair1_ID.append(bpa.split(",")[3])
        pair2_ID.append(bpa.split(",")[6])
        W_C.append(bpa.split(",")[9].strip("\n"))

    for res_n in res_ID:
        if res_n in pair1_ID:
            res_n_pairindex1 = pair1_ID.index(res_n)
            if W_C[res_n_pairindex1] == "Y":
                ss.append("(")
            elif W_C[res_n_pairindex1] == "N":
                ss.append("[")
        elif res_n in pair2_ID:
            res_n_pairindex2 = pair2_ID.index(res_n)
            if W_C[res_n_pairindex2] == "Y":
                ss.append(")")
            elif W_C[res_n_pairindex2] == "N":
                ss.append("]")
        elif (res_n not in pair1_ID) and (res_n not in pair2_ID):
            ss.append(".")

    [B_mean, B_std, B_mean1, B_std1, nucleotide_B_factor_nor1, nucleotide_name1, outlier_index] = b_normal(RNAID_chain[0:4] + ".cif", RNAID_chain[5:7].strip())

    n = 0
    for j in ss:
        if n in outlier_index:
            n += 1
            continue
        elif j == "(" or j == ")":
            ss_list.append("0 0 1 ")      # Watson-Creek paired
        elif j == "[" or j == "]":
            ss_list.append("0 1 0 ")      # Non-watson Creek paired
        elif j == ".":
            ss_list.append("1 0 0 ")      # unpaired
        n += 1
    return ss_list

def chain_lon(cha_l):
    global long1
    long1 = 0
    if cha_l <= 50:
        long1 = "0 0 0 1 "
    elif 50 < cha_l <= 100:
        long1 = "0 0 1 0 "
    elif 100 < cha_l <= 200:
        long1 = "0 1 0 0 "
    elif 200 < cha_l:
        long1 = "1 0 0 0 "
    return long1


def size_RNA(RNAID_chain):
    [RNA_name_chain, N_nucleotide, nucleotide_index, nucleotide_name, nucleotide_B_factor, nucleotide_coordinate] = \
        cif_read(RNAID_chain[0:4].lower() + ".cif", RNAID_chain[5:7].strip())
    dis_lar = 0
    num = len(nucleotide_coordinate)
    for i in range(num):
        for j in range(num):
            if i == j:
                continue
            elif i != j:
                dis = math.sqrt((nucleotide_coordinate[i][0] - nucleotide_coordinate[j][0]) ** 2
                                + (nucleotide_coordinate[i][1] - nucleotide_coordinate[j][1]) ** 2
                                + (nucleotide_coordinate[i][2] - nucleotide_coordinate[j][2]) ** 2)
                if dis_lar <= dis:
                    dis_lar = dis
    return dis_lar


def rg(RNAID_chain):
    [RNA_name_chain, N_nucleotide, nucleotide_index, nucleotide_name, nucleotide_B_factor, nucleotide_coordinate] = \
        cif_read(RNAID_chain[0:4].lower() + ".cif", RNAID_chain[5:7].strip())
    f = open(RNAID_chain[0:4] + ".cif", "r")
    rna_chain = RNAID_chain[5:7].strip()
    ''' Computed center of mass '''
    nucleotide_name = []
    nucleotide_coordinate = []
    for ato in f:
        ato = ato.split()
        if ato == "\n":
            continue
        elif len(ato) <= 3:
            continue
        elif ato[0] == 'END':
            break
        elif ((ato[0] == 'ATOM') and (ato[3] == "\"C1\'\"") and (ato[5] in ['A', 'U', 'C', 'G']) and ato[18] == rna_chain):
            coordinate_x = float(ato[10])
            coordinate_y = float(ato[11])
            coordinate_z = float(ato[12])
            nucleotide_coordinate.append([coordinate_x, coordinate_y, coordinate_z])
            nucleotide_name.append(ato[5])
    f.close()
    M = 0
    for rg_i in nucleotide_name:
        if rg_i == "A":
            M = M + 347  # A Molecular weight
        elif rg_i == "G":
            M = M + 363  # G Molecular weight
        elif rg_i == "C":
            M = M + 323  # C Molecular weight
        elif rg_i == "U":
            M = M + 324  # U Molecular weight
    Xc = 0
    Yc = 0
    Zc = 0
    for rg_i1 in range(len(nucleotide_coordinate)):
        if nucleotide_name[rg_i1] == "A":
            Xc = Xc + 347 * float(nucleotide_coordinate[rg_i1][0])
            Yc = Yc + 347 * float(nucleotide_coordinate[rg_i1][1])
            Zc = Zc + 347 * float(nucleotide_coordinate[rg_i1][2])
        elif nucleotide_name[rg_i1] == "G":
            Xc = Xc + 363 * float(nucleotide_coordinate[rg_i1][0])
            Yc = Yc + 363 * float(nucleotide_coordinate[rg_i1][1])
            Zc = Zc + 363 * float(nucleotide_coordinate[rg_i1][2])
        elif nucleotide_name[rg_i1] == "C":
            Xc = Xc + 323 * float(nucleotide_coordinate[rg_i1][0])
            Yc = Yc + 323 * float(nucleotide_coordinate[rg_i1][1])
            Zc = Zc + 323 * float(nucleotide_coordinate[rg_i1][2])
        elif nucleotide_name[rg_i1] == "U":
            Xc = Xc + 324 * float(nucleotide_coordinate[rg_i1][0])
            Yc = Yc + 324 * float(nucleotide_coordinate[rg_i1][1])
            Zc = Zc + 324 * float(nucleotide_coordinate[rg_i1][2])
    Xc = Xc / M
    Yc = Yc / M
    Zc = Zc / M

    ''' Calculate the moment of inertia '''
    Jz = 0
    for rg_j in range(len(nucleotide_coordinate)):
        r = math.sqrt((nucleotide_coordinate[rg_j][0] - Xc) ** 2
                      + (nucleotide_coordinate[rg_j][1] - Yc) ** 2
                      + (nucleotide_coordinate[rg_j][2] - Zc) ** 2)
        if nucleotide_name[rg_j] == "A":
            Jz = Jz + 347 * r ** 2
        elif nucleotide_name[rg_j] == "G":
            Jz = Jz + 363 * r ** 2
        elif nucleotide_name[rg_j] == "C":
            Jz = Jz + 323 * r ** 2
        elif nucleotide_name[rg_j] == "U":
            Jz = Jz + 324 * r ** 2

    ''' radius of gyration '''
    RG = "%.4f" % (float(math.sqrt(Jz / M)/N_nucleotide))
    return RG


if __name__ == '__main__':

    rna_id_chain = '4LNT_RA'   # Change to the target PDB file
    [RNA_name_chain, N_nucleotide, nucleotide_index, nucleotide_name, nucleotide_B_factor, nucleotide_coordinate] = \
        cif_read(rna_id_chain[0:4].lower() + ".cif", rna_id_chain[5:7].strip())
    [B_mean, B_std, B_mean1, B_std1, nucleotide_B_factor_nor1, nucleotide_name1, outlier_index] = b_normal(
        rna_id_chain[0:4].lower() + ".cif", rna_id_chain[5:7].strip())

    nucleotide = nucleotide_name1
    label = nucleotide_B_factor_nor1

    allfea_win_siz = open(rna_id_chain + "_all_feature.csv", "w")
    size_RNA1 = 0
    one = []

    ''' one-hot encoding '''
    for t1 in nucleotide:
        one.append(letter(t1))

    ''' FRI '''
    [FRI_list3, FRI_list4, FRI_list5, FRI_list6, FRI_list7, FRI_list8, FRI_list9] = FRI(rna_id_chain)

    ''' Topological centrality '''
    [du_centrality_list, jujixishu_list, jiejin_centrality_list, zhongjie_centrality_list] = complex_network(rna_id_chain, CN_cutoff=23)

    ''' ASA '''
    ASA1_list = ASA(rna_id_chain)

    ''' SS '''
    ss1_list = ss(rna_id_chain)

    ''' Chain length '''
    num_res = len(nucleotide)
    long1 = chain_lon(num_res)

    ''' RNA size '''
    size_RNA1 = str("%.4f" % float((size_RNA(rna_id_chain)-171.8095)/70.82437))

    ''' radius of gyration '''
    RG = rg(rna_id_chain)


    allfea_win_siz.writelines('base1' + "," + 'base2' + "," + 'base3' + "," + 'base4' + "," + 'FIn=3' + "," + 'FIn=4' + "," + 'FIn=5' + "," + 'FIn=6' + "," + 'FIn=7' + "," + 'FIn=8' + "," + 'FIn=9' + "," + 'degree' + "," + 'cluster' + "," + 'closeness' + "," + 'betweennesses' + "," + 'ASA' + "," + 'SS1' + "," + 'SS2' + "," + 'SS3' + "," + 'chain_long1' + "," + 'chain_long2' + "," + 'chain_long3' + "," + 'chain_long4' + "," + 'size' + "," + 'RG' + "," + 'labels' + "\n")
    for l1 in range(len(nucleotide)):
        allfea_win_siz.writelines(sin1 + "," for sin1 in one[l1].split())
        allfea_win_siz.writelines(FRI_list3[l1])  # FRI
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(FRI_list4[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(FRI_list5[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(FRI_list6[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(FRI_list7[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(FRI_list8[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(FRI_list9[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(du_centrality_list[l1])    # Topological centrality
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(jujixishu_list[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(jiejin_centrality_list[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(zhongjie_centrality_list[l1])
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(ASA1_list[l1])                # ASA
        allfea_win_siz.write(",")
        allfea_win_siz.writelines(sin2 + "," for sin2 in ss1_list[l1].split())                 # SS
        for sin3 in long1.split():                 # Chain length
            allfea_win_siz.write(sin3)
            allfea_win_siz.write(",")
        allfea_win_siz.write(str(size_RNA1))                    # RNA size
        allfea_win_siz.write(",")
        allfea_win_siz.write(str(RG))                           # radius of gyration
        allfea_win_siz.write(",")
        allfea_win_siz.write(label[l1])                # label
        allfea_win_siz.write("\n")

