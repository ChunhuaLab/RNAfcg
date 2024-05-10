def cif_read(filename, rna_chain):
    cif = open(filename, 'r')
    nucleotide_index = []
    nucleotide_name = []
    nucleotide_coordinate = []
    nucleotide_B_factor = []
    RNA = cif.readline()
    RNA_name = RNA[5:9]

    for ato in cif:
        ato = ato.split()
        if ato == "\n":
            continue
        elif len(ato) <= 3:
            continue
        elif ato[0] == 'END':
            break
        elif ((ato[0] == 'ATOM') and (ato[19] == "\"C1\'\"")  and (ato[17] in ['A', 'U', 'C', 'G']) and ato[18] == rna_chain):
            coordinate_x = float(ato[10])
            coordinate_y = float(ato[11])
            coordinate_z = float(ato[12])
            nucleotide_coordinate.append([coordinate_x, coordinate_y, coordinate_z])
            nucleotide_B_factor.append(float(ato[14]))
            nucleotide_index.append(int(ato[16]))
            nucleotide_name.append(ato[17])
    cif.close()
    N_nucleotide = len(nucleotide_index)
    RNA_name_chain = RNA_name + "_" + rna_chain
    return RNA_name_chain, N_nucleotide, nucleotide_index, nucleotide_name, nucleotide_B_factor, nucleotide_coordinate


