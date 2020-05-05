import sys
import numpy as np
import argparse

def add_to_clonotype(clonotype_to_CDR3, clonotype_id, TCRab, CDR3):
    try:
        clonotype_to_CDR3[clonotype_id]
    except KeyError:
        clonotype_to_CDR3[clonotype_id] = {}
    try:
        clonotype_to_CDR3[clonotype_id][(TCRab, CDR3)]
    except KeyError:
        clonotype_to_CDR3[clonotype_id][(TCRab, CDR3)] = 1
    return clonotype_to_CDR3

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-a", default=None, help="clonotypeAB|clonotypeGD|CDR3b|CDR3g")
    parser.add_argument("-c", default=None, help="filtered clonotype annotation csv")
    parser.add_argument("-b", default=None, help="CBC whitelist")
    args = parser.parse_args()

    p_app = args.a
    f_clonotypes = args.c
    f_whitelist = args.b

    clonotype_to_CDR3 = {}
    all_cdr3 = {}
    with open(f_clonotypes, "r") as f:
        for no, line in enumerate(f):
            if no == 0:
                continue
            
            row = line.strip("\r\n").split(",")
            cbc = row[0][:len(row[0])-2]  # row[0]
            TCRab = row[5]
            full_length = row[10]
            productive = row[11]
            CDR3 = row[12]
            clonotype_id = row[16]
            clonotype_consensus = row[17]
            
            if "consensus" not in clonotype_consensus:
                continue
            
            clonotype_to_CDR3 = add_to_clonotype(clonotype_to_CDR3, clonotype_id, TCRab, CDR3)

            if p_app == "CDR3b":
                if TCRab == "TRB" and full_length == "True" and productive == "True":
                    try:
                        all_cdr3[CDR3].append(cbc)
                    except KeyError:
                        all_cdr3[CDR3] = [cbc]
            elif p_app == "CDR3g":
                if TCRab == "TRG" and full_length == "True" and productive == "True":
                    try:
                        all_cdr3[CDR3].append(cbc)
                    except KeyError:
                        all_cdr3[CDR3] = [cbc]
            elif p_app == "clonotypeAB":
                if TCRab in ("TRA", "TRB") and full_length == "True" and productive == "True":
                    try:
                        if cbc not in all_cdr3[clonotype_id]:
                            all_cdr3[clonotype_id].append(cbc)
                    except KeyError:
                        all_cdr3[clonotype_id] = [cbc]
            elif p_app == "clonotypeGD":
                if TCRab in ("TRG", "TRD") and full_length == "True" and productive == "True":
                    try:
                        if cbc not in all_cdr3[clonotype_id]:
                            all_cdr3[clonotype_id].append(cbc)
                    except KeyError:
                        all_cdr3[clonotype_id] = [cbc]

    
    index_cbc = np.loadtxt(f_whitelist, dtype=str).tolist()
    #print ",".join([""] + [ cbc[:-2] for cbc in index_cbc ])
    print ",".join([""] + [ cbc for cbc in index_cbc ])
    
    index_cdr3 = []
    for k, v in sorted(all_cdr3.items(), key=lambda v: len(v[1]), reverse=True):
        index_cdr3.append(k)

    uniq_cdr3 = {}
    
    for cdr3 in index_cdr3:
        data = np.zeros(len(index_cbc), dtype=int)
        try:
            cbcs = all_cdr3[cdr3]
            for cbc in cbcs:
                try:
                    data[index_cbc.index(cbc)] += 1
                except ValueError:
                    pass
        except KeyError:
            continue
        
        if p_app in ("CDR3b", "CDR3g"):
            print ",".join(map(str, [cdr3] + data.tolist()))
        elif p_app in ("clonotypeAB", "clonotypeGD"):
            
            name = [ "%s:%s" % (k, v) for k, v in sorted(clonotype_to_CDR3[cdr3].keys()) ]
            name = "|".join(name)

            try:
                uniq_cdr3[name]
            except KeyError:
                uniq_cdr3[name] = 1
                print ",".join(map(str, [name] + data.tolist()))
    
            
