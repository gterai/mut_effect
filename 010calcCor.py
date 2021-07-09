# coding: shift_jis

import sys
import re
import numpy as np
args = sys.argv

if len(args) != 3:
    print ("usage:", args[0], "[aaindex1] [mut2titer.txt]")
    exit(0);

def main():
    id2desc = {}
    aaindex = read_aaindex(args[1], id2desc)
    mut2titer = {}
    id2cor = {}
    try:
        with open(args[2], 'r') as fh:
            for line in fh.readlines():
                line = line.strip()
                (mut, titer) = line.split()

                mut2titer[mut] = float(titer)
    except FileNotFoundError as e:
        print("Fine not found", e)
    except Exception as e: 
        print(e)
    
    for id in aaindex.keys():
        x = np.empty(0)
        y = np.empty(0)
        for mut in mut2titer.keys():
            diff = 0;
            if(mut != "WT"):
                m = re.match(r'^([A-Z])(\d+)([A-Z])$', mut)
                (aa_before, pos, aa_after) = m.groups() 
                if(aaindex[id][aa_before] >= 1e100 or aaindex[id][aa_before] >= 1e100):
                    continue
                
                diff = aaindex[id][aa_before] - aaindex[id][aa_after]
            
            x = np.append(x,diff)
            y = np.append(y,mut2titer[mut])
            
        #print(x,y)
        #print (np.corrcoef(x,y)[0][1])
        cor = np.corrcoef(x,y)[0][1]
        #print (id, x, y, file=sys.stderr)
        #continue
        if(np.isnan(cor)):
            #skip
            continue
        else:
            id2cor[id] = cor
            #print (cor)
        #    for id in aaindex:
        #        print(id)
    #exit(1)
    ranked = sorted(id2cor.items(), key = lambda x:abs(x[1]), reverse=True)
    for item in ranked:
        print(item[0], item[1], id2desc[item[0]], sep="\t")
    
def read_aaindex(file, id2d):
    dic = {}
    id = "";
    try:
        with open(file, 'r') as fh:
            #for line in fh.readlines():
            line = fh.readline()
            while(line):
                line = line.strip()
                #print(line)
                
                if(line[0] == 'H' and line[1] == " "): 
                    a = line.split(" ")
                    id = a[1]
                    #print(id)
                elif(line[0] == 'D' and line[1] == " "): 
                    #a = line.split(" ")
                    id2d[id] = line[2:]
                    #print(id)
                elif(line[0] == "I" and line[1] == " "):
                    tmp_aa_list = line.split()
                    aa_pref = []
                    aa_suff = []
                    for aa_pair in tmp_aa_list[1:]:
                        if(aa_pair[1] != '/'):
                            print("Unexpected amino acid pair ({0})".format(aa_pair), file=sys.stderr)
                            sys.exit(1)
                        else:
                            aa_pref.append(aa_pair[0])
                            aa_suff.append(aa_pair[2])
                    aa_list = aa_pref + aa_suff
                    #print(aa_list)
                    
                    line1 = fh.readline()
                    line2 = fh.readline()
                    line1 = line1.strip()
                    line2 = line2.strip()
                    #print(line1)
                    #print(line2)
                    
                    val1 = line1.split()
                    val2 = line2.split()
                    
                    val = val1 + val2
                    
                    #print (aa_list)
                    #print (id, val)
                    dic[id] = {}
                    for i in range(len(aa_list)):
                        if(val[i] == "NA"):
                            val[i] = 1e100
                        dic[id][aa_list[i]] = float(val[i])

                    #(aa_list, val)
                    #print(dic[id])
                    
                line = fh.readline()

    except FileNotFoundError as e:
        print("Fine not found", e)
    except Exception as e: 
        print(e)

    return dic

main()
