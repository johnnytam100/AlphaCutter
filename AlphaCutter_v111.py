import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from argparse import ArgumentParser
import glob
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from scipy import stats
import pcmap
import itertools
import networkx as nx
import re
import os

debug = 0

if debug == 0:

    parser = ArgumentParser()
    parser.add_argument('--loop_min', type=int)
    parser.add_argument('--helix_min', type=int)
    parser.add_argument('--domain_min', type=int)
    parser.add_argument('--fragment_min', type=int)
    parser.add_argument('--pLDDT_min', type=float)
    parser.add_argument('--local_contact_range', type=int)
    parser.add_argument('--domain_out', action='store_true')
    parser.add_argument('--single_out', action='store_true')
    args = parser.parse_args()

    loop_min_len = args.loop_min
    helix_min_len = args.helix_min
    domain_min_len = args.domain_min
    if all([loop_min_len==None,
             helix_min_len==None]):
        fragment_min_len = 0
        local_contact_range = 0
    else:
        fragment_min_len = args.fragment_min
        local_contact_range = args.local_contact_range
    pLDDT_min = args.pLDDT_min
    domain_out = args.domain_out
    single_out = args.single_out
    all_pdb = '*.pdb'

if debug == 1:

    loop_min_len = None
    helix_min_len = None
    domain_min_len = 0
    fragment_min_len = 5
    pLDDT_min = 0
    local_contact_range = 5
    domain_out = True
    single_out = True
    all_pdb = '*_v2.pdb'
    #all_pdb = './break_pdb/*v2.pdb'
    #all_pdb = 'AF-A5DP36-F1-model_v2.pdb'        # break, res443-444
    #all_pdb = 'AF-A0A248AFK6-F1-model_v2.pdb'    # break, res1900-1901
    #all_pdb = 'AF-A5GIN1-F1-model_v2.pdb'
    #all_pdb = 'ma-*-????.pdb'
    #all_pdb = "AF-A0A0E3SVE7-F1-model_v2.pdb"
    #all_pdb = "AF-A0A0G2KTI4-F1-model_v2.pdb"
    #all_pdb = "AF-A0A0A2JW91-F1-model_v2.pdb"


def idx2range(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


pdb_df_groupby = np.nan

def get_helix_loop_list(dssp_secondary_struc, pdb_df_groupby, contact_list, min_len):
    # loop list
    loop_list = []

    for idx, rows in pdb_df_groupby.iterrows():

        # collect consective secondary structure in a list
        secondary_struc_list_consec = []

        start_idx = idx

        for num in range(min_len):
            
            current_idx = start_idx + num

            if current_idx < pdb_df_groupby.index.max():

                current_secondary_struc = pdb_df_groupby.at[current_idx, "secondary_struc"]

                secondary_struc_list_consec.append(current_secondary_struc)

        # check if all elements in secondary_struc_list_consec are loop
        secondary_struc_list_consec_boo = [i == dssp_secondary_struc for i in secondary_struc_list_consec]

        # print(idx, secondary_struc_list_consec)

        # append residue range of loops
        if all([len(secondary_struc_list_consec_boo) > 0,
                all(secondary_struc_list_consec_boo),
                current_idx <= pdb_df_groupby.index.max()]):
              
                loop_list.append([start_idx, current_idx])

    # Remove loop with any non-local contact
    loop_list_contact_boo = []

    for loop in loop_list:

        all_no_contact_boo = all(contact_list[loop[0]: loop[1]])

        loop_list_contact_boo.append(all_no_contact_boo)

    loop_list = [i for (i, v) in zip(loop_list, loop_list_contact_boo) if v]    # remove elements according to boolean list

    return loop_list


def merge_list(input_list, pdb_df_groupby):

    # merge loop list (i.e. merge overlapping residue range of loops to get range of long loops)
    output_list = []

    idx = 0

    # As long as there is next loop in loop list...
    while idx + 1 <= len(input_list):

        # define start position (old, to keep) and end position (new, to be updated)
        start_old = input_list[idx][0]
        end_new = input_list[idx][-1]

        # if reaching the last index, append last loop directly
        if idx + 1 == len(input_list):

            output_list.append([start_old, end_new])

            break

        else:
            # As long as the end position of current loop is in the range of the next loop...
            while end_new in range(input_list[idx+1][0], input_list[idx+1][-1]):

                # update end position
                end_new = input_list[idx+1][-1]

                # index + 1, to consider the next loop for merging
                idx += 1      

                # avoid list out of range error, break this while loop if idx + 1 == len(input_list)
                if idx + 1 == len(input_list):

                    break

            # append the new range of the long loop from start_old to end_new
            output_list.append([start_old, end_new])

        # index + 1, to consider the next loop for next merging cycle
        idx += 1

        # avoid list out of range error, break this while loop if idx + 1 == len(input_list)
        if idx == len(input_list):

            break

        # avoid repeatingly adding loop lists at the C-terimal, break this while loop if end_new already reached the max residue number
        if end_new == pdb_df_groupby.index.max():

            break

    return output_list



# three2one_dict
three2one_dict = {
    "GLY":"G",
    "ALA":"A",
    "VAL":"V",
    "LEU":"L",
    "ILE":"I",
    "PRO":"P",
    "PHE":"F",
    "TRP":"W",
    "SER":"S",
    "THR":"T",
    "MET":"M",
    "CYS":"C",
    "TYR":"Y",
    "GLU":"E",
    "ASP":"D",
    "GLN":"Q",
    "ASN":"N",
    "HIS":"H",
    "ARG":"R",
    "LYS":"K"
}


# summary df
summary_df_all = pd.DataFrame(columns=["PDB","chainID","type","start","end","pLDDT_mean"])


# Collect sequences from all *.pdb in current directory
filepath_list = []

for filepath in glob.iglob(all_pdb):
  filepath_list.append(filepath)

filepath_list.sort()

for pdb in filepath_list:

    print('Processing... ' + pdb)

    # Open .pdb to dataframe
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb)
    pdb_df = ppdb.df['ATOM']

    # individual chains df to dict
    pdb_df_dict = { chain_id : pdb_df.loc[pdb_df["chain_id"]==chain_id] for chain_id in pdb_df["chain_id"].unique()}
    pdb_df_dict_groupby = {}

    is_looping = True

    # separate, save chains to individual pdb
    for chain_id, chain_df in pdb_df_dict.items():

        # separate, save
        ppdb.df['ATOM'] = chain_df
        pdb_sep = pdb.replace(".pdb","") + "_chain" + chain_id + ".pdb"
        ppdb.to_pdb(path=pdb_sep, records=['ATOM'], gz=False, append_newline=True)

        # DSSP
        p = PDBParser()
        structure = p.get_structure(pdb_sep, pdb_sep)
        model = structure[0]
        dssp = DSSP(model, pdb_sep, dssp='mkdssp')

        secondary_struc_list = []

        chain_df["secondary_struc"] = np.nan

        for key in list(dssp.keys()):

            # check if res_num_dssp matches residue in pdb_df
            res_num_dssp = dssp[key][0]
            res_dssp = dssp[key][1]
            secondary_struc = dssp[key][2]

            chain_df_res_three = chain_df.loc[chain_df["residue_number"] == res_num_dssp]["residue_name"].iloc[0]
            chain_df_res_one = three2one_dict[chain_df_res_three]

            if chain_df_res_one != res_dssp:

                print("Mismatch between DSSP and PDB residue detected! Res_num =",
                      res_num_dssp,
                      "DSSP_res:",
                      res_dssp,
                      "PDB_res:",
                      chain_df_res_one) 
                
                is_looping = False
                
                break
              
            else:

                chain_df.loc[chain_df['residue_number'] == res_num_dssp, 'secondary_struc'] = secondary_struc

        if not is_looping:

            break

        # groupby residue
        pdb_df_dict_groupby[chain_id] = chain_df.groupby('residue_number').agg(lambda x: stats.mode(x)[0][0])


    if not is_looping:

        # remove chains .pdb
        for chain_id in pdb_df_dict.keys():
            os.remove(pdb.replace(".pdb","") + "_chain" + chain_id + ".pdb")

        continue


    # contact map df
    contact_map_dict = {}

    for chain_id1, chain_df1 in pdb_df_dict_groupby.items():

        res_list1 = [ chain_id1 + str(res_num1) for res_num1 in chain_df1.index ]

        res_list2_concat = []

        for chain_id2, chain_df2 in pdb_df_dict_groupby.items():

            res_list2 = set(chain_df2.index)

            res_list2 = [chain_id2 + str(res_num2) for res_num2 in res_list2]

            res_list2_concat += res_list2
      
        contact_map_dict[chain_id1] = pd.DataFrame(0, index=res_list1, columns=res_list2_concat)

    # contact map dfs to dict
    for chain_id1, chain_df1 in pdb_df_dict_groupby.items():

        pdb1 = pdb.replace(".pdb","") + "_chain" + chain_id1 + ".pdb"

        for chain_id2, chain_df2 in pdb_df_dict_groupby.items():

            pdb2 = pdb.replace(".pdb","") + "_chain" + chain_id2 + ".pdb"

            pcmap_dict = pcmap.contactMap(pdb1, pdb2)

            # assign contact map
            for i in pcmap_dict['data']:

                res_num = (chain_id1 + i['root']['resID']).replace(" ","")

                res_contact_num_list = []
                
                for j in i['partners']:

                    res_contact_num_list.append((chain_id2 + j['resID']).replace(" ",""))

                for res_contact_num in res_contact_num_list:

                    contact_map_dict[chain_id1].at[res_num, res_contact_num] = 1

    # remove chains .pdb
    for chain_id in pdb_df_dict.keys():
        os.remove(pdb.replace(".pdb","") + "_chain" + chain_id + ".pdb")

    # all-to-all contact map
    contact_map_list = [contact_map for chain_id, contact_map in contact_map_dict.items()]
    contact_map_df = pd.concat(contact_map_list)

    # contact list dict
    contact_list_dict = {}

    for (chain_id, chain_df), (chain_id, contact_map) in zip(pdb_df_dict_groupby.items(), contact_map_dict.items()):

        contact_list = []

        for idx, rows in chain_df.iterrows():

            # drop list (the index list for local contact -> plus and minus local_contact_range)
            if idx + local_contact_range > chain_df.index.max():

                drop_list = [chain_id + str(i) for i in range(idx - local_contact_range, idx)]

            elif idx - local_contact_range < chain_df.index.min():

                drop_list = [chain_id + str(i) for i in range(idx, idx + local_contact_range)]

            else:

                drop_list = [chain_id + str(i) for i in range(idx - local_contact_range, idx + local_contact_range)]

            # print(drop_list)
            
            # get contact map for 
            res_num = chain_id + str(idx)
            contact_map_current = contact_map.loc[res_num]

            # drop  local contact plus and minus local_contact_range
            contact_map_current = contact_map_current.drop(drop_list)

            # check if completely no contact for other position
            current_contact = (contact_map_current == 0).all()

            # append boolean
            contact_list.append(current_contact)

          # print(idx, helix_contact_list_consec)

        contact_list_dict[chain_id] = contact_list


    ############
    ''' Loop '''
    ############

    if loop_min_len == None:

        loop_list_dict = { chain_id : [] for chain_id, chain_df in pdb_df_dict_groupby.items() }
        loop_list_merged_dict = loop_list_dict

    else:

        loop_list_dict = { chain_id : get_helix_loop_list("-", chain_df, contact_list_dict[chain_id], loop_min_len) for chain_id, chain_df in pdb_df_dict_groupby.items() }
        loop_list_merged_dict = { chain_id : merge_list(get_helix_loop_list("-", chain_df, contact_list_dict[chain_id], loop_min_len), chain_df) for chain_id, chain_df in pdb_df_dict_groupby.items() }

    #############
    ''' Helix '''
    #############

    if helix_min_len == None:

        helix_list_dict = { chain_id : [] for chain_id, chain_df in pdb_df_dict_groupby.items() }
        helix_list_merged_dict = helix_list_dict

    else:

        helix_list_dict = { chain_id : get_helix_loop_list("H", chain_df, contact_list_dict[chain_id], helix_min_len) for chain_id, chain_df in pdb_df_dict_groupby.items() }
        helix_list_merged_dict = { chain_id : merge_list(get_helix_loop_list("H", chain_df, contact_list_dict[chain_id], helix_min_len), chain_df) for chain_id, chain_df in pdb_df_dict_groupby.items() }

    ##############
    ''' Domain '''
    ##############

    loop_helix_list_merged_dict = { chain_id : sorted(loop_list + helix_list) for (chain_id, loop_list), (chain_id, helix_list) in zip(loop_list_merged_dict.items(), helix_list_merged_dict.items()) }

    domain_list_dict = {}

    for (chain_id, chain_df), (chain_id, loop_helix_list_merged) in zip(pdb_df_dict_groupby.items(), loop_helix_list_merged_dict.items()):

        # domain list (i.e. get the range of domain to keep)
        chain_df["loop_helix_boo"] = False

        for i in loop_helix_list_merged:

            for idx in range(i[0], i[1]+1):

                chain_df.at[idx, "loop_helix_boo"] = True

        domain_index_list = chain_df.loc[chain_df["loop_helix_boo"] == False].index.tolist()
        domain_list = [list(i) for i in list(idx2range(domain_index_list))]

        # Avoid removing loop/helix with length < loop_min_len at the C-terminal
        if len(domain_list) > 0:

            if all([loop_min_len == None,
                    helix_min_len == None]):

                pass

            elif all([loop_min_len == None,
                      helix_min_len != None]):
                              
                if chain_df.index.max() - domain_list[-1][-1] < helix_min_len:
                  
                    domain_list[-1][-1] = chain_df.index.max()
            
            elif all([loop_min_len != None,
                      helix_min_len == None]):
                              
                if chain_df.index.max() - domain_list[-1][-1] < loop_min_len:
                  
                    domain_list[-1][-1] = chain_df.index.max()

            else:
              
                if any([chain_df.index.max() - domain_list[-1][-1] < loop_min_len,
                        chain_df.index.max() - domain_list[-1][-1] < helix_min_len]):

                    domain_list[-1][-1] = chain_df.index.max()

        domain_list_dict[chain_id] = domain_list

    #########################################
    ''' Merge contacting domain fragments '''
    #########################################

    domain_dict_dict = {}
    domain_key_contact_df_list = []

    for chain_id, domain_list in domain_list_dict.items():

        # domain pairwise contact df
        domain_dict = {}
        count = 1
        for domain in domain_list:
            domain_dict["domain_" + chain_id + str(count)] = domain
            count +=1


        domain_dict_dict[chain_id] = domain_dict
        

    domain_dict_keys = [k for chain_id, domain_dict in domain_dict_dict.items() for k, v in domain_dict.items()]

    domain_key_contact_df = pd.DataFrame(0, index=domain_dict_keys, columns=domain_dict_keys)


    # assign contact
    for chain_id1, domain_dict1 in domain_dict_dict.items():

        for chain_id2, domain_dict2 in domain_dict_dict.items():

            for k1, v1 in domain_dict1.items():

                for k2, v2 in domain_dict2.items():

                    if v1 != v2:
                        
                        # domain1 and domain2 contact map

                        domain1_domain2_contact_map_df = contact_map_df.loc[chain_id1 + str(v1[0]) : chain_id1 + str(v1[1]) , chain_id2 + str(v2[0]) : chain_id2 + str(v2[1])]

                        # if domain1 makes contacts with domain2
                        if not (domain1_domain2_contact_map_df == 0).all().all():

                            # assign
                            domain_key_contact_df.at[k1, k2] = 1


    # get isolated network
    G = nx.DiGraph(domain_key_contact_df.values)
    G = G.to_undirected()
    components= nx.connected_components(G)
    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    domain_list_final = []

    for x in S:
        idx_list = list(x.nodes)
        isolated_network = domain_key_contact_df.iloc[idx_list,:].index.tolist()
        domain_list_final.append(isolated_network)


    # domain to range
    domain_dict_dict_format = { domain: [ chain_id + str(range[0]), chain_id + str(range[1])] for chain_id, domain_dict in domain_dict_dict.items() for domain, range in domain_dict.items()}
    domain_list_final_dict = []

    for merged_domain in domain_list_final:

        merged_domain_dict = {}

        for k in pdb_df_dict.keys():

            merged_domain_dict[k] = []

        for domain in merged_domain:

            domain_range = domain_dict_dict_format[domain]

            start = re.findall('\d+|\D+',domain_range[0])
            end = re.findall('\d+|\D+',domain_range[1]) 

            merged_domain_dict[start[0]].append(start[1])
            merged_domain_dict[end[0]].append(end[1])

        domain_list_final_dict.append(merged_domain_dict)


    ##################
    ''' Cut Domain '''
    ##################

    domain_df_list_single = []

    for merged_domain in domain_list_final_dict:

        domain_df_list = []

        for chain_id, domain in merged_domain.items():

            if len(domain) > 0:

                for idx in range(0, len(domain), 2):

                    domain_df_full = pdb_df_dict[chain_id]
                    domain_df_cut = domain_df_full.loc[ ( domain_df_full["residue_number"] >= int(domain[idx]) ) & ( domain_df_full["residue_number"] <= int(domain[idx+1])) ]

                    # domain cut length
                    domain_df_cut_chain_res_unique = domain_df_cut[['chain_id', 'residue_number']].drop_duplicates()
                    domain_cut_len = domain_df_cut_chain_res_unique.shape[0]

                    # fragment length cutoff
                    if domain_cut_len >= fragment_min_len:

                        # append
                        domain_df_list.append(domain_df_cut)

        
        if len(domain_df_list) > 0:

            # merge
            domain_df = pd.concat(domain_df_list)

            # domain length
            domain_df_chain_res_unique = domain_df[['chain_id', 'residue_number']].drop_duplicates()
            domain_len = domain_df_chain_res_unique.shape[0]

            # domain pLDDT mean
            domain_pLDDT_mean = domain_df["b_factor"].mean()

            # filter
            if all([domain_pLDDT_mean >= pLDDT_min,
                    domain_len >= domain_min_len]) :
                    
                # load
                ppdb.read_pdb(pdb)

                # assign domain df to original ppdb object
                domain_df_out = domain_df.drop(["secondary_struc"], axis=1)
                ppdb.df['ATOM'] = domain_df_out

                # append to single PDB output
                domain_df_list_single.append(domain_df_out)

                # range for naming
                output_range_list = []
                for chain_id in domain_df_chain_res_unique['chain_id'].unique():
                    start = domain_df_chain_res_unique.loc[domain_df_chain_res_unique['chain_id'] == chain_id ]['residue_number'].min()
                    end = domain_df_chain_res_unique.loc[domain_df_chain_res_unique['chain_id'] == chain_id ]['residue_number'].max()
                    output_range_list.append(chain_id + str(start) + "-" + str(end))

                # save
                if domain_out == True:

                    output_range = "_".join(output_range_list)
                    ppdb.to_pdb(path=pdb.replace('.pdb','') + '_AFCT-OUT_domain_' + output_range + '.pdb', 
                            records=None, 
                            gz=False, 
                            append_newline=True)

    
    # Single PDB output (i.e. merge all non-contacting domains)
    if all([single_out == True,
            len(domain_df_list_single) > 0]):

        # merge single PDB output
        domain_df_single = pd.concat(domain_df_list_single)

        # domain length
        domain_df_single_chain_res_unique = domain_df_single[['chain_id', 'residue_number']].drop_duplicates()

        # range for naming
        output_range_list = []
        for chain_id in domain_df_single_chain_res_unique['chain_id'].unique():
            start = domain_df_single_chain_res_unique.loc[domain_df_single_chain_res_unique['chain_id'] == chain_id ]['residue_number'].min()
            end = domain_df_single_chain_res_unique.loc[domain_df_single_chain_res_unique['chain_id'] == chain_id ]['residue_number'].max()
            output_range_list.append(chain_id + str(start) + "-" + str(end))

        # sort atom
        domain_df_single = domain_df_single.sort_values("atom_number")

        # load
        ppdb.read_pdb(pdb)

        # assign domain df to original ppdb object
        ppdb.df['ATOM'] = domain_df_single

        # save
        output_range = "_".join(output_range_list)
        ppdb.to_pdb(path=pdb.replace('.pdb','') + '_AFCT-OUT_single_' + output_range + '.pdb', 
                records=None, 
                gz=False, 
                append_newline=True)


    ##################
    ''' Summary df '''
    ##################

    # start, end
    loop_start_list = []
    loop_end_list = []
    loop_chainid_list = []

    helix_start_list = []
    helix_end_list = []
    helix_chainid_list = []

    domain_start_list = []
    domain_end_list = []
    domain_chainid_list = []

    for chain_id, loop_list_merged in loop_list_merged_dict.items():

        if len(loop_list_merged) > 0:

            for i in loop_list_merged:
                loop_start_list.append(i[0])
                loop_end_list.append(i[1])
                loop_chainid_list.append(chain_id)

    for chain_id, helix_list_merged in helix_list_merged_dict.items():

        if len(helix_list_merged) > 0:

            for i in helix_list_merged:
                helix_start_list.append(i[0])
                helix_end_list.append(i[1])
                helix_chainid_list.append(chain_id)

    for chain_id, domain_list in domain_list_dict.items():

        if len(domain_list) > 0:

            for i in domain_list:
                domain_start_list.append(i[0])
                domain_end_list.append(i[1])
                domain_chainid_list.append(chain_id)

    # pLDDT mean
    loop_pLDDT_mean_list = []
    helix_pLDDT_mean_list = []
    domain_pLDDT_mean_list = []

    for chain_id, loop_list_merged in loop_list_merged_dict.items():

        if len(loop_list_merged) > 0:

            for loop in loop_list_merged:
                loop_df = pdb_df.loc[ ( pdb_df["residue_number"] >= loop[0] ) & ( pdb_df["residue_number"] <= loop[1] )]
                loop_pLDDT_mean_list.append(loop_df["b_factor"].mean())

    for chain_id, helix_list_merged in helix_list_merged_dict.items():

        if len(helix_list_merged) > 0:

            for helix in helix_list_merged:
                helix_df = pdb_df.loc[ ( pdb_df["residue_number"] >= helix[0] ) & ( pdb_df["residue_number"] <= helix[1] )]
                helix_pLDDT_mean_list.append(helix_df["b_factor"].mean())

    for chain_id, domain_list in domain_list_dict.items():

        if len(domain_list) > 0:

            for domain in domain_list:
                domain_df = pdb_df.loc[ ( pdb_df["residue_number"] >= domain[0] ) & ( pdb_df["residue_number"] <= domain[1] )]
                domain_pLDDT_mean_list.append(domain_df["b_factor"].mean())


    loop_start_end_df = pd.DataFrame({"PDB": pdb.replace('.pdb',''),
                            "chainID": loop_chainid_list,
                            "type": "loop",
                            "start": loop_start_list,
                            "end": loop_end_list,
                            "pLDDT_mean": loop_pLDDT_mean_list})

    helix_start_end_df = pd.DataFrame({"PDB": pdb.replace('.pdb',''),
                            "chainID": helix_chainid_list,
                            "type": "helix",
                            "start": helix_start_list,
                            "end": helix_end_list,
                            "pLDDT_mean": helix_pLDDT_mean_list})

    domain_start_end_df = pd.DataFrame({"PDB": pdb.replace('.pdb',''),
                            "chainID": domain_chainid_list,
                            "type": "domain",
                            "start": domain_start_list,
                            "end": domain_end_list,
                            "pLDDT_mean": domain_pLDDT_mean_list})

    summary_df = pd.concat([domain_start_end_df,
                            helix_start_end_df,
                            loop_start_end_df])

    # append
    summary_df = summary_df.sort_values("start")
    summary_df_all = pd.concat([summary_df_all,summary_df])

# assign
summary_df_all.reset_index(drop=True, inplace=True)
summary_df_all.insert(5, "len", summary_df_all["end"] - summary_df_all["start"] + 1)

# add seq
seq_list = []
summary_df_all.sort_values(['PDB', 'chainID', 'start'], ascending=[True, True, True], inplace=True)

for pdb in summary_df_all["PDB"].unique():

    summary_df_all_pdb = summary_df_all.loc[summary_df_all["PDB"] == pdb]

    ppdb.read_pdb(pdb + ".pdb")
    pdb_df = ppdb.df['ATOM']

    for chain_id in summary_df_all_pdb["chainID"].unique():

        summary_df_all_pdb_one = summary_df_all_pdb.loc[summary_df_all_pdb["chainID"]==chain_id]

        for idx, row in summary_df_all_pdb_one.iterrows():

            start_idx = pdb_df.loc[(pdb_df["chain_id"]==chain_id) & (pdb_df["residue_number"]==row["start"])].index[0]
            end_idx = pdb_df.loc[(pdb_df["chain_id"]==chain_id) & (pdb_df["residue_number"]==row["end"])].index[0]
                   
            seq_df = ppdb.amino3to1()
            seq = seq_df.loc[seq_df["chain_id"]==chain_id].loc[int(start_idx):int(end_idx)+1]['residue_name'].tolist()
            seq = ''.join(seq)
            seq_list.append(seq)

summary_df_all["seq"] = seq_list

summary_df_all.reset_index(drop=True, inplace=True)

# save fasta
fasta = []

for index, row in summary_df_all.iterrows():
  fasta.append(">"+row["PDB"]+"_"+row["chainID"]+str(int(row["start"]))+"-"+str(int(row["end"]))+"_"+row["type"]+"\n"+row["seq"]+"\n")

with open("AFCT-OUT_output.fasta", "w") as output:
  output.write(''.join(fasta))

# save summary
summary_df_all.to_csv("AFCT-OUT_summary.csv")
