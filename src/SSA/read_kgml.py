"""functions to read kgml files"""
#%%
import xml.etree.ElementTree as ET
import pandas as pd
from . import ReactionNetwork
import requests
import time
#%%
#return a dataframe of reactions
def get_reactionData(filepath):
    with open(filepath, 'r') as f:
        kgml = f.read()

    root = ET.fromstring(kgml)

    dict_id = IDdict(filepath)

    reactions = [child for child in root if child.tag == 'reaction']

    reaction_list = []
    for reac in reactions:
        reacName=reac.attrib['name']
        sub_id_list = [sub.attrib['id'] for sub in reac.findall('substrate')]
        pro_id_list = [pro.attrib['id'] for pro in reac.findall('product')]
        sub_name_list = [dict_id[sub_id] for sub_id in sub_id_list]
        pro_name_list = [dict_id[pro_id] for pro_id in pro_id_list]
        reaction_list.append([reacName,sub_name_list, pro_name_list])

        direction = reac.attrib['type']
        if direction=='reversible':
            reaction_list.append([reacName+'_2',pro_name_list,sub_name_list])

    return reaction_list

def make_network_from_kgml(filepath,info=True):
    reaction_list=get_reactionData(filepath)
    if len(reaction_list)==0:
        raise Exception('empty reaction_list')
    network=ReactionNetwork.ReactionNetwork(reaction_list,info=info)
    df_reaction=network.to_df()

    if df_reaction.index.duplicated().sum()>0:
        # same name, not identical reactions
        # print('reaction name is not identical,', df_reaction.index.duplicated().sum())
        raise Exception('reaction name is not identical')
    else:
        return network

def IDdict(filepath):
    dict_id=dict()
    with open(filepath, 'r') as f:
        kgml = f.read()
    root = ET.fromstring(kgml)

    entries = [child for child in root if child.tag == 'entry']
    for entry in entries:
        dict_id[entry.attrib['id']]=entry.attrib['name']
    return dict_id

#convert cpdname to id in kgml
#cpd names need to be identified by IDs in kgml file
def cpdname2id(cpdname, filepath):
    sglcpd_id_list=[]
    dict_id=IDdict(filepath)
    for elem in dict_id.items():
        if cpdname in elem[1].split(' '):
            return elem[0]
    
    raise Exception(f'{cpdname} does not exist in {filepath}')


def complete_reaction(reaction):
    # get reaction data from kegg api
    # reaction=['rn:R01070_2',['cpd:C00118'],['cpd:C05378']]

    reaction_id=reaction[0][:9]
    lhs,rhs=reaction[1],reaction[2]
    lhs=[cpd[-6:] for cpd in lhs]
    rhs=[cpd[-6:] for cpd in rhs]

    api='http://rest.kegg.jp/get/'+reaction_id
    form_data=requests.get(api).text
    time.sleep(1)
    for line in form_data.splitlines():
        if 'EQUATION' in line:
            form_line=line
            break
    else:
        print('api data is not correct')
        1/0

    kegg_lhs,kegg_rhs=form_line[12:].split('<=>') #remove 'EQUATION' part

    kegg_lhs_cpds=kegg_lhs.replace(' ','').split('+')
    kegg_rhs_cpds=kegg_rhs.replace(' ','').split('+')
    for cpd in kegg_lhs_cpds:
        if cpd[0]!='C':
            # coefficient of the reaction
            n = int(cpd[:cpd.find('C')])
            cpdname = cpd[cpd.find('C'):]
            kegg_lhs_cpds.remove(cpd)
            kegg_lhs_cpds.extend([cpdname]*n)
    for cpd in kegg_rhs_cpds:
        if cpd[0]!='C':
            # coefficient of the reaction
            n = int(cpd[:cpd.find('C')])
            cpdname = cpd[cpd.find('C'):]
            kegg_rhs_cpds.remove(cpd)
            kegg_rhs_cpds.extend([cpdname]*n)
    
    # direction...  True: => False: <+
    if (not set(lhs).isdisjoint(set(kegg_lhs_cpds))) and (not set(rhs).isdisjoint(set(kegg_rhs_cpds))):
        direction=True
    elif (not set(lhs).isdisjoint(set(kegg_rhs_cpds))) and (not set(rhs).isdisjoint(set(kegg_lhs_cpds))):
        direction=False
    else:
        print('direction cannot be determined')
        1/0
    if direction:
        reaction_complete=[reaction[0],kegg_lhs_cpds,kegg_rhs_cpds]
    else:
        reaction_complete=[reaction[0],kegg_rhs_cpds,kegg_lhs_cpds]

    return reaction_complete
#%%
if __name__ == '__main__':
    org = 'hsa'
    pathway_id=f'{org}01200'
    path_pathway = f'{DATADIR}/network/pathways_kgml/{org}/{pathway_id}.xml'