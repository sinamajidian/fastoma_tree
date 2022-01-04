#!/usr/bin/env python
# coding: utf-8

# In[129]:


#!/usr/bin/python3
import os
import numpy as np
from sys import argv
from datetime import datetime
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ete3
#import dendropy
from dendropy import Tree
datasets_address="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/"
# oma_database_address = datasets_address + "OmaServer.h5"
omaID_address = datasets_address+"oma-species.txt"
#bird6ID_address = datasets_address+"info.tsv"
uniprot_table_address = datasets_address+"proteomes_taxonomy_Mammalia_40674_refyes.tab"


# In[130]:


def read_omaID_file(omaID_address):
    
    omaID_file = open(omaID_address,'r')

    taxonID_omaID={}
    omaID_taxonID={}

    # omaID_scienceFull={}
    # scienceFull_omaID={}
    scientific_taxonID={}
    taxonID_scientific={}

    #  !!! limitation  ignoring strains and isolate
    for line in omaID_file:
        line_strip = line.strip()
        if line_strip.startswith('#'):
            pass
            #header_lines_list.append(line_strip)
        else:
            line_parts = line_strip.split('\t')

            omaID = line_parts[0]
            taxonID = line_parts[1]
            taxonID_omaID[taxonID] = omaID
            omaID_taxonID[omaID] = taxonID

            scientific = line_parts[2]
            taxonID_scientific[taxonID] = scientific
            scientific_taxonID[scientific] = taxonID

    omaID_file.close()

    print("-- The map for OMA taxonID of",len(taxonID_omaID),"records have read.") 
    
    return (taxonID_omaID,omaID_taxonID,taxonID_scientific,scientific_taxonID)


def read_taxonID_uniprot(uniprot_table_address):

    taxonID_scientific = {}
    scientific_taxonID = {}
    common_taxonID = {}
    taxonID_uniprot = {}

    uniprot_file = open(uniprot_table_address,'r')
    for line in uniprot_file:
        line_strip = line.strip()
        if line_strip.startswith('Pro'):
            pass
        else:
            line_parts = line_strip.split('\t')
            uniprot = line_parts[0]
            taxonID = line_parts[2]

            taxonID_uniprot[uniprot]=taxonID
            
            line_parts_1= line_parts[1].split("(")
            if len(line_parts_1) > 1:

                scientific = line_parts_1[0]
                common = line_parts_1[1].strip()[:-1]
            else:
                scientific = line_parts_1[0]
                common = line_parts_1[0]
            
            scientific_taxonID[taxonID] = scientific
            taxonID_scientific[scientific] = taxonID 

            common_taxonID[taxonID] = common

    print("-- The map for uniprot taxonID of",len(scientific_taxonID),"records have read.") 
    return (scientific_taxonID,taxonID_scientific,common_taxonID,taxonID_uniprot)



def read_taxonID_uniprot_temp(uniprot_table_address):

    taxonID_scientific = {}
    scientific_taxonID = {}
    common_taxonID = {}
    taxonID_uniprot = {}
    
    taxonID_scientificUnder={}
    
    
    uniprot_file = open(uniprot_table_address,'r')
    for line in uniprot_file:
        line_strip = line.strip()
        if line_strip.startswith('Pro'):
            pass
        else:
            line_parts = line_strip.split('\t')
            uniprot = line_parts[0]
            taxonID = line_parts[2]

            taxonID_uniprot[uniprot]=taxonID
            
            line_parts_1= line_parts[1].split("(")
            if len(line_parts_1) > 1:

                scientific = line_parts_1[0]
                common = line_parts_1[1].strip()[:-1]
            else:
                scientific = line_parts_1[0]
                common = line_parts_1[0]
            
            scientific_taxonID[taxonID] = scientific
            taxonID_scientific[scientific] = taxonID 

            common_taxonID[taxonID] = common
            scientificUnder='_'.join(scientific.split(" ")[:2])
            taxonID_scientificUnder[scientificUnder]=taxonID

    print("-- The map for uniprot taxonID of",len(scientific_taxonID),"records have read.") 
    return (taxonID_scientificUnder)


# In[131]:


taxonID_scientificUnder = read_taxonID_uniprot_temp(uniprot_table_address)


(scientific_taxonID,taxonID_scientific,common_taxonID,taxonID_uniprot)=read_taxonID_uniprot(uniprot_table_address)

(taxonID_omaID,omaID_taxonID,taxonID_scientific,scientific_taxonID)= read_omaID_file(omaID_address)


# In[ ]:





# In[ ]:





# # Reading FastOMA tree - collapsed

# In[148]:


project_folder_root="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/mammal/mml27/"
fast_tree_address=project_folder_root+"fastoma_core/_100_msa_concatanated.txt.contree_collapsed_94_0.001"

print(fast_tree_address)
print(round(os.path.getsize(fast_tree_address)/1000),"kb")


# In[149]:


fast_tree= ete3.Tree(fast_tree_address) # ,format=0
print("length of tree is ",len(fast_tree))


# In[150]:


fast_tree.write()


# In[151]:


list_species=[]
for node in fast_tree.traverse(strategy="postorder"):
    if node.is_leaf() :    
        list_species.append(node.name)
        
print("number of leaves, species, is:",len(list_species))


# In[152]:



list_species_tax=[]
list_species_omaID=[]
no_list=[]
for node in fast_tree.traverse(strategy="postorder"):
    if node.is_leaf() :    
        node_name = node.name
        list_species_omaID.append(node_name)
        if node_name in omaID_taxonID:
            
            node_name_tax = omaID_taxonID[node_name]
            node.name = node_name_tax            
            list_species_tax.append(node_name_tax)
        else:
            no_list.append(node_name)
            
print("number of leaves, species, is:",len(list_species_omaID),len(list_species_tax),len(no_list))


# In[153]:


fast_tree.write()


# In[154]:


print(len(fast_tree))
fast_tree.prune(list_species_tax,preserve_branch_length=True)
print(len(fast_tree))


# In[155]:


fast_tree.write()


# In[ ]:





# In[156]:


# list_sepcies_mammal = [10020, 10029, 10036, 10090, 10116, 10141, 10160, 10181, 103600, 116960, 118797, 127582, 132908, 13616, 1706337, 185453, 1868482, 191816, 230844, 2715852, 27622, 29073, 29088, 29139, 30522, 30538, 30611, 310752, 32536, 336983, 34884, 37032, 37293, 379532, 38626, 391180, 39432, 40151, 419612, 43179, 43346, 46360, 51298, 56216, 59463, 59472, 59479, 60711, 61621, 61622, 61853, 72004, 77932, 885580, 89673, 9258, 9305, 9365, 9402, 9407, 9483, 9515, 9531, 9541, 9544, 9545, 9555, 9568, 9595, 9597, 9598, 9601, 9606, 9615, 9627, 9643, 9646, 9669, 9678, 9685, 9696, 9708, 9713, 9739, 9749, 9755, 9764, 9770, 9785, 9796, 9823, 9838, 9880, 9886, 9913, 9915, 9925, 9940, 9986, 9995]
# #outgroups=[8364, 9031] # frog chicken
# list_sepcies_str=[str(i) for i in list_sepcies_mammal]

# print(len(fast_tree))
# fast_tree.prune(list_sepcies_str)
# print(len(fast_tree))

# fast_tree.set_outgroup("9258") 
# # frog 8364
# # 9258  Ornithorhynchus anatinus (Duckbill platypus)


# # Reading  NCBI tree

# In[157]:



ncbi_taxon_list = list_species_tax


ncbi = ete3.NCBITaxa()  # first time download in ~/.etetoolkit/
ncbi_sub_tree = ncbi.get_topology(ncbi_taxon_list)

print(len(ncbi_sub_tree))


# In[ ]:





# In[158]:


ncbi_sub_tree.write()


# # RF distance

# In[159]:


# pron fastoma


# print(len(fast_tree))
# fast_tree.prune(list_sepcies_str)
# print(len(fast_tree))



# In[160]:


def parition_editor(partitions_input,intersect_leaves):
    pivot_species= sorted(intersect_leaves)[0]
    partitions_ed=[]
    for parition in  partitions_input:
        #parition=set(parition)
        if len(parition)<len(intersect_leaves)/2 and len(parition)>1:  # removing partitions with length of 0/1 
            partitions_ed.append(parition)
        if len(parition)>len(intersect_leaves)/2 and len(parition)<len(intersect_leaves)-1:
            partitions_ed.append(tuple(intersect_leaves-set(parition)))
        if len(parition)==len(intersect_leaves)/2 and pivot_species not in parition: # make sure alwyas one half is there
            partitions_ed.append(tuple(intersect_leaves-set(parition)))
            
    # partitions = [e for e in edges if n > len(e) > 1]
    # partitions = list(map(lambda p: p if len(p) <= n/2 else all_species -set(p), partitions)      
    # if len(p) < n/2 or len(p) == n/2 and pivot_species in p: return p
    # else: return all_species - set(p)   
    return partitions_ed
        
        


# In[162]:


fast_tree.set_outgroup("9913")  # bovin


# In[163]:


fast_tree.write()


# In[161]:


fast_tree.set_outgroup("9913")  # bovin



intersection_sets= set(ncbi_taxon_list)

out_fast_ncbi = fast_tree.robinson_foulds(ncbi_sub_tree)
(rf_fast_ncbi, maxparts_fast_ncbi, common_attrs, partitions_fast_raw, partitions_ncbi_raw, discard_t1, discard_t2)=out_fast_ncbi 
#partitions_fast = set(parition_editor(partitions_fast_raw, intersection_sets))
#partitions_ncbi = set(parition_editor(partitions_ncbi_raw, intersection_sets))

print("RF distance is %s over a total of %s" %(rf_fast_ncbi, maxparts_fast_ncbi))
print("in fast:",len(partitions_fast),", in ncbi:",len(partitions_ncbi)) # , "Common Partitions",len(common_attrs_paper),)
print("in both fast and ncbi:", len(partitions_ncbi & partitions_fast))
print("only in ncbi, not in fast:", len(partitions_ncbi-partitions_fast)) # order not sure  http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#robinson-foulds-distance
print("only in fast, not in ncbi:", len(partitions_fast-partitions_ncbi))



# In[146]:


no_list=[]
for node in fast_tree.traverse(strategy="postorder"):
    if node.is_leaf() :           
        #fast_tree_leaves.append(node.name)
        #old_node_name
        if node.name in common_taxonID:
            node_name = common_taxonID[node.name]
            node.name=node_name + "_"+node.name
        elif node.name=='30608' : node.name="Lesser mouse lemur_30608"
        else:
            no_list.append(node.name)
        
print(no_list)


# In[ ]:





# In[125]:


no_list=[]
for node in ncbi_sub_tree.traverse(strategy="postorder"):
    if node.is_leaf() :           
        if node.name in common_taxonID:
            node_name = common_taxonID[node.name]
            node.name=node_name + "_"+node.name
        elif node.name=='30608' : node.name="Lesser mouse lemur_30608"
        else:
            no_list.append(node.name)
        
print(no_list)


# In[126]:


ncbi_sub_tree.write()


# In[127]:


fast_tree.write()


# ## fast tree with oma species

# In[53]:


# ##oma_standalon backbone species names
# import pyoma.browser.db as db
# project_folder_root="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/mammal/old/mml6/"
# oma_database_address=project_folder_root+"omamer_database/oma_path/OmaServer.h5"
# oma_db = db.Database(oma_database_address);
# print(oma_db.tax.as_dict())
#dic1['children'][0]['children']

dic_oma={'S0005':'rabbit', 'S0003':'human', 'S0006':'redfox', 'S0004':'opossum', 'S0001':'chicken','S0002':'frog'}
#{'rabbit':'S0005', 'human':'S0003', 'redfox': 'S0006', 'opossum':'S0004', 'chicken':'S0001','frog':'S0002'}


# In[123]:


project_folder_root="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/mammal/old/mml6/"
fast_tree_address=project_folder_root+"fastoma_core/_80_msa_concatanated.txt.contree_collapsed_94_0.001"
fast_tree= ete3.Tree(fast_tree_address) # ,format=0
print("length of tree is ",len(fast_tree))


# In[124]:


no_list=[]
for node in fast_tree.traverse(strategy="postorder"):
    if node.is_leaf() :           
        #fast_tree_leaves.append(node.name)
        #old_node_name
        if node.name in dic_oma:
            node_name = dic_oma[node.name]
            node.name=node_name + "_"+node.name
        elif node.name in common_taxonID:
            node_name = common_taxonID[node.name]
            node.name=node_name + "_"+node.name
        else:
            no_list.append(node.name)
        
print(no_list)


# In[125]:


fast_tree.write()


# In[ ]:





# In[ ]:





# ## Standard OMA  -slow

# In[107]:



slow_tree_address="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/mammal/old/mml_v1/marker_genes/tree/"

#project_folder = project_folder_root +"proteome/" # fastoma_out
slow_tree_address=slow_tree_address+"/_100_msa_concatanated.txt.contree"  # _collapsed_brnch_0.01# 1e-05" # 
slow_tree_address += "_collapsed2_95"

slow_tree= ete3.Tree(slow_tree_address) # ,format=0
print(len(slow_tree))


# In[108]:


no_list=[]
slow_species_taxon=[]
for node in slow_tree.traverse(strategy="postorder"):
    if node.is_leaf() :           
        if node.name in omaID_taxonID:
            node_name = omaID_taxonID[node.name]
            node.name=  node_name 
            slow_species_taxon.append(node_name)
        else:
            no_list.append(node.name)
        
print(no_list)


# In[ ]:





# In[109]:


fast_tree= ete3.Tree(fast_tree_address) # ,format=0
print("length of tree is ",len(fast_tree))


intersection_slow_fast=set(slow_species_taxon) & set(list_sepcies_str)

print(len(fast_tree),len(slow_tree))
fast_tree.prune(intersection_slow_fast)
slow_tree.prune(intersection_slow_fast)
print(len(fast_tree),len(slow_tree))


# In[ ]:





# In[111]:


fast_tree.write()


# In[118]:


slow_tree.write()


# In[ ]:





# In[120]:



ncbi_taxon_list = list_sepcies_str


ncbi = ete3.NCBITaxa()  # first time download in ~/.etetoolkit/
ncbi_sub_tree = ncbi.get_topology(ncbi_taxon_list)

print(len(ncbi_sub_tree))


# In[121]:


slow_tree.set_outgroup("9258") 


out_slow_ncbi = slow_tree.robinson_foulds(ncbi_sub_tree)
(rf_slow_ncbi, maxparts_slow_ncbi, common_attrs, partitions_slow_raw, partitions_ncbi_raw, discard_t1, discard_t2)=out_slow_ncbi 
partitions_slow = set(parition_editor(partitions_slow_raw, intersection_sets))
partitions_ncbi = set(parition_editor(partitions_ncbi_raw, intersection_sets))

print("RF distance is %s over a total of %s" %(rf_slow_ncbi, maxparts_slow_ncbi))
print("in slow:",len(partitions_slow),", in ncbi:",len(partitions_ncbi)) # , "Common Partitions",len(common_attrs_paper),)
print("in both slow and ncbi:", len(partitions_ncbi & partitions_slow))
print("only in ncbi, not in slow:", len(partitions_ncbi-partitions_slow)) # order not sure  http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#robinson-foulds-distance
print("only in slow, not in ncbi:", len(partitions_slow-partitions_ncbi))


# In[122]:



intersection_sets= set(list_sepcies_str)

out_fast_ncbi = fast_tree.robinson_foulds(ncbi_sub_tree)
(rf_fast_ncbi, maxparts_fast_ncbi, common_attrs, partitions_fast_raw, partitions_ncbi_raw, discard_t1, discard_t2)=out_fast_ncbi 
partitions_fast = set(parition_editor(partitions_fast_raw, intersection_sets))
partitions_ncbi = set(parition_editor(partitions_ncbi_raw, intersection_sets))

print("RF distance is %s over a total of %s" %(rf_fast_ncbi, maxparts_fast_ncbi))
print("in fast:",len(partitions_fast),", in ncbi:",len(partitions_ncbi)) # , "Common Partitions",len(common_attrs_paper),)
print("in both fast and ncbi:", len(partitions_ncbi & partitions_fast))
print("only in ncbi, not in fast:", len(partitions_ncbi-partitions_fast)) # order not sure  http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#robinson-foulds-distance
print("only in fast, not in ncbi:", len(partitions_fast-partitions_ncbi))



# In[ ]:





# In[112]:


print("only in fastOMA:", len((partitions_fast - partitions_ncbi) - partitions_slow))
print("only in fastOMA & slow:",len((partitions_fast & partitions_slow) - partitions_ncbi))
print("only in fastOMA & NCBI:",len((partitions_fast & partitions_ncbi) - partitions_slow))
print("only in slow & NCBI:",len((partitions_slow & partitions_ncbi) -  partitions_fast))
print("only in fastOMA & NCBI & slow:",len((partitions_fast & partitions_ncbi) & partitions_slow))
print("only in slow:",len((partitions_slow -partitions_fast) - partitions_ncbi))
print("only in NCBI:",len((partitions_ncbi - partitions_fast) - partitions_slow))


# In[ ]:





# ### three way comparison

# In[114]:


# fast_tree= ete3.Tree(fast_treee_address_out, format=1)
# fast_tree.prune(intersect_leaves)

# fast_tree.set_outgroup(outgroup) # CALPUG MOUS  BUCCAP CHLMAC

out_paper_ncbi = bird_paper_tree.robinson_foulds(ncbi_sub_tree) # , expand_polytomies = True  #,polytomy_size_limit=20  ,unrooted_trees=True # , expand_polytomies = True
(rf_paper_ncbi, maxparts_paper_ncbi, common_attrs, partitions_paper_raw, partitions_ncbi_raw, discard_t1, discard_t2)=out_paper_ncbi 
partitions_paper = set(parition_editor(partitions_paper_raw, intersect_leaves))
partitions_ncbi = set(parition_editor(partitions_ncbi_raw, intersect_leaves))

out_fast_ncbi = fast_tree.robinson_foulds(ncbi_sub_tree)  # ,unrooted_trees=True # , expand_polytomies = True
(rf_fast_ncbi, maxparts_fast_ncbi, common_attrs, partitions_fast_raw, partitions_ncbi_raw2, discard_t1, discard_t2)=out_fast_ncbi 
partitions_fast = set(parition_editor(partitions_fast_raw, intersect_leaves))
partitions_ncbi2 = set(parition_editor(partitions_ncbi_raw2, intersect_leaves))
assert len(partitions_ncbi2-partitions_ncbi)==0
out_fast_paper = fast_tree.robinson_foulds(bird_paper_tree)  # ,unrooted_trees=True # , expand_polytomies = True
(rf_fast_paper, maxparts_fast_paper, common_attrs, partitions_fast_raw2, partitions_paper_raw2, discard_t1, discard_t2)=out_fast_paper
partitions_paper2 = set(parition_editor(partitions_paper_raw2, intersect_leaves))
partitions_fast2 = set(parition_editor(partitions_fast_raw2, intersect_leaves))
assert len(partitions_fast2-partitions_fast)==0
assert len(partitions_paper2-partitions_paper)==0


values_comparison=[len((partitions_fast - partitions_ncbi) - partitions_paper),
                   len((partitions_fast & partitions_paper) - partitions_ncbi2),
                   len((partitions_fast & partitions_ncbi) - partitions_paper),
                   len((partitions_paper & partitions_ncbi) -  partitions_fast),
                   len((partitions_fast & partitions_ncbi) &  partitions_paper),
                   len((partitions_paper - partitions_fast) - partitions_ncbi),
                   len((partitions_ncbi - partitions_fast) - partitions_paper) ]
print(values_comparison)


# In[ ]:





# In[ ]:


sorted(tresh_list)


# In[ ]:





# In[ ]:





# #  collapase   branch length & support value

# In[46]:


# dendropy

import dendropy

fast_tree_raw_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/mammal/"

fast_tree_raw_address += "mml27/fastoma_core/_100_msa_concatanated.txt.contree"

tree = dendropy.Tree.get_from_path(fast_tree_raw_address,"newick") # newick nexus

print(len(tree))

support_treshold = 94
for idx, node in enumerate(tree):
    if not node.is_leaf() and node.label:
        support_value = int(node.label)
        if support_value <= support_treshold:
            edge= node.edge   #print(edge.head_node.child_nodes())
            edge.collapse()
print(len(tree))

threshold=1e-3 

for e in tree.postorder_edge_iter():    
    if e.length is None or (e.length <= threshold) and e.is_internal():
        e.collapse()
print(len(tree))

fast_treee_address_out=fast_tree_raw_address+"_collapsed_"+str(support_treshold)+"_"+str(threshold)
tree.write_to_path(fast_treee_address_out, "newick")
print(fast_treee_address_out)


# In[ ]:





# # ****

# # others

# ### plot different collapsing  support value
# 

# In[ ]:



tresh=60 
project_folder_root="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3b/"
tresh_list=[60,80,90,95,100]
tresh_list_name=["_collapsed_"+str(i) for i in tresh_list] 
tresh_list_name=[""]+tresh_list_name

values_comparison_all=[]
for tresh in tresh_list_name:
    fast_tree_address=project_folder_root+"iqtree5/_100_msa_concatanated.txt_copy1.contree"+tresh

    fast_tree= ete3.Tree(fast_tree_address)
    #print("len fast tree",len(fast_tree))


    fast_tree.prune(intersect_leaves)
    #print("len fast tree after prune",len(fast_tree))

    fast_tree.set_outgroup(outgroup) # CALPUG MOUS  BUCCAP  CHLMAC

    out_paper_ncbi = bird_paper_tree.robinson_foulds(ncbi_sub_tree) # , expand_polytomies = True  #,polytomy_size_limit=20  ,unrooted_trees=True # , expand_polytomies = True
    (rf_paper_ncbi, maxparts_paper_ncbi, common_attrs, partitions_paper_raw, partitions_ncbi_raw, discard_t1, discard_t2)=out_paper_ncbi 
    partitions_paper = set(parition_editor(partitions_paper_raw, intersect_leaves))
    partitions_ncbi = set(parition_editor(partitions_ncbi_raw, intersect_leaves))

    out_fast_ncbi = fast_tree.robinson_foulds(ncbi_sub_tree)  # ,unrooted_trees=True # , expand_polytomies = True
    (rf_fast_ncbi, maxparts_fast_ncbi, common_attrs, partitions_fast_raw, partitions_ncbi_raw2, discard_t1, discard_t2)=out_fast_ncbi 
    partitions_fast = set(parition_editor(partitions_fast_raw, intersect_leaves))
    partitions_ncbi2 = set(parition_editor(partitions_ncbi_raw2, intersect_leaves))
    assert len(partitions_ncbi2-partitions_ncbi)==0
    out_fast_paper = fast_tree.robinson_foulds(bird_paper_tree)  # ,unrooted_trees=True # , expand_polytomies = True
    (rf_fast_paper, maxparts_fast_paper, common_attrs, partitions_fast_raw2, partitions_paper_raw2, discard_t1, discard_t2)=out_fast_paper
    partitions_paper2 = set(parition_editor(partitions_paper_raw2, intersect_leaves))
    partitions_fast2 = set(parition_editor(partitions_fast_raw2, intersect_leaves))
    assert len(partitions_fast2-partitions_fast)==0
    assert len(partitions_paper2-partitions_paper)==0


    values_comparison=[len((partitions_fast - partitions_ncbi) - partitions_paper),
                       len((partitions_fast & partitions_paper) - partitions_ncbi2),
                       len((partitions_fast & partitions_ncbi) - partitions_paper),
                       len((partitions_paper & partitions_ncbi) -  partitions_fast),
                       len((partitions_fast & partitions_ncbi) &  partitions_paper),
                       len((partitions_paper - partitions_fast) - partitions_ncbi),
                       len((partitions_ncbi - partitions_fast) - partitions_paper) ]
    #sum1=values_comparison[0]+values_comparison[4]+values_comparison[5]+values_comparison[6]
    print(values_comparison) # ,sum1
    values_comparison_all.append(values_comparison)
values_comparison_all=np.array(values_comparison_all)
#print(values_comparison_all)

#t_all=np.transpose(values_comparison_all)
lenged1=['fastOMA','fastOMA & paper',  'fastOMA & NCBI',
  'paper & NCBI', 'fastOMA & NCBI & paper',  'paper', 'NCBI']

import pandas as pd
df = pd.DataFrame(data=values_comparison_all)
df.columns=lenged1   #ind = [0,60,80,90,95,100] # ,"BirdPaper"
df.index=["no"]+tresh_list
df


# ### plot different collapsing branch length

# In[ ]:



#tresh=60 
project_folder_root="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3b/"
tresh_list=[1e-05,0.0001,0.001,0.005,0.01,0.1]
tresh_list_name=["_collapsed_brnch_"+str(i) for i in tresh_list] 
tresh_list_name=[""]+tresh_list_name

values_comparison_all=[]
for tresh in tresh_list_name:
    fast_tree_address=project_folder_root+"iqtree5/_100_msa_concatanated.txt_copy1.contree"+tresh

    fast_tree= ete3.Tree(fast_tree_address)
    #print("len fast tree",len(fast_tree))


    fast_tree.prune(intersect_leaves)
    #print("len fast tree after prune",len(fast_tree))

    fast_tree.set_outgroup(outgroup) # CALPUG MOUS  BUCCAP

    out_paper_ncbi = bird_paper_tree.robinson_foulds(ncbi_sub_tree) # , expand_polytomies = True  #,polytomy_size_limit=20  ,unrooted_trees=True # , expand_polytomies = True
    (rf_paper_ncbi, maxparts_paper_ncbi, common_attrs, partitions_paper_raw, partitions_ncbi_raw, discard_t1, discard_t2)=out_paper_ncbi 
    partitions_paper = set(parition_editor(partitions_paper_raw, intersect_leaves))
    partitions_ncbi = set(parition_editor(partitions_ncbi_raw, intersect_leaves))

    out_fast_ncbi = fast_tree.robinson_foulds(ncbi_sub_tree)  # ,unrooted_trees=True # , expand_polytomies = True
    (rf_fast_ncbi, maxparts_fast_ncbi, common_attrs, partitions_fast_raw, partitions_ncbi_raw2, discard_t1, discard_t2)=out_fast_ncbi 
    partitions_fast = set(parition_editor(partitions_fast_raw, intersect_leaves))
    partitions_ncbi2 = set(parition_editor(partitions_ncbi_raw2, intersect_leaves))
    assert len(partitions_ncbi2-partitions_ncbi)==0
    out_fast_paper = fast_tree.robinson_foulds(bird_paper_tree)  # ,unrooted_trees=True # , expand_polytomies = True
    (rf_fast_paper, maxparts_fast_paper, common_attrs, partitions_fast_raw2, partitions_paper_raw2, discard_t1, discard_t2)=out_fast_paper
    partitions_paper2 = set(parition_editor(partitions_paper_raw2, intersect_leaves))
    partitions_fast2 = set(parition_editor(partitions_fast_raw2, intersect_leaves))
    assert len(partitions_fast2-partitions_fast)==0
    assert len(partitions_paper2-partitions_paper)==0


    values_comparison=[len((partitions_fast - partitions_ncbi) - partitions_paper),
                       len((partitions_fast & partitions_paper) - partitions_ncbi2),
                       len((partitions_fast & partitions_ncbi) - partitions_paper),
                       len((partitions_paper & partitions_ncbi) -  partitions_fast),
                       len((partitions_fast & partitions_ncbi) &  partitions_paper),
                       len((partitions_paper - partitions_fast) - partitions_ncbi),
                       len((partitions_ncbi - partitions_fast) - partitions_paper) ]
    print(values_comparison)
    values_comparison_all.append(values_comparison)
values_comparison_all=np.array(values_comparison_all)
#print(values_comparison_all)

#t_all=np.transpose(values_comparison_all)
lenged1=['fastOMA','fastOMA & paper',  'fastOMA & NCBI',
  'paper & NCBI', 'fastOMA & NCBI & paper',  'paper', 'NCBI']

import pandas as pd
df = pd.DataFrame(data=values_comparison_all)
df.columns=lenged1   #ind = [0,60,80,90,95,100] # ,"BirdPaper"
df.index=["no"]+tresh_list
df


# ## read ncbi tree

# In[138]:


list_oma=["OTOGA", "GORGO", "MACNE", "BOVIN", "MOUSE"]
ncbi_taxon_list= [omaID_taxonID[i] for i in list_oma]




ncbi = ete3.NCBITaxa()  # first time download in ~/.etetoolkit/
ncbi_sub_tree = ncbi.get_topology(ncbi_taxon_list)

print(len(ncbi_sub_tree))


# In[139]:


no_list=[]
for node in ncbi_sub_tree.traverse(strategy="postorder"):
    if node.is_leaf() :           
        if node.name in taxonID_omaID:
            node_name = taxonID_omaID[node.name]
            node.name=node_name
            


# In[140]:


ncbi_sub_tree.write()


# In[ ]:


#10 # '((((((((MANLE:1,MACNE:1)1:1,COLAP:1)1:1,(PANTR:1,GORGO:1)1:1)1:1,SAIBB:1)1:1,CARSF:1)1:1,OTOGA:1)1:1,MOUSE:1)1:1,BOVIN:1);'
#5  # '((((MACNE:1,GORGO:1)1:1,OTOGA:1)1:1,MOUSE:1)1:1,BOVIN:1)myroot:0;'

