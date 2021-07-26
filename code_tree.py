#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/python3

import numpy as np
from sys import argv

from datetime import datetime
#import concurrent.futures

# import ast
#  for development 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import ete3


# In[2]:


def read_taxonID_map(omaID_address,bird6ID_address):
    
    omaID_file = open(omaID_address,'r')
    taxonID_omaID={}
    omaID_taxonID={}
    for line in omaID_file:
        line_strip = line.strip()
        if line_strip.startswith('#'):
            pass
            #header_lines_list.append(line_strip)
        else:
            line_parts=line_strip.split('\t')

            omaID = line_parts[0]
            taxonID = int(line_parts[1])
            taxonID_omaID[taxonID]=omaID
            omaID_taxonID[omaID]=taxonID
            
    omaID_file.close()
        
    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- The map for taxonID omaID of",len(taxonID_omaID),"records have read.") 
    
    
    bird6ID_file = open(bird6ID_address,'r')
    taxonID_bird6ID={}
    bird6ID_taxonID= {}
    for line in bird6ID_file:
        line_strip = line.strip()
        if line_strip.startswith('Or'):
            pass
            line_strip1=line_strip
            #header_lines_list.append(line_strip)
        else:
            line_parts=line_strip.split('\t')

            bird6ID = line_parts[6]
            taxonID = int(line_parts[8])
            taxonID_bird6ID[taxonID]=bird6ID
            bird6ID_taxonID[bird6ID]=taxonID
            
    bird6ID_file.close()
        
    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- The map for taxonID bird6ID of",len(taxonID_bird6ID),"records have read.") 
  
    return (taxonID_omaID,taxonID_bird6ID,omaID_taxonID,bird6ID_taxonID)
    
            


# In[3]:



datasets_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/"
# oma_database_address = datasets_address + "OmaServer.h5"
# hog_og_map_address = datasets_address + "hog_og_map.dic"
omaID_address = datasets_address+"oma-species.txt"
bird6ID_address = datasets_address+"info.tsv"
project_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3a/ST/f4_100S/" 

bird_tree_address="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3a/hogmapX/out4/_tree_filt_7.txt"


# In[4]:


(taxonID_omaID,taxonID_bird6ID,omaID_taxonID,bird6ID_taxonID) = read_taxonID_map(omaID_address,bird6ID_address)
# left is the ky
# taxonID_omaID  {taxonID:omaID, .. }
# taxonID_omaID[1111]=ABABA

taxonID_map = {**taxonID_omaID,**taxonID_bird6ID }

print(len(taxonID_omaID),len(taxonID_bird6ID),len(taxonID_map))
#taxonID_list= list(taxonID_map.keys())


# # Reading FastOMA bird  tree

# In[5]:


bird_tree_raw= ete3.Tree(bird_tree_address)
len(bird_tree_raw)
#bird_tree.write()


# In[6]:


bird_tree_leaves_omaID_bird6ID=[]
for node in bird_tree_raw.traverse(strategy="postorder"):
    if node.is_leaf() : # why ?
        bird_tree_leaves_omaID_bird6ID.append(node.name)
    
print("bird_hog_tree_leaf_count",len(bird_tree_leaves_omaID_bird6ID))
bird_tree_leaves_omaID_bird6ID[:3]


# In[7]:


bird_tree_leaves_taxonID = []
for i3 in bird_tree_leaves_omaID_bird6ID:
    if i3 in omaID_taxonID:
        taxonID=omaID_taxonID[i3]
    if i3  in bird6ID_taxonID:
        taxonID= bird6ID_taxonID[i3]
    # if it was in the both I save the bird ID    
    bird_tree_leaves_taxonID.append(taxonID)

bird_tree_leaves_taxonID_unq=list(set(bird_tree_leaves_taxonID))

print(len(bird_tree_leaves_taxonID),len(bird_tree_leaves_taxonID_unq))
bird_tree_leaves_taxonID_unq[:3]


# In[8]:


both_bird_oma_taxonID=[]
for i in bird_tree_leaves_taxonID_unq:
    if bird_tree_leaves_taxonID.count(i)>1:
        both_bird_oma_taxonID.append(i)
print(both_bird_oma_taxonID)

both_bird_oma_omaID=[]
for i in both_bird_oma_taxonID:
    both_bird_oma_omaID.append(taxonID_omaID[i])

print(both_bird_oma_omaID)

bird_tree_leaves_omaID_bird6ID_uniq=[i for i in bird_tree_leaves_omaID_bird6ID if i not in both_bird_oma_omaID]
print(len(bird_tree_leaves_omaID_bird6ID),len(bird_tree_leaves_omaID_bird6ID_uniq))


bird_tree= ete3.Tree(bird_tree_address)
bird_tree.prune(bird_tree_leaves_omaID_bird6ID_uniq)

print(len(bird_tree_raw),len(bird_tree))


# In[9]:


#bird_tree.write()


# # Reading  NCBI tree

# In[10]:


ncbi = ete3.NCBITaxa()  # first time download in ~/.etetoolkit/
ncbi_sub_tree = ncbi.get_topology(bird_tree_leaves_taxonID_unq)


# In[11]:


## change the node name from NCBI taxon id (integer) to   omaID (5-letter) or bird6ID (6-letter)

for node in ncbi_sub_tree.traverse(strategy="postorder"):
    node.name2 = node.name
    if node.is_leaf() : # why ? and int(node.name) in taxonID_map
        node.name = taxonID_map[int(node.name)]

#current_time = datetime.now().strftime("%H:%M:%S")
#print(current_time, 
print("The NCBI taxanomy is read and the leaves names changed to OMA/bird6 ID containing")
print(len(ncbi_sub_tree)) 
#print(ncbi_sub_tree.get_ascii(attributes=["name"]))     

#ncbi_sub_tree.write() 


# In[ ]:





# In[ ]:





# In[ ]:





# In[36]:


ncbi_sub_tree.write() 


# In[ ]:





# In[13]:


#bird_tree.write()


# # compareing fastOMA bird and NCBI
# 

# In[41]:


bird_tree2=bird_tree
bird_tree2.set_outgroup("XENLA")

#bird_tree2.rerootatedge(edge_XENLA)

out = bird_tree2.robinson_foulds(ncbi_sub_tree) # , expand_polytomies = True  #,polytomy_size_limit=20  ,unrooted_trees=True # , expand_polytomies = True

##  grep "min_comparison" /work/FAC/FBM/DBC/cdessim2/default/smajidi1/software/miniconda3/lib/python3.8/site-packages/ete3/coretype/tree.py

(rf, max_parts, common_attrs, edges1, edges2, discard_t1, discard_t2)=out 
# print(out[0],out[1])
# print(len(list(out[2])),len(list(out[3])),len(list(out[4])),out[5],out[6])
print("RF distance is %s over a total of %s" %(rf, max_parts))

print("Len of common_attrs",len(common_attrs),", len edges1",len(edges1),", len edges2",len(edges2))

e2_1=edges2 - edges1
print("Partitions in ncbi no in fastOMA:", len(e2_1)) # order not sure  http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#robinson-foulds-distance
e1_2=edges1 - edges2
print("Partitions in fastOMA no in ncbi:", len(e1_2))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# # Reading bird paper tree

# In[26]:


bird_paper_tree_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/bird_paper_tree.txt"
bird_paper_tree_raw= ete3.Tree(bird_paper_tree_address, format=1)
#len(bird_paper_tree_raw)


bird_paper_tree_leaves=[]
for node in bird_paper_tree_raw.traverse(strategy="postorder"):
    if node.is_leaf() : # why ?
        bird_paper_tree_leaves.append(node.name)
print(len(bird_paper_tree_leaves))
bird_paper_tree_leaves[:3]



# In[27]:



bird_SCINAME_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/"+"SCINAME_all_pure_under.txt"

def read_bird_SCINAME(bird_SCINAME_address):
    
    bird_SCINAME_file = open(bird_SCINAME_address,'r')
    bird_SCINAME_list = []
    
    for line in bird_SCINAME_file:
        line_strip = line.strip()
        bird_SCINAME_list.append(line_strip)
        
    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- The bird science name  for ",len(bird_SCINAME_list),"records have read.") 
    return bird_SCINAME_list


bird_SCINAME_list = read_bird_SCINAME(bird_SCINAME_address)


# In[28]:


#bird_paper_tree_raw.prune(bird_SCINAME_list)

notin_paper_tree= []
for i in bird_SCINAME_list:
    if i not in bird_paper_tree_leaves: 
        notin_paper_tree.append(i)
print(notin_paper_tree)
#notin_paper_tree=["Cercotrichas_coryphaeus","Corvus_cornix","Eolophus_roseicapillus","Nannopterum_auritus","Nannopterum_brasilianus","Nannopterum_harrisi","Urile_pelagicus"]


# In[29]:



bird_SCINAME_list_filt=[i for i in bird_SCINAME_list if i not in notin_paper_tree]
print(len(bird_SCINAME_list),len(bird_SCINAME_list_filt) )
    


# In[30]:


bird_paper_tree= ete3.Tree(bird_paper_tree_address, format=1)

bird_paper_tree.prune(bird_SCINAME_list_filt)
len(bird_paper_tree)


# In[33]:


bird_paper_tree_nodes =[]


for node in bird_paper_tree.traverse(strategy="postorder"):
    node.name2 = node.name
    if node.is_leaf() : # why ? and int(node.name) in taxonID_map
        
        node_name = node.name
        node_name_split= node_name.split('_')
        six_letter_name= ''.join([i[:3].upper() for i in node_name_split])
        if len(six_letter_name)<5:
            print(node_name)
        else:
            node.name= six_letter_name
            #print(six_letter_name)
            bird_paper_tree_nodes.append(six_letter_name)

len(bird_paper_tree_nodes)


# In[35]:


bird_paper_tree.write()


# In[ ]:





# In[ ]:





# In[ ]:


# compareing bird paper and NCBI


# In[ ]:





# In[40]:



out = bird_paper_tree.robinson_foulds(ncbi_sub_tree) # , expand_polytomies = True  #,polytomy_size_limit=20  ,unrooted_trees=True # , expand_polytomies = True


(rf, max_parts, common_attrs, edges1, edges2, discard_t1, discard_t2)=out 
# print(out[0],out[1])
# print(len(list(out[2])),len(list(out[3])),len(list(out[4])),out[5],out[6])
print("RF distance is %s over a total of %s" %(rf, max_parts))

#print("Len of common_attrs",len(common_attrs),", len edges1",len(edges1),", len edges2",len(edges2))

e2_1=edges2 - edges1
print("Partitions in ncbi not in bird_paper:", len(e2_1)) # order not sure  http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#robinson-foulds-distance
e1_2=edges1 - edges2
print("Partitions in bird_paper not in ncbi:", len(e1_2))


# In[ ]:





# In[43]:




out = bird_paper_tree.robinson_foulds(bird_tree2) # , expand_polytomies = True  #,polytomy_size_limit=20  ,unrooted_trees=True # , expand_polytomies = True


(rf, max_parts, common_attrs, edges1, edges2, discard_t1, discard_t2)=out 
# print(out[0],out[1])
# print(len(list(out[2])),len(list(out[3])),len(list(out[4])),out[5],out[6])
print("RF distance is %s over a total of %s" %(rf, max_parts))

#print("Len of common_attrs",len(common_attrs),", len edges1",len(edges1),", len edges2",len(edges2))

e2_1=edges2 - edges1
print("Partitions in fastOMA not in bird_paper:", len(e2_1)) # order not sure  http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#robinson-foulds-distance
e1_2=edges1 - edges2
print("Partitions in bird_paper not in fastOMA:", len(e1_2))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


list_i_XENLA_0=[]
list_i_XENLA_1=[]
for i in edges1:
    list_i=list(i)
    if len(list_i)>2:
        print(len(list_i))
    for j0 in list_i[0]:
        if j0=='XENLA':
            list_i_XENLA_0.append(list_i)
    for j1 in list_i[1]:
        if j1=='XENLA':
            list_i_XENLA_1.append(list_i)


# In[ ]:


edge_XENLA=list(list_i_XENLA_1[0])[1]


# In[ ]:


#list_i_XENLA_1


# In[ ]:


#(rf, rf_max, common_attrs, names, edges_t1, edges_t2, discarded_edges_t1) =out

print("RF distance is %s over a total of %s" %(rf, max_parts))
#print("Partitions in tree2 that were not found in tree1:", edges1 - edges2)

e2_1=edges2 - edges1
print("Partitions in tree1 that were not found in tree2:", len(e2_1))


# In[ ]:


print(len(edges2),len(edges1))
e2_1=edges2 - edges1

# for i in e2_1:
#     #print(len(i))
#     for k in i:
#         len(k)
e2_1_list=list(e2_1)
len(e2_1_list[0]), len(list(e2_1_list[0])[0]), len(list(e2_1_list[0])[1]),list(e2_1_list[0])[1]


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


result= bird_tree2.compare(ncbi_sub_tree)


# In[ ]:


result[“rf”] = robinson-foulds distance between the two trees. (average of robinson-foulds distances if target tree contained duplication and was split in several subtrees)
result[“max_rf”] = Maximum robinson-foulds distance expected for this comparison
result[“norm_rf”] = normalized robinson-foulds distance (from 0 to 1)
result[“effective_tree_size”] = the size of the compared trees, which are pruned to the common shared nodes.
result[“ref_edges_in_source”] = compatibility score of the target tree with respect to the source tree (how many edges in reference are found in the source)
result[“source_edges_in_ref”] = compatibility score of the source tree with respect to the reference tree (how many edges in source are found in the reference)
result[“source_subtrees”] = number of subtrees in the source tree (1 if do not contain duplications)
result[“common_edges”] = a set of common edges between source tree and reference
result[“source_edges”] = the set of edges found in the source tree
result[“ref_edges”] = the set of edges found in the reference tree
result[“treeko_dist”] = TreeKO speciation distance for comparisons including duplication nodes.


# In[ ]:


result["norm_rf"]


# In[ ]:


t1 = ete3.Tree('(((a,(b,i)),c),((e, f), g));')
#t1 = ete3.Tree('(((a,c),b), ((e, f), g));')
t2 = ete3.Tree('(((a,c),b), ((e, f), g));')


# In[ ]:


t1.write(), t2.write()


# In[ ]:


#t1.compare(t2)
out= t1.robinson_foulds(t2)

(rf, max_parts, common_attrs, edges1, edges2, discard_t1, discard_t2)=out 
# print(out[0],out[1])
# print(len(list(out[2])),len(list(out[3])),len(list(out[4])),out[5],out[6])
print("RF distance is %s over a total of %s" %(rf, max_parts))

print("Len of common_attrs",len(common_attrs),", len edges1",len(edges1),", len edges2",len(edges2))

e2_1=edges2 - edges1
print("Partitions in tr2 no in tr1:", len(e2_1)) # order not sure  http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#robinson-foulds-distance
e1_2=edges1 - edges2
print("Partitions in tr1 no in tr2:", len(e1_2))


# In[ ]:


e2_1,e1_2


# In[ ]:


edges1


# In[ ]:


edges2


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


bird_tree_leaves_omaID_bird6ID_uniq_set=set(bird_tree_leaves_omaID_bird6ID_uniq)
bird_paper_tree_nodes_set=set(bird_paper_tree_nodes)

mine_bird_paper= list(bird_paper_tree_nodes_set.intersection(bird_tree_leaves_omaID_bird6ID_uniq_set))
print(len(bird_tree_leaves_omaID_bird6ID_uniq),len(mine_bird_paper))



# In[ ]:


bird_paper_tree=bird_paper_tree_raw

bird_paper_tree.prune(mine_bird_paper)

print(len(bird_paper_tree_raw),len(bird_paper_tree))


# In[ ]:


nodes=["sss","sss","f"]


# In[ ]:


name2node = dict([ [n, None] for n in nodes if type(n) is str])


# In[ ]:


name2node


# In[ ]:


len(bird_paper_tree_nodes),len(set(bird_paper_tree_nodes))


# In[ ]:


for i in list(set(bird_paper_tree_nodes)):
    if bird_paper_tree_nodes.count(i)>3:
        a=1#print(i)


# In[ ]:




len(bird_paper_tree_nodes_raw),len(set(bird_paper_tree_nodes_raw))


# In[ ]:



for i in bird_paper_tree_nodes_raw:
    node_name_split= i.split('_')
    six_letter_name= ''.join([i[:3].upper() for i in node_name_split])
    if six_letter_name=="CHLFLA" :
        print(i)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


bird_paper_tree.write()


# In[ ]:





# In[ ]:





# In[ ]:




