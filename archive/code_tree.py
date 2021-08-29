#!/usr/bin/env python
# coding: utf-8

# In[80]:


#!/usr/bin/python3
import numpy as np
from sys import argv
from datetime import datetime
#import concurrent.futures
# import ast
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ete3
import dendropy
from  dendropy import Tree
datasets_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/"
# oma_database_address = datasets_address + "OmaServer.h5"
omaID_address = datasets_address+"oma-species.txt"
bird6ID_address = datasets_address+"info.tsv"

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

(taxonID_omaID,taxonID_bird6ID,omaID_taxonID,bird6ID_taxonID) = read_taxonID_map(omaID_address,bird6ID_address)
# left is the ky
# taxonID_omaID  {taxonID:omaID, .. }
# taxonID_omaID[1111]=ABABA
taxonID_map = {**taxonID_omaID,**taxonID_bird6ID }
print(len(taxonID_omaID),len(taxonID_bird6ID),len(taxonID_map))
#taxonID_list= list(taxonID_map.keys())


# # Reading FastOMA bird  tree

# In[81]:



project_folder_root="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3b/" # v3b/
project_folder = project_folder_root +"hogmapX/"
fast_tree_address=project_folder_root+"iqtree5/_100_msa_concatanated.txt_copy1.contree"  # _collapsed_brnch_0.01# 1e-05" # 
#"iqtree5/_100_msa_concatanated.txt_copy1.contree"
#"iqtree4_1/_100_msa_concatanated.txt.contree"
#"iqtree5/_100_msa_concatanated.txt_copy1.contree"
#"iqtree11/border_out_parition.nex.contree"
#_100_msa_concatanated.txt_copy1.contree" # _collapsed_95


fast_tree_address 


# In[82]:


fast_tree_raw= ete3.Tree(fast_tree_address) # ,format=0
print("len fast tree",len(fast_tree_raw))
fast_tree_leaves_omaID_bird6ID=[]
for node in fast_tree_raw.traverse(strategy="postorder"):
    if node.is_leaf() : # why ?
        fast_tree_leaves_omaID_bird6ID.append(node.name)
print("bird_hog_tree_leaf_count",len(fast_tree_leaves_omaID_bird6ID))
fast_tree_leaves_omaID_bird6ID[:3]

fast_tree_leaves_taxonID = []
for i3 in fast_tree_leaves_omaID_bird6ID:
    if i3 in omaID_taxonID:
        taxonID=omaID_taxonID[i3]
    if i3  in bird6ID_taxonID:
        taxonID= bird6ID_taxonID[i3]
    fast_tree_leaves_taxonID.append(taxonID) # if it was in the both I save the bird ID    

fast_tree_leaves_taxonID_unq=list(set(fast_tree_leaves_taxonID))
print(len(fast_tree_leaves_taxonID),len(fast_tree_leaves_taxonID_unq))


both_bird_oma_taxonID=[]
for i in fast_tree_leaves_taxonID_unq:
    if fast_tree_leaves_taxonID.count(i)>1:
        both_bird_oma_taxonID.append(i)
print(both_bird_oma_taxonID)
both_bird_oma_omaID=[]
for i in both_bird_oma_taxonID:
    both_bird_oma_omaID.append(taxonID_omaID[i])
print(both_bird_oma_omaID)

fast_tree_leaves_omaID_bird6ID_uniq=[i for i in fast_tree_leaves_omaID_bird6ID if i not in both_bird_oma_omaID]
print(len(fast_tree_leaves_omaID_bird6ID),len(fast_tree_leaves_omaID_bird6ID_uniq))
fast_tree= ete3.Tree(fast_tree_address)
fast_tree.prune(fast_tree_leaves_omaID_bird6ID_uniq)
print(len(fast_tree_raw),len(fast_tree))


# # Reading  NCBI tree

# In[83]:


ncbi = ete3.NCBITaxa()  # first time download in ~/.etetoolkit/
ncbi_sub_tree = ncbi.get_topology(fast_tree_leaves_taxonID_unq)
## change the node name from NCBI taxon id (integer) to   omaID (5-letter) or bird6ID (6-letter)
for node in ncbi_sub_tree.traverse(strategy="postorder"):
    node.name2 = node.name
    if node.is_leaf() : # why ? and int(node.name) in taxonID_map
        node.name = taxonID_map[int(node.name)]

print("The NCBI taxanomy is read and the leaves names changed to OMA/bird6 ID containing")
print(len(ncbi_sub_tree)) 
#ncbi_sub_tree.write() 


# # Reading bird paper tree

# In[84]:


bird_paper_tree_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/bird_paper_tree.txt"
bird_paper_tree_raw= ete3.Tree(bird_paper_tree_address, format=1)
#len(bird_paper_tree_raw)

bird_paper_tree_leaves=[]
for node in bird_paper_tree_raw.traverse(strategy="postorder"):
    if node.is_leaf() : # why ?
        bird_paper_tree_leaves.append(node.name)
print(len(bird_paper_tree_leaves))
bird_paper_tree_leaves[:3]

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
#bird_paper_tree_raw.prune(bird_SCINAME_list)

notin_paper_tree= []
for i in bird_SCINAME_list:
    if i not in bird_paper_tree_leaves: 
        notin_paper_tree.append(i)
print(notin_paper_tree)
#notin_paper_tree=["Cercotrichas_coryphaeus","Corvus_cornix","Eolophus_roseicapillus","Nannopterum_auritus","Nannopterum_brasilianus","Nannopterum_harrisi","Urile_pelagicus"]

bird_SCINAME_list_filt=[i for i in bird_SCINAME_list if i not in notin_paper_tree]
print(len(bird_SCINAME_list),len(bird_SCINAME_list_filt) )    
bird_paper_tree= ete3.Tree(bird_paper_tree_address, format=1)
bird_paper_tree.prune(bird_SCINAME_list_filt)


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


# # Comapring intersection all trees

# In[85]:


ncbi_sub_tree_nodes =[]
for node in ncbi_sub_tree.traverse(strategy="postorder"):
    if node.is_leaf() : 
        ncbi_sub_tree_nodes.append(node.name)
len(ncbi_sub_tree_nodes),len(set(ncbi_sub_tree_nodes+bird_paper_tree_nodes))

fast_tree_leaves_omaID_bird6ID_uniq_set=set(fast_tree_leaves_omaID_bird6ID_uniq)
bird_paper_tree_nodes_set=set(bird_paper_tree_nodes)
fast_tree_paper_tree_intersection= bird_paper_tree_nodes_set.intersection(fast_tree_leaves_omaID_bird6ID_uniq_set)
print(len(fast_tree_paper_tree_intersection))

ncbi_fast_tree_paper_tree_intersection= fast_tree_paper_tree_intersection.intersection(set(ncbi_sub_tree_nodes))
print(len(ncbi_fast_tree_paper_tree_intersection))
ncbi_sub_tree.prune(ncbi_fast_tree_paper_tree_intersection) # fast_tree_paper_tree_intersection 

fast_tree.prune(fast_tree_paper_tree_intersection)
bird_paper_tree.prune(fast_tree_paper_tree_intersection) # fast_tree_paper_tree_intersection 
len(fast_tree),len(bird_paper_tree), len(ncbi_sub_tree)


# In[108]:


# # out_folder="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3b/"
# # #node.name

# fast_tree.write(format=1, outfile=out_folder+"tree_fast_tree_intersection.nw")
# bird_paper_tree.write(format=1, outfile=out_folder+"tree_bird_paper_tree_intersection.nw")
# ncbi_sub_tree.write(format=1, outfile=out_folder+"tree_ncbi_sub_tree_intersection.nw")


# In[86]:


fast_tree.set_outgroup("CALPUG") # CALPUG MOUS  BUCCAP

out_paper_ncbi = bird_paper_tree.robinson_foulds(ncbi_sub_tree) # , expand_polytomies = True  #,polytomy_size_limit=20  ,unrooted_trees=True # , expand_polytomies = True
(rf_paper_ncbi, maxparts_paper_ncbi, common_attrs, partitions_paper, partitions_ncbi, discard_t1, discard_t2)=out_paper_ncbi 
print("RF distance is %s over a total of %s" %(rf_paper_ncbi, maxparts_paper_ncbi))
print("Partitions: ")
print("in paper:",len(partitions_paper),", in ncbi:",len(partitions_ncbi)) # , "Common Partitions",len(common_attrs_paper),)
print("in both paper and ncbi:", len(partitions_ncbi & partitions_paper))
print("only in ncbi, not in paper:", len(partitions_ncbi-partitions_paper)) # order not sure  http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#robinson-foulds-distance
print("only in paper, not in ncbi:", len(partitions_paper-partitions_ncbi))


print("-"*15)
out_fast_ncbi = fast_tree.robinson_foulds(ncbi_sub_tree)  # ,unrooted_trees=True # , expand_polytomies = True
(rf_fast_ncbi, maxparts_fast_ncbi, common_attrs, partitions_fast, partitions_ncbi2, discard_t1, discard_t2)=out_fast_ncbi 
print("RF distance is %s over a total of %s" %(rf_fast_ncbi, maxparts_fast_ncbi))
print("Partitions: ")
print("in fast:",len(partitions_fast),", in ncbi:",len(partitions_ncbi2)) 
print("in both fast and ncbi:", len(partitions_fast & partitions_ncbi2 ))
print("only in ncbi, not in fast:", len(partitions_ncbi2 - partitions_fast))
print("only in fast, not in ncbi:", len(partitions_fast- partitions_ncbi2))

assert len(partitions_ncbi2-partitions_ncbi)==0

print("-"*15)
out_fast_paper = fast_tree.robinson_foulds(bird_paper_tree)  # ,unrooted_trees=True # , expand_polytomies = True
(rf_fast_paper, maxparts_fast_paper, common_attrs, partitions_fast, partitions_paper, discard_t1, discard_t2)=out_fast_paper
print("RF distance is %s over a total of %s" %(rf_fast_ncbi, maxparts_fast_ncbi))
print("Partitions: ")
print("in fast:",len(partitions_fast),", in ncbi:",len(partitions_paper)) 
print("in both fast and ncbi:", len(partitions_fast & partitions_paper))
print("only in ncbi, not in fast:", len(partitions_paper - partitions_fast))
print("only in fast, not in ncbi:", len(partitions_fast- partitions_paper))


# In[87]:


print("only in fastOMA:", len((partitions_fast-partitions_ncbi2)-partitions_paper))
print("only in fastOMA & paper:",len((partitions_fast&partitions_paper)-partitions_ncbi2))
print("only in fastOMA & NCBI:",len((partitions_fast&partitions_ncbi2)- partitions_paper))
print("only in paper & NCBI:",len((partitions_paper&partitions_ncbi2)-  partitions_fast))
print("only in fastOMA & NCBI & paper:",len((partitions_fast&partitions_ncbi2)& partitions_paper))
print("only in paper:",len((partitions_paper-partitions_fast)-partitions_ncbi2 ))
print("only in NCBI:",len((partitions_paper-partitions_fast)-partitions_ncbi2 ))


# # plot different collapsing 

# In[ ]:



# val=[len((edges_fast-edges_ncbi2)-edges_paper),len((edges_fast&edges_paper)-edges_ncbi2),
#      len((edges_fast&edges_ncbi2)- edges_paper),len((edges_paper&edges_ncbi2)-  edges_fast),
#      len((edges_fast&edges_ncbi2)& edges_paper),
#      len((edges_paper-edges_fast)-edges_ncbi2 ),len((edges_paper-edges_fast)-edges_ncbi2) ]
# print(val)
# # branch lenght
# t_all=[ [131, 135, 8, 6, 433, 124, 124], # no
#        [131, 135, 8, 6, 433, 124, 124],  # 0.00001
#        [131, 135, 8, 6, 433, 124, 124],  # 0.0001
#        [129, 135, 8, 6, 433, 124, 124],  # 0.001
#        [53, 108, 5, 9, 430, 151, 151],   # 0.005
#         [22, 69, 5, 17, 422, 190, 190],  #0.01
#         [1, 0, 0, 80, 359, 259, 259]]  #0.1
# ind = ["no",0.00001,0.0001,0.001,0.005,0.01,0.1]
# # ind = [0,60,80,90,95,100] # ,"BirdPaper"
# # # support values
# # t_all=[[131, 135, 8, 6, 433, 124, 124], # np
# #        [122, 130, 8, 6, 433, 129, 129], # 60
# #        [111, 129, 8, 6, 433, 130, 130], # 80
# #        [92, 127, 8, 6, 433, 132, 132],# 90
# #        [80, 124, 6, 7, 432, 135, 135], # 95
# #        [63, 120, 6, 7, 432, 139, 139]]

# import numpy as np
# t_all=np.array(t_all)
# t_all_=np.transpose(t_all)
# len(t_all_)

# # t_all=[[131, 135, 8, 6, 433, 124, 124], # np
# #        [122, 130, 8, 6, 433, 129, 129], # 60
# #        [111, 129, 8, 6, 433, 130, 130], # 80
# #        [92, 127, 8, 6, 433, 132, 132],# 90
# #        [80, 124, 6, 7, 432, 135, 135], # 95
# #        [63, 120, 6, 7, 432, 139, 139]]
# lenged1=['fastOMA','fastOMA & paper',  'fastOMA & NCBI',
#   'paper & NCBI', 'fastOMA & NCBI & paper',  'paper', 'NCBI']

# import pandas as pd
# df = pd.DataFrame(data=t_all)
# df.columns=lenged1   #ind = [0,60,80,90,95,100] # ,"BirdPaper"
# df.index=ind
# df


# In[ ]:


# ind=[str(i) for i in ind]
# width = 0.20
# palet=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
# fig = plt.figure()
# ax = fig.add_axes([0,0,1,1])
# ax.bar(ind, t_all_[0], width, color=palet[0])
# ax.bar(ind, t_all_[1], width,bottom=t_all_[0], color=palet[1])
# ax.bar(ind, t_all_[2], width,bottom=t_all_[0]+t_all_[1], color=palet[2])
# ax.bar(ind, t_all_[3], width,bottom=t_all_[0]+t_all_[1]+t_all_[2], color=palet[3])
# ax.bar(ind, t_all_[4], width,bottom=t_all_[0]+t_all_[1]+t_all_[2]+t_all_[3], color=palet[4])
# ax.bar(ind, t_all_[5], width,bottom=t_all_[0]+t_all_[1]+t_all_[2]+t_all_[3]+t_all_[4], color=palet[5])
# ax.bar(ind, t_all_[6], width,bottom=t_all_[0]+t_all_[1]+t_all_[2]+t_all_[3]+t_all_[4]+t_all_[5], color=palet[6])

# ax.set_xlabel('Support value threshold')
# ax.set_ylabel('Number of partitions')
# #ax.set_title('Collapase edges with support value < threshold ')
# ax.legend(labels=lenged1,loc=3)

# plt.show()
# plt.savefig(project_folder_root+"collapse2.pdf")


# In[ ]:


ind=["only in fastOMA","only in fastOMA & paper",
      "only in fastOMA & NCBI","only in fastOMA & NCBI & paper",
      "only in paper","only in NCBI"]
ind


# # collapsing low support values

# In[ ]:


# dendropy

# '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3b/iqtree5/collpased.tree_collapsed'
print(bi_tree_address)
tree = dendropy.Tree.get_from_path(fast_treee_address,"newick") # newick nexus
print(len(tree))
support_treshold=95
for idx, node in enumerate(tree):
    if not node.is_leaf() and node.label:
        support_value=int(node.label)
        if support_value<support_treshold:
            edge= node.edge   #print(edge.head_node.child_nodes())
            edge.collapse()
print(len(tree))

fast_treee_address_out=fast_treee_address+"_collapsed2_"+str(support_treshold)
tree.write_to_path(fast_treee_address_out, "newick")
print(fast_treee_address_out)


# # collapsing branch length

# In[ ]:


# dendropy

print(fast_tree_address)
fast_tree_address="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3b/iqtree5/_100_msa_concatanated.txt_copy1.contree"
print(fast_tree_address)

tree = dendropy.Tree.get_from_path(fast_tree_address,"newick") # newick nexus
print(len(tree))
len_all=[]
threshold=0.01
#for idx, node in enumerate(tree):
for e in tree.postorder_edge_iter():    
#     len_edg=e.length
#     if len_edg:
#         len_all.append(len_edg)
    if e.length is None or (e.length <= threshold) and e.is_internal():
        e.collapse()


# In[ ]:


# dendropy

fast_treee_address_out=fast_treee_address+"_collapsed_brnch_"+str(threshold)+"_"
tree.write_to_path(fast_treee_address_out, "newick")
print(fast_treee_address_out)


# In[ ]:


plt.hist(len_all,bins=500)
plt.show()


# # PNAS paper

# In[ ]:


project_folder_root="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/" # v3b/
pnas_tree_address=project_folder_root+"tree_pnas.1813206116.sd02.txt"  # _collapsed_brnch_0.01# 1e-05" # 

#tree2 = dendropy.Tree.get(file=open(pnas_tree_address, "r"), schema="nexus")
tree2 = dendropy.TreeList.get(path=pnas_tree_address, schema='nexus')
print(tree2[0].as_ascii_plot())
#print(tree2[0].as_ascii_plot())
print(tree2[10])


# # Standard OMA  -slow

# In[12]:



slow_tree_address="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v4a/iqtree4_1/_100_msa_concatanated.txt.contree"
slow_tree= ete3.Tree(slow_tree_address) # ,format=0
print("len slow tree",len(slow_tree))
slow_tree.prune(fast_tree_paper_tree_intersection)
print(len(slow_tree))


# In[38]:


slow_tree.set_outgroup("CALPUG") # CALPUG MOUS  BUCCAP

print("-"*15)
out_slow_ncbi = slow_tree.robinson_foulds(ncbi_sub_tree)  # ,unrooted_trees=True # , expand_polytomies = True
(rf_slow_ncbi, maxparts_slow_ncbi, common_attrs, partitions_slow, partitions_ncbi, discard_t1, discard_t2)=out_slow_ncbi

print("RF distance is %s over a total of %s" %(rf_slow_ncbi, maxparts_slow_ncbi))
print("Partitions: ")
print("in ncbi:",len(partitions_ncbi),", in slow:",len(partitions_slow)) 
print("in both ncbi and slow:", len(partitions_ncbi & partitions_slow))
print("only in ncbi, not in slow:", len(partitions_ncbi - partitions_slow))
print("only in slow, not in ncbi:", len(partitions_slow- partitions_ncbi))


# In[19]:


slow_tree.set_outgroup("CALPUG") # CALPUG MOUS  BUCCAP

print("-"*15)
out_slow_fast = slow_tree.robinson_foulds(fast_tree)  # ,unrooted_trees=True # , expand_polytomies = True
(rf_slow_fast, maxparts_slow_fast, common_attrs, partitions_slow, partitions_fast, discard_t1, discard_t2)=out_slow_fast

print("RF distance is %s over a total of %s" %(rf_slow_fast, maxparts_slow_fast))
print("Partitions: ")
print("in fast:",len(partitions_fast),", in slow:",len(partitions_slow)) 
print("in both fast and slow:", len(partitions_fast & partitions_slow))
print("only in fast, not in slow:", len(partitions_fast - partitions_slow))
print("only in slow, not in fast:", len(partitions_slow- partitions_fast))


# # coloring based on partitions differences in ncbi/paper/fast

# In[40]:





# In[103]:


#partitions_paper-partitions_fast)-partitions_ncbi2
part_only_ncbi= list((partitions_ncbi2 - partitions_fast)-partitions_paper)

part_only_ncbi_ancestor=[]
for i in part_only_ncbi:
    ancestor = ncbi_sub_tree.get_common_ancestor(i)
    part_only_ncbi_ancestor.append(ancestor.name)    
print(len(part_only_ncbi_ancestor),part_only_ncbi_ancestor[:2])

#ancestor
#ancestor.children, ancestor.children[0].children

part_ncbi_fast_notpaper= list((partitions_ncbi2 & partitions_fast)-partitions_paper)
part_ncbi_fast_notpaper_ancestor=[]
for i in part_ncbi_fast_notpaper:
    ancestor = ncbi_sub_tree.get_common_ancestor(i)
    part_ncbi_fast_notpaper_ancestor.append(ancestor.name)    
print(len(part_ncbi_fast_notpaper_ancestor),part_ncbi_fast_notpaper_ancestor[:2])


part_ncbi_fast_paper= list((partitions_ncbi2 & partitions_fast)& partitions_paper)
part_ncbi_fast_paper_ancestor=[]
for i in part_ncbi_fast_paper:
    ancestor = ncbi_sub_tree.get_common_ancestor(i)
    part_ncbi_fast_paper_ancestor.append(ancestor.name)    
print(len(part_ncbi_fast_paper_ancestor),part_ncbi_fast_paper_ancestor[:2])



part_ncbi_notfast_paper= list((partitions_ncbi2 & partitions_paper)- partitions_fast)
part_ncbi_notfast_paper_ancestor=[]
for i in part_ncbi_notfast_paper:
    ancestor = ncbi_sub_tree.get_common_ancestor(i)
    part_ncbi_notfast_paper_ancestor.append(ancestor.name)    
print(len(part_ncbi_notfast_paper_ancestor),part_ncbi_notfast_paper_ancestor[:2])


# In[ ]:





# In[104]:



ncbi_sub_tree_part= ncbi_sub_tree
for node in ncbi_sub_tree.traverse(strategy="postorder"):
    
    if not node.is_leaf() : # why ? and int(node.name) in taxonID_map
        if node.name in part_only_ncbi_ancestor:
            node.support = .2
        elif node.name in part_ncbi_fast_notpaper_ancestor:
            node.support = .6
        elif node.name in part_ncbi_notfast_paper_ancestor:
            node.support = .4
        elif node.name in part_ncbi_fast_paper_ancestor:
            node.support = .8
        else:
            node.support = 0
            print(node.children)
        


# In[106]:


out_folder="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v3b/"
#node.name

ncbi_sub_tree_part.write(format=1, outfile=out_folder+"ncbi_colored.nw")
ncbi_sub_tree_part.write(format=1, outfile=out_folder+"ncbi_colored.nw")


# In[ ]:





# In[ ]:





# In[ ]:


# write ete3 tree into file
#Et.write(format=1, outfile="new_tree.nw")

#tree2 = dendropy.Tree.get_from_path("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/ncbi.tree","newick") # newick nexus
print(len(tree2))


# In[ ]:


ancestor = t.get_common_ancestor("C", "J", "B")


# In[ ]:




#edges_fast-edges_ncbi2

all1=[]
support_treshold=95
for idx, node in enumerate(tree2):
#     if node.is_leaf():
#         node1=node
        
    if not node.is_leaf() and node.label:
        a=[str(i.taxon)[1:-1] for i in node.adjacent_nodes()]

        all1+=a
        for i in set_all:
            if i in a:
                node.label=str(10)
                #print(i)
                


# In[ ]:


len(set(all1)),len(set_all)


# In[ ]:


#set_all.intersection(set(all1))


# In[ ]:





# In[ ]:





# In[ ]:


print(tree2)


# In[ ]:


node2.taxon


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


#ncbi_sub_tree.write()


# In[ ]:


out_fast_ncbi = fast_tree.robinson_foulds(ncbi_sub_tree)  # ,unrooted_trees=True # , expand_polytomies = True
(rf_fast_ncbi, maxparts_fast_ncbi, common_attrs, edges_fast, edges_ncbi2, discard_t1, discard_t2)=out_fast_ncbi 

print("RF distance is %s over a total of %s" %(rf_fast_ncbi, maxparts_fast_ncbi))
print("Partitions: ")
print("in fast:",len(edges_fast),". in ncbi:",len(edges_ncbi2)) 
print("in both fast and ncbi:", len(edges_fast &edges_ncbi2 ))
print("only in ncbi, not in fast:", len(edges_fast-edges_ncbi2))
print("only in fast, not in ncbi:", len(edges_ncbi2-edges_fast))


# In[ ]:


set1=edges_fast-edges_ncbi2
set1_list=list(set1)
print(set1_list[1])


# In[ ]:


ancestor = t.get_common_ancestor("C", "J", "B")


# In[ ]:





# In[ ]:


print(set1_list[1])


# In[ ]:



for i in set1:
    print(len(i))


# In[ ]:





# In[ ]:


i


# In[ ]:


ncbi.mrca( each partiito) 


# In[ ]:





# In[ ]:





# In[ ]:




