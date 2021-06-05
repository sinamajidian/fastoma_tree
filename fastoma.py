#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/python3
from sys import argv
import pyoma.browser.db as db
import pyoma.browser.models as mod
import zoo.wrappers.aligners.mafft as mafft
import zoo.wrappers.treebuilders.fasttree as fasttree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, SingleLetterAlphabet
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/OmaServer.h5"

omamer_output_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v2a/HUMAN_q2.omamer"  #v1d/omamer_out_Eukaryota.fa"  #v1f/AEGTS.omamer"
query_protein_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v2a/HUMAN_q2.fa" #v1d-omamer/query3.fa"

# The species name of query is the name of the file; for HUMAN_q.fa will be "HUMAN_q"

# oma_output_address = argv[1] 
# oma_database_address= argv[2]  # address  "../smajidi1/fastoma/archive/OmaServer.h5"


# In[2]:


################### Parsing  query protein  ########
#####################################################

query_protein_records = list(SeqIO.parse(query_protein_address, "fasta")) 
print(len(query_protein_records))


# In[3]:


################### Parsing omamer's output  ########
#####################################################
omamer_output_file = open(omamer_output_address,'r');

list_query_protein_name= []
list_query_inferred_hog= []

for line in omamer_output_file:
    line_strip=line.strip()
    if not line_strip.startswith('qs'):
        line_split= line_strip.split("\t")        
        list_query_protein_name.append(line_split[0])
        list_query_inferred_hog.append(line_split[1])
print("number of proteins in omamer output",len(list_query_inferred_hog),list_query_inferred_hog)


# In[4]:


# ############# Extracting the unique HOG list  ########
# #####################################################
# list_inferred_hog_unique = list(set(list_query_inferred_hog))


# if len(list_inferred_hog_unique) == len(list_query_inferred_hog):
#     list_inferred_hog_unique = list_query_inferred_hog
#     list_idx_query_of_uniq_hog = list(range(len(list_inferred_hog_unique)))
# else: #   remove those queries of reptead hogs
    
#     list_idx_query_of_uniq_hog=[1,3]
    
#     print(" to be coded")
# #     list_idx_query_per_uniq_hog=[]      
# #     for hog_unique in list_inferred_hog_unique:

# #     each element is not a list anymore, each element is an integer 
# #         indx_list=[i for i, x in enumerate(list_query_inferred_hog) if x == hog_unique]  # [[0], [1,2,3], ]
# #         list_idx_query_per_uniq_hog.append(indx_list)

# list_idx_query_of_uniq_hog=[1,3]

# unique_num=len(list_idx_query_of_uniq_hog)
# print(unique_num)


# In[5]:


#         list_query_protein_name
#         list_query_inferred_hog
        
#         query_protein_records
        
#         query_protein= list_query_protein_name_filtered[item_idx]

        
#         ### filtering   query recrod

        
        
list_query_inferred_hog_filtered = list_query_inferred_hog
list_query_protein_name_filtered = list_query_protein_name

query_protein_records_filtered= query_protein_records

unique_num=len(list_query_protein_name_filtered)


# In[ ]:





# In[6]:


############ Extracting the most frequent OG  ########
#####################################################

oma_db = db.Database(oma_database_address)
print("OMA data is parsed and its release name is :", oma_db.get_release_name())
mostFrequent_OG_list=[]


for  item_idx in range(unique_num):
    
    hog_id= list_query_inferred_hog_filtered[item_idx]
    hog_members = oma_db.member_of_hog_id(hog_id, level = None)                # members of the input hog_id as objects
    proteins_id_hog = [hog_member[0] for hog_member in hog_members]              # the protein IDs of hog members
    proteins_object_hog = [db.ProteinEntry(oma_db, pr) for pr in proteins_id_hog]  # the protein IDs of hog members
    OGs_correspond_proteins = [pr.oma_group for pr in proteins_object_hog]
    mostFrequent_OG = max(set(OGs_correspond_proteins), key = OGs_correspond_proteins.count)
    mostFrequent_OG_list.append(mostFrequent_OG)

    #index_query_protein=list_idx_query_of_uniq_hog[item_idx]
    #query_protein= list_query_protein_name[index_query_protein]
    query_protein= list_query_protein_name_filtered[item_idx]
    
    print("For query",query_protein )
    print("the most Frequent_OG is:",mostFrequent_OG, "with frequency of",OGs_correspond_proteins.count(mostFrequent_OG),
          "out of", len(OGs_correspond_proteins),"\n")


# In[ ]:


########## Combine proteins of OG with queries ##################
#################################################################

query_species_name = ''.join(query_protein_address.split("/")[-1].split(".")[:-1]) #'/work/fastoma/v2a/HUMAN_q.fa'


seqRecords_all=[]
for  item_idx in range(unique_num):
    mostFrequent_OG = mostFrequent_OG_list[item_idx]
    OG_members = oma_db.oma_group_members(mostFrequent_OG)
    proteins_object_OG = [db.ProteinEntry(oma_db, pr) for pr in OG_members]  # the protein IDs of og members
    seqRecords_OG= [SeqRecord(Seq(pr.sequence), str(pr.genome.uniprot_species_code), '', '') for pr in proteins_object_OG ] # covnert to biopython objects
    
    #query_protein_record_id = list_idx_query_of_uniq_hog[item_idx]
    #index_query_protein=list_idx_query_of_uniq_hog[item_idx]
    query_protein= list_query_protein_name_filtered[item_idx]
    
    print("For query index",item_idx,"name",query_protein)
    seqRecords_query =  query_protein_records_filtered[item_idx] #query_protein_records[index_query_protein]
    
    seqRecords_query_edited = SeqRecord(Seq(str(seqRecords_query.seq)), query_species_name, '', '')
    seqRecords =seqRecords_OG + [seqRecords_query_edited]
    seqRecords_all.append(seqRecords)
    
    print("length of OG",mostFrequent_OG,"was",len(seqRecords_OG),",now is",len(seqRecords),"\n")
    
print("number of OGs",len(seqRecords_all))


# In[ ]:


############## MSA  ##############
##################################

result_maf2_all=[]
for  item_idx in range(unique_num):
    seqRecords=seqRecords_all[item_idx]
        
    wrapper_maf = mafft.Mafft(seqRecords,datatype="PROTEIN")
    result_maf1 = wrapper_maf()
    time_taken_maf = wrapper_maf.elapsed_time  # 
    print("time elapsed for multiple sequence alignment: ",time_taken_maf)

    result_maf2 = wrapper_maf.result
    result_maf2_all.append(result_maf2)

# for  item_idx in range(unique_num):
#     result_maf2=result_maf2_all[item_idx]
#     out_name_msa=omamer_output_address+"_out_msa"+str(item_idx)+".txt"
#     handle_msa_fasta = open(out_name_msa,"w")
#     SeqIO.write(result_maf2, handle_msa_fasta,"fasta")
#     handle_msa_fasta.close()


# In[ ]:


############## Concatante alignments  ##############
####################################################

alignments= result_maf2_all
all_labels = set(seq.id for aln in alignments for seq in aln)
print(len(all_labels))

# Make a dictionary to store info as we go along
# (defaultdict is convenient -- asking for a missing key gives back an empty list)
concat_buf = defaultdict(list)

# Assume all alignments have same alphabet
alphabet = alignments[0]._alphabet

for aln in alignments:
    length = aln.get_alignment_length()
    # check if any labels are missing in the current alignment
    these_labels = set(rec.id for rec in aln)
    missing = all_labels - these_labels

    # if any are missing, create unknown data of the right length,
    # stuff the string representation into the concat_buf dict
    for label in missing:
        new_seq = UnknownSeq(length, alphabet=alphabet)
        concat_buf[label].append(str(new_seq))

    # else stuff the string representation into the concat_buf dict
    for rec in aln:
        concat_buf[rec.id].append(str(rec.seq))

# Stitch all the substrings together using join (most efficient way),
# and build the Biopython data structures Seq, SeqRecord and MultipleSeqAlignment
msa = MultipleSeqAlignment(SeqRecord(Seq(''.join(seq_arr), alphabet=alphabet), id=label)
                            for (label, seq_arr) in concat_buf.items())

out_name_msa=omamer_output_address+"_out_msa_concatanated.txt"
handle_msa_fasta = open(out_name_msa,"w")
SeqIO.write(msa, handle_msa_fasta,"fasta")
handle_msa_fasta.close()
    
print(len(msa),msa.get_alignment_length()) 


# In[ ]:


############## Tree inference  ###################
##################################################

wrapper_tree=fasttree.Fasttree(msa,datatype="PROTEIN")
result_tree1 = wrapper_tree()

time_taken_tree = wrapper_tree.elapsed_time 
time_taken_tree

result_tree2 = wrapper_tree.result
tree_nwk=str(result_tree2["tree"])
print(len(tree_nwk))



out_name_tree=omamer_output_address+"_tree_02.txt"
file1 = open(out_name_tree,"w")
file1.write(tree_nwk)
file1.close() 


# In[ ]:


tree_nwk
