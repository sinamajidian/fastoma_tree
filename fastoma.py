#!/usr/bin/env python
# coding: utf-8

# In[31]:


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

#my_hog_id = argv[1]     # hog id output of OMAmer 
#oma_database = argv[2]  # address  "../smajidi1/fastoma/archive/OmaServer.h5"
# query_protein= argv[3] # 


# In[32]:


oma_output_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v1d-omamer/omamer_out_Eukaryota_que3.fa"
oma_output_file = open(oma_output_address,'r');

list_query_protein_name= []
list_query_inferred_hog= []

for line in oma_output_file:
    line_strip=line.strip()
    if not line_strip.startswith('qs'):
        line_split= line_strip.split("\t")        
        list_query_protein_name.append(line_split[0])
        list_query_inferred_hog.append(line_split[1])


# In[33]:


#print(list_query_inferred_hog)
list_inferred_hog_unique = list(set(list_query_inferred_hog))
list_idx_query_per_uniq_hog=[]

for hog_unique in list_inferred_hog_unique:
    indx_list=[i for i, x in enumerate(list_query_inferred_hog) if x == hog_unique]
    list_idx_query_per_uniq_hog.append(indx_list)
#print(list_idx_query_per_uniq_hog)


# In[34]:


oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/OmaServer.h5"
oma_db = db.Database(oma_database_address)
print("OMA data is parsed and its release name is :", oma_db.get_release_name())
mostFrequent_OG_list=[]
unique_num=len(list_inferred_hog_unique)

for  item_idx in range(unique_num):
    
    hog_id= list_inferred_hog_unique[item_idx]
    hog_members = oma_db.member_of_hog_id(hog_id, level = None)                # members of the input hog_id as objects
    proteins_id_hog = [hog_member[0] for hog_member in hog_members]              # the protein IDs of hog members
    proteins_object_hog = [db.ProteinEntry(oma_db, pr) for pr in proteins_id_hog]  # the protein IDs of hog members
    OGs_correspond_proteins = [pr.oma_group for pr in proteins_object_hog]
    mostFrequent_OG = max(set(OGs_correspond_proteins), key = OGs_correspond_proteins.count)
    mostFrequent_OG_list.append(mostFrequent_OG)
    
    print("For queries",[list_query_protein_name[i]  for i in list_idx_query_per_uniq_hog[item_idx]])
    print("the most Frequent_OG is:",mostFrequent_OG, "with frequency of",OGs_correspond_proteins.count(mostFrequent_OG),
          "out of", len(OGs_correspond_proteins),"\n")


# In[35]:


query_protein= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v1d-omamer/query3.fa"
query_protein_records = list(SeqIO.parse(query_protein, "fasta")) 
print(len(query_protein_records))


# In[66]:


seqRecords_all=[]
for  item_idx in range(unique_num):
    mostFrequent_OG = mostFrequent_OG_list[item_idx]

    OG_members = oma_db.oma_group_members(mostFrequent_OG)
    proteins_object_OG = [db.ProteinEntry(oma_db, pr) for pr in OG_members]  # the protein IDs of og members
    seqRecords_OG= [SeqRecord(Seq(pr.sequence), str(pr.omaid), '', '') for pr in proteins_object_OG ] # covnert to biopython objects
    # genome.uniprot_species_code
    seqRecords_queries = []
    for i in  list_idx_query_per_uniq_hog[item_idx]:
        seqRecords_queries.append(query_protein_records[i])
        
    seqRecords =seqRecords_OG + seqRecords_queries
    print(mostFrequent_OG,len(seqRecords),len(seqRecords_OG))
    seqRecords_all.append(seqRecords)


# In[67]:


############## MSA  ##############
##################################

result_maf2_all=[]
for  item_idx in range(unique_num):
    seqRecords=seqRecords_all[item_idx]
        
    wrapper_maf = mafft.Mafft(seqRecords)
    result_maf1 = wrapper_maf()
    time_taken_maf = wrapper_maf.elapsed_time  # 
    print("time elapsed for multiple sequence alignment: ",time_taken_maf)

    result_maf2 = wrapper_maf.result
    result_maf2_all.append(result_maf2)

for  item_idx in range(unique_num):
    result_maf2=result_maf2_all[item_idx]
    out_name_msa="out_msa_"+str(item_idx)+".txt"
    handle_msa_fasta = open(out_name_msa,"w")
    SeqIO.write(result_maf2, handle_msa_fasta,"fasta")
    handle_msa_fasta.close()


# In[68]:


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
print(len(msa),msa.get_alignment_length()) 


# In[63]:


############## Tree inference  ###################
##################################################

wrapper_tree=fasttree.Fasttree(msa)
result_tree1 = wrapper_tree()

time_taken_tree = wrapper_tree.elapsed_time  # 
time_taken_tree

result_tree2 = wrapper_tree.result
tree_nwk=str(result_tree2["tree"])
print(len(tree_nwk))



out_name_tree="tree_01.txt"
file1 = open(out_name_tree,"w")
file1.write(tree_nwk)
file1.close() 


# In[ ]:




