#!/usr/bin/env python

from sys import argv
import pyoma.browser.db as db
import pyoma.browser.models as mod
import zoo.wrappers.aligners.mafft as mafft
import zoo.wrappers.treebuilders.fasttree as fasttree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#my_hog_id = argv[1]     # hog id output of OMAmer 
#oma_database = argv[2]  # address  "../smajidi1/fastoma/archive/OmaServer.h5"
# query_protein= argv[3] # 

oma_output_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v1e/omamer_out_Eukaryota_que3.fa"

oma_output_file = open(oma_output_address,'r');

list_query_protein_name= []
list_query_inferred_hog= []

for line in oma_output_file:
    line_strip=line.strip()
    if not line_strip.startswith('qs'):
        line_split= line_strip.split("\t")
        
        list_query_protein_name.append(line_split[0])
        list_query_inferred_hog.append(line_split[1])

print(list_query_inferred_hog)

list_inferred_hog_unique = list(set(list_query_inferred_hog))
#list_inferred_hog_unique
list_idx_query_per_uniq_hog=[]

for hog_unique in list_inferred_hog_unique:
    indx_list=[i for i, x in enumerate(list_query_inferred_hog) if x == hog_unique]
    list_idx_query_per_uniq_hog.append(indx_list)

print(list_idx_query_per_uniq_hog)

oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/OmaServer.h5"
query_protein= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v1d/HUMAN48905.fa"

oma_db = db.Database(oma_database_address)
print("OMA data is parsed and its release name is :", oma_db.get_release_name())

#hog_id = "HOG:A0563565.5a.6a.2b"
mostFrequent_OG_list=[]

for  hog_idx, hog_id  in enumerate(list_inferred_hog_unique):

    hog_members = oma_db.member_of_hog_id(hog_id, level = None)                # members of the input hog_id as objects
    proteins_id_hog = [hog_member[0] for hog_member in hog_members]              # the protein IDs of hog members

    proteins_object_hog = [db.ProteinEntry(oma_db, pr) for pr in proteins_id_hog]  # the protein IDs of hog members

    OGs_correspond_proteins = [pr.oma_group for pr in proteins_object_hog]

    mostFrequent_OG = max(set(OGs_correspond_proteins), key = OGs_correspond_proteins.count)

    mostFrequent_OG_list.append(mostFrequent_OG)
    
    print("For queries",[list_query_protein_name[i]  for i in list_idx_query_per_uniq_hog[hog_idx]])
    print("the most Frequent_OG is:",mostFrequent_OG, "with frequency of",OGs_correspond_proteins.count(mostFrequent_OG),
          "out of", len(OGs_correspond_proteins),"\n")

for mostFrequent_OG in mostFrequent_OG_list:
    OG_members = oma_db.oma_group_members(mostFrequent_OG)
    proteins_object_OG = [db.ProteinEntry(oma_db, pr) for pr in OG_members]  # the protein IDs of og members
    
    seqRecords_OG= [SeqRecord(Seq(pr.sequence), str(pr.omaid), '', '') for pr in proteins_object_OG ] # covnert to biopython objects

query_protein_record = SeqIO.read(query_protein, "fasta")
#seqRecords=seqRecords_OG+ [query_protein_record]

query_prot = str(query_protein_record.seq)
seqRecords =seqRecords_OG+ [SeqRecord(Seq(query_prot), 'query', '', '') ]

print(len(seqRecords))

#for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
#    print(seq_record.id)
#    print(repr(seq_record.seq))
#    print(len(seq_record))

wrapper_maf = mafft.Mafft(seqRecords)
result_maf1 = wrapper_maf()
time_taken_maf = wrapper_maf.elapsed_time  # 
print("time elapsed for multiple sequence alignment: ",time_taken_maf)


result_maf2 = wrapper_maf.result
#result_again
#str(result_again)
handle_msa_fasta = open("out_msa_2.txt","w")
SeqIO.write(result_maf2, handle_msa_fasta,"fasta")
handle_msa_fasta.close()

wrapper_tree=fasttree.Fasttree(result_maf2)
result_tree1 = wrapper_tree()
time_taken_tree = wrapper_tree.elapsed_time  # 
print("time elapsed for tree inference: ", time_taken_tree)

result_tree2 = wrapper_tree.result
tree_nwk=str(result_tree2["tree"])
file1 = open("tree_2.txt","w")
file1.write(tree_nwk)
file1.close() #to change file access modes
  