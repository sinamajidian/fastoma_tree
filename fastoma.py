
from sys import argv
import pyoma.browser.db as db
import pyoma.browser.models as mod
import zoo.wrappers.aligners.mafft as mafft
import zoo.wrappers.treebuilders.fasttree as fasttree

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

hog_id = "HOG:A0563565.5a.6a.2b"
oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/OmaServer.h5"
query_protein= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/v1c/HUMAN91530.fa"

oma_db = db.Database(oma_database_address)
print("OMA data is parsed and its release name is :", oma_db.get_release_name())

hog_members = oma_db.member_of_hog_id(hog_id, level = None)                # members of the input hog_id as objects
proteins_id_hog = [hog_member[0] for hog_member in hog_members]              # the protein IDs of hog members

proteins_object_hog = [db.ProteinEntry(oma_db, pr) for pr in proteins_id_hog]  # the protein IDs of hog members

OGs_correspond_proteins = [pr.oma_group for pr in proteins_object_hog]

mostFrequent_OG = max(set(OGs_correspond_proteins), key = OGs_correspond_proteins.count)

print(mostFrequent_OG, OGs_correspond_proteins.count(mostFrequent_OG), len(OGs_correspond_proteins))



OG_members = oma_db.oma_group_members(mostFrequent_OG)
proteins_object_OG = [db.ProteinEntry(oma_db, pr) for pr in OG_members]  # the protein IDs of hog members

#seqRecords_raw= [(pr.sequence, pr.entry_nr) for pr in proteins_object_OG ]
seqRecords_OG= [SeqRecord(Seq(pr.sequence), str(pr.entry_nr), '', '') for pr in proteins_object_OG ] # covnert to biopython objects


query_protein_record = SeqIO.read(query_protein, "fasta")

query_prot = str(query_protein_record.seq)
seqRecords =seqRecords_OG+ [SeqRecord(Seq(query_prot), 'query', '', '') ]

print(len(seqRecords))


wrapper_maf = mafft.Mafft(seqRecords)
result_maf1 = wrapper_maf()
time_taken_maf = wrapper_maf.elapsed_time  # 
print("time elapsed for multiple sequence alignment: ",time_taken_maf)



result_maf2 = wrapper_maf.result

handle_msa_fasta = open("out_msa.txt","w")
SeqIO.write(result_maf2, handle_msa_fasta,"fasta")
handle_msa_fasta.close()


wrapper_tree=fasttree.Fasttree(result_maf2)
result_tree1 = wrapper_tree()
time_taken_tree = wrapper_tree.elapsed_time  # 
print("time elapsed for tree inference: ", time_taken_tree)

result_tree2 = wrapper_tree.result
tree_nwk=str(result_tree2["tree"])
file1 = open("tree2.txt","w")
file1.write(tree_nwk)
file1.close() #to change file access modes
  

