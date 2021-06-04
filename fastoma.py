
#!/usr/bin/python3


# 05/26/2021

from sys import argv


my_hog_id = argv[1]     # hog id output of OMAmer 
oma_database = argv[2]  # address  "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/OmaServer.h5"

root_level = argv[3]   # "Homininae"

import pyoma.browser.db
db = pyoma.browser.db.Database(oma_database)
db.get_release_name()



my_hog_list=db.hog_members_from_hog_id(my_hog_id,root_level)
prot_ids=[i[0] for i in my_hog_list]

og_list=[]
for prot_id in prot_ids:
   og_id=db.entry_by_entry_nr(prot_id)["OmaGroup"]
   og_list.append(og_id)

#print(og_list)
mostFrequent_og=max(set(og_list), key=og_list.count)
print(mostFrequent_og, og_list.count(mostFrequent_og), len(og_list))



