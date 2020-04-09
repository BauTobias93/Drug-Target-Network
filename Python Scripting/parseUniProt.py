from Bio import SeqIO

uniprot = SeqIO.index("./data/uniprot_sprot.dat", "swiss")
with open("selected.dat", "wb") as out_handle:
    #for acc in ["P33487", "P19801", "P13689", "Q8JZQ5", "Q9TRC7"]:
    for acc in ["P00533"]:
        out_handle.write(uniprot.get_raw(acc))

print("Finished parsing")