genome_size = 267

for MB in range(genome_size):
    if MB % 3 == 0:
        print("indpalmata.genomes.addNewMutation(m2, 0.0, ", MB, ":", MB, ");", sep="")