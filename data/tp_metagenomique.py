import gzip
with gzip.open("amplicon.fasta.gz", "rt") as f:
    for line in f:
        print(line.strip())