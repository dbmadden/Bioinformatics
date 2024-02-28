# load libraries at the start of the script
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("DECIPHER")
library(Biostrings)
library(DECIPHER)
library(Biostrings)
library(msa) # missing
library(seqinr) # missing

# first i am setting the folder i will be working in 

setwd("/Users/iamdeanooo/Documents/GitHub/Bioinformatics/")
# then i am reading the fasta file in the folder and giving it a variable
sequences <- readDNAStringSet("sequences.fasta")

# this is singling out whichever sequence I want
sequences$Homo_sapiens_20

myAAseqs <- AAStringSet(sequences) # this line of code isn't working
myAlignment <- msa(sequences, "ClustalW") # need to pass the 'sequences' variable, not amino acids
myAlignment

# after running my alignment i can see homo_sapiens_6 has two mismatches and a deletion
# before the ... which is AT when most of the other sequences end in CA

# these lines of code weren't working so I removed them by making them into comments
# alignment <- pairwiseAlignment("sequences.fasta")
# msaIdentity("sequences.fasta")
# pairwiseAlignment("sequences.fasta")
# sequences <- pairwiseAlignment(sequences)
print(myAlignment)
# print(identi)
# this next line isn't working. Make sure to look at the help page for the function to make sure that you're using the right variable format for the function
writeXStringSet("AATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGG
                GCTGGGCATGAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAG
                CAACCTCAAACAGACACCATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCA
                AGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAA
                GGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCT
                CTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAG
                TCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAG
                TGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGA
                GCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTT
                CTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGG") # missing the quotes at the end of the string
# i am now converting my alignment to seqinr                
                
dis_align <- msaConvert(myAlignment, type = "seqinr::alignment")

# now i am running a distance alignment
dismat <- dist.alignment(dis_align)

dna_sequences <- sequences 
protein.sequence <- translate(DNAString(dna_sequences)) # missing the final parenthesis
print(protein.sequences)
dna_sequence <- ("AATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGG
                GCTGGGCATGAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAG
                CAACCTCAAACAGACACCATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCA
                AGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAA
                GGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCT
                CTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAG
                TCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAG
                TGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGA
                GCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTT
                CTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGG")


                