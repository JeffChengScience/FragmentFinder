# Searching Private Transposon

library("Biostrings") 

#Generate the all the different concatonating fragments sizes
x= "ggcttctctcatgagaagtcttttttatttaaaataaatataaaataaaatagaggctataaatagcctctattttatgtgagaaatccctaaataaaaagatgccagtgtgctggtaacaggttggctgataagtccccggtctgacagaagcaaacttaagagtgtgttgatagtgcagtatcttaaaattttgtataataggaattgaagttaaattagatgctaaaaatttgtaattaagaaggagtgattacatgaacaaaaatataaaatattctcaaaactttttaacgagtgaaaaagtactcaaccaaataataaaacaattgaatttaaaagaaaccgataccgtttacgaaattggaacaggtaaagggcatttaacgacgaaactggctaaaataagtaaacaggtaacgtctattgaattagacagtcatctattcaacttatcgtcagaaaaattaaaactgaatactcgtgtcactttaattcaccaagatattctacagtttcaattccctaacaaacagaggtataaaattgttgggagtattccttaccatttaagcacacaaattattaaaaaagtggtttttgaaagccatgcgtctgacatctatctgattgttgaagaaggattctacaagcgtaccttggatattcaccgaacactagggttgctcttgcacactcaagtctcgattcagcaattgcttaagctgccagcggaatgctttcatcctaaaccaaaagtaaacagtgtcttaataaaacttacccgccataccacagatgttccagataaatattggaagctatatacgtactttgtttcaaaatgggtcaatcgagaatatcgtcaactgtttactaaaaatcagtttcatcaagcaatgaaacacgccaaagtaaacaatttaagtaccgttacttatgagcaagtattgtctatttttaatagttatctattatttaacgggaggaaataattctatgagtcgcttttgtaaatttggaaagttacacgttactaaaggcataaaaataagaagcctgcaaatgcaggcttcttatttttatg"
Len <-nchar(x)+1
Bank_of_Frags<-c()
substring(x, 1, Len)
for (i in 1:Len){
  NewX<-substring(x,1,Len- i )
  Len_NewX<-nchar(NewX)+1
  for (j in 1:Len_NewX){
    NewNewX<-substring(NewX,1+j,Len_NewX)
    if (nchar(NewNewX)>29 & nchar(NewNewX)<31) {
      Bank_of_Frags<-c(Bank_of_Frags,NewNewX)
    }
    
  }
}

print(Bank_of_Frags)

# Time to align

Len_of_Frags<-length(Bank_of_Frags)
Len_of_Frags
Align_Bank<-c()

for ( k in 1:Len_of_Frags){
  fasta_file=readDNAStringSet(filepath = "Original_trimmed.fq", format="fastq", nrec=-1L, skip=0L, seek.first.rec = FALSE, use.names = TRUE)
  fasta_sequence_forward = fasta_file
  DNA<-Bank_of_Frags[k]
  sigma_motif = DNAString(x=DNA, start = 1, nchar = NA)
  result =vmatchPattern(sigma_motif, fasta_sequence_forward, max.mismatch = 1, fixed = FALSE)
  Align_Bank<-c(Align_Bank, result)
}

print("Done")

View(Align_Bank)

# Check and collapse all the data to which some should have aligned.

data_len = length(Align_Bank)
for (l in 1:data_len){
  Store = unlist(lapply(Align_Bank[[l]]@ends,function(x)
    ifelse(is.null(x),0,x)
  ))
  Store = Store[!Store == 0]
}

Store
