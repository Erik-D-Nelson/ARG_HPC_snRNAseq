### MNT add 06Apr2021 =======
  # Check Read 1 files for discrepant read lengths, as seen with 2937 (interactively):
  # Should be exactly 28bp = 16 [BC] + 12 [UMI]

library(ShortRead)
library(jaffelab)

mmFASTQs <- "/dcl01/ajaffe/data/lab/singleCell/mouse_10x/FASTQ/"

R1files <- c(sapply(c("2937", "3182", "3185",
                      "2985", "2986", "2987", "2988"),
                    function(x){list.files(paste0(mmFASTQs,x),
                                           pattern="R1")})
)

for(i in R1files){
  cat(paste0("Checking R1 length distribution for: ", i, "\n"))
  temp.R1s <- readFastq(paste0(mmFASTQs, ss(i,"_",1), "/", i), withIds=F)
  print(head(sread(temp.R1s), n=4))
  print(table(width(sread(temp.R1s))))
  rm(temp.R1s)
  cat("\n\n")
}

# Checking R1 length distribution for: 2985_S1_L001_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 CAACGATTCAGCTTNCCCCTTTAAGAGT
# [2]    28 CTGAGGCGTTACACNGGGCACTAATACG
# [3]    28 TAACTTCCACTGTCNGTTTTAAATTCAA
# [4]    28 AGGTCATTCGCTTTNTTTGTCAGGTATG
# 
# 28 
# 54995033 
# 
# 
# Checking R1 length distribution for: 2985_S1_L002_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 CTGGACGGTCTTGNGGGCGTAGCTCTTG
# [2]    28 CAAGACTAGCATGNTCCGAAAATCACCT
# [3]    28 ATAGACCGTATACNGAATTATGTGCTGA
# [4]    28 TCTGCCATCATCTNTCTAGAGTCCATGT
# 
# 28 
# 53507493 
# 
# 
# Checking R1 length distribution for: 2985_S1_L003_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 CCTTGTGAGCGATGCACCGACGTCCACG
# [2]    28 TGAGGAGCAGACAAGCATCACGTACGAT
# [3]    28 CAAGACTCACAGTATCATCGTTTTAGCG
# [4]    28 CATTGAGTCGGTAGAGGAGCCTCTGCGT
# 
# 28 
# 54547132 
# 
# 
# Checking R1 length distribution for: 2985_S1_L004_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 NGGGTCCAGGTTATAGTAAAGGCCCTGT
# [2]    28 NCATTTCGTGCCCACACTCGGTCCTTTT
# [3]    28 NTCTACCAGAGCCTGATCAGTAAAATCT
# [4]    28 NGGGAGTGTCTCCCTACCGCACTTGGAT
# 
# 28 
# 54136725 
# 
# 
# Checking R1 length distribution for: 2986_S1_L001_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 CCACANTTCATCGCTCTACCGGAGTACG
# [2]    28 TCCCANTTCAGCCTCTCATAAATGTCCG
# [3]    28 CATTGNTGTACCTGTAAAGTTATTGGTT
# [4]    28 AGATGNAAGCCTAGGATTAAGTGCCTTA
# 
# 28 
# 63186062 
# 
# 
# Checking R1 length distribution for: 2986_S1_L002_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 GGGACNAAGAGAAGGTGAACAACGGGAG
# [2]    28 AGTTCNCGTCAGTCCGACACTCCTTGAG
# [3]    28 TTACTNTAGACCCTTGACACCGAGCCCA
# [4]    28 AATTCNTGTGGACTAGGACAAGCAGTCC
# 
# 28 
# 62103788 
# 
# 
# Checking R1 length distribution for: 2986_S1_L003_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 NCTCAGCGTTCCTAAGCTGTGCCCCATC
# [2]    28 NTTACCGTCTCAATCTCGATCATTGCCG
# [3]    28 NAGCGTTGTACCCGCAGAGGCTCGCTCG
# [4]    28 NTGTCTTGTGCATACTCGACAATTTCAA
# 
# 28 
# 62805236 
# 
# 
# Checking R1 length distribution for: 2986_S1_L004_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 NTGGGAGAGGTCACTTTTGGTTGTACGT
# [2]    28 NCCTAATTCAACGAGGAAAAGTCATAGA
# [3]    28 NCACGCGTCATGCGGCCGGTACCGTCAC
# [4]    28 NATTTCCAGGCCCAAAGCCGCGCAAGAA
# 
# 28 
# 62394929 
# 
# 
# Checking R1 length distribution for: 2987_S2_L001_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 TTGGGATTCGAGCCNGCCGAAATTGTTC
# [2]    28 GTTACAGACAATAGNAGGCTACCATGTG
# [3]    28 AGACCATTCGTGTGNTGATAGGGTGACC
# [4]    28 TTGTTGTCAGTCCCNAAACCCCGCGCGT
# 
# 28 
# 68041070 
# 
# 
# Checking R1 length distribution for: 2987_S2_L002_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 AGTACTGTCAGGCNATTACAATCAACCC
# [2]    28 CCCTAACGTGGTANGGGGAGACGTTCGA
# [3]    28 GTGATGTTCCTAANCGTTTGTAACACCT
# [4]    28 CTGGTCTAGCCGTNATTGCTCCGGAGCT
# 
# 28 
# 66246613 
# 
# 
# Checking R1 length distribution for: 2987_S2_L003_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 CATACCCAGGGCAAGGAGCTCGGAATAA
# [2]    28 AACACACAGTATGAACACAGGGTCTTTG
# [3]    28 GTTGCGGAGTCATGCTCATCATTTAGTC
# [4]    28 GACTCTCAGCTGGCCTATGCTGAATCTA
# 
# 28 
# 67535680 
# 
# 
# Checking R1 length distribution for: 2987_S2_L004_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 NGCCACGCACTGCGTGGCGCGGTCGATG
# [2]    28 NAACGATCAGGACGATTCCATTTGCAGT
# [3]    28 NGTCCACAGGGTGAGGACTTCCAGTTAA
# [4]    28 NGGGAGAGTCAGCGTCCTTTATCAAACA
# 
# 28 
# 67040887 
# 
# 
# Checking R1 length distribution for: 2988_S2_L001_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 GAATCNTTCGGTAACTGTAGAAGAGGAG
# [2]    28 ACCCTNACATTAGGAACGGAATCCCGCA
# [3]    28 TGCTCNAGTTCTTAGGGCCATAAAGCTT
# [4]    28 CTAACNTGTCAGACTTTGGAAAAACCGG
# 
# 28 
# 58117351 
# 
# 
# Checking R1 length distribution for: 2988_S2_L002_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 AGAACNTCATCCAACAAAGAAGGAACCG
# [2]    28 ACGGANGGTTAAACAGGTATTAGGCGAA
# [3]    28 GGGTANAGTGATGGCATAACCCCGTATC
# [4]    28 GCCTGNTAGAAGATCTCTTTTCCCTGGG
# 
# 28 
# 57100162 
# 
# 
# Checking R1 length distribution for: 2988_S2_L003_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 NAATTCCTCGTTCGCTCGTTGGTCTCGT
# [2]    28 NGGTTACGTACGGGATGTCGCCACACAA
# [3]    28 NCTATCAGTATCAGCTTCCTACCGCGCA
# [4]    28 NCTCCCACATAGGTTCCCAGCTCAGTCT
# 
# 28 
# 57764333 
# 
# 
# Checking R1 length distribution for: 2988_S2_L004_R1_001.fastq.gz
# DNAStringSet object of length 4:
#   width seq
# [1]    28 NATACCCAGGTTCACTTAAGCCCGAGAA
# [2]    28 NTGTTTGTCCGCACTTTTAATCGCTTGA
# [3]    28 NTCCACATCAGCGGAACTTGGGGTTCTT
# [4]    28 NCGGTCGCATAACTCGTCTTTCAAAACA
# 
# 28 
# 57396714 



sessionInfo()
