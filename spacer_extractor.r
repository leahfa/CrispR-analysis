library(tools)			#basic library for filename manipulation
library(Biostrings)		#library for fasta files manipulation

#Hardcode all arrays and their repeat sequences:
FULL_A <- "CCTCTATATTCAGCACACAGTTTGTATCAATTACTGAGGATTTAAACAGAGGCGTTCGTCTGTAACCGTAGACCGAGCTGTTGATGCCGGTCGACCCAGAGCTTCAACCCAACTAGGGTTCGTCTGTAACCTGTTGCTCTTGAATGTAACATGACGACGAAATTCACTGGCGCGCCCGCTTCAACCCAACTAGGGTTCGTCTGTAACGGTGCTGTCGATGGTGATAGACTCGTCCCCGTCCACGCTTCAACCCAACTAGGGTTCGTCTGTAACTCGTCGTCGGGCGTGGGTTCGACTACGATTTCCGCGCTTCAACCCAACTAGGGTTCGTCTGTAACAGATTTTCAAGCGCGAGTTCGGTATTTGACTGCTGCTTCAACCCAACTAGGGTTCGTCTGTAACCCGGGGTTTACGACGTGAAGAAGGTGCGTCCAGAAGCTTCAACCCAACTAGGGTTCGTCTGTAACGCCCTGGTCGTCGACCGGCCCGAAGACTGTCACGGCTTCAACCCAACTAGGGTTCGTCTGTAACCTGCGACCCTCCGGAATTGGGGTCTGCCCCAACGTAGCTTCAACCCAACTAGGGTTCGTCTGTAACAAGCGTTGTAGAACAACATGGCCAGCTTCCAGACGCTTCAACCCAACTAGGGTTCGTCTGTAACTCTCACAACCCCCCCAAATCCTTCCCCCCATGCCGTGTGCTTCAACCCAACTAGGGTTCGTCTGTAACACCTCTGGATAAATCTCCCCAATCTGCCGGAGGACCGGCTTCAACCCAACTAGGGTTCGTCTGTAACGATATGCGTTCCGGTCGGGCAACATTTCGAGACCGTGCTTCAACCCAACTAGGGTTCGTCTGTAACAATTTCACATGGATTTCATGGATTGTCCCCTCGTCAAGCTTCAACCCAACTAGGGTTCGTCTGTAACACGTCATGCATGCAACGCCGATTGTTGCATTCGCGGCTTCAACCCAACTAGGGTTCGTCTGTAACACGCGAAGGCCGGGTCTGCGGAGTACATGGGTCTGGTCGCTTCAACCCAACTAGGGTTCGTCTGTAACGACGCTCGAAGCCGGCGCGCTCGTCCTCGCCGACAAAGGCTTCAACCCAACTAGGGTTCGTCTGTAACTAGATGGTTCATACGTTACTCCCTCAATAATCTTCGCTTCAACCCAACTAGGGTTCGTCTGTAACGTCCAAATCGGGCAGGTGTACGACCTTCGCGTCGAAGCTTCAACCCAACTAGGGTTCGTCTGTAACCGGCGGGCCGTCGGGTCTGCGACACCTCCCACTCCCCGCTTCAACCCAACTAGGGTTCGTCTGTAACCCGACGTTAGATTTTGTGACTCGTTTTGTACCTCCGCTTCAACCCAACTAGGGTTCGTCTGTAACATTTAGCAGAGTCAGTCATTCCGAACCACCGTCCGCGCTTCAACCCAACTAGGGTTCGTCTGTAACTGTCGAAGGCGACTACGTCTTAGCCGTCCATCCGTGGCTTCAACCCAACTAGGGTTCGTCTGTAACACCATCGACGGCACGGTCGTGGGCGACCTCGTACACGCTTCAACCCAACTAGGGTTCGTCTGTAACCGACTGTCGCGTCTTCGATTGGGGTTTTGGGTCCGCTTCAACCCAACTAGGGTTCGTCTGTAACCAAGAGCGTGCGCTCGGAAGTCAACGCAACGCGCAAGGCACACGACCTGACGTGGCTTCAACCCAACTAGGGTTCGTCTGTAACCATGGCCGATTTGAGTGCCACAAGTACTCTACAGGCGGTTCGCACTCAAGTGTTTCCGTCGACCTCCAATACCCCTCAACCCCCCGGGGGGTCGACGGAA"
REVERSE_A <- toString(reverseComplement(DNAString(FULL_A)))
REPEAT_A <- "GCTTCAACCCAACTAGGGTTCGTCTGTAAC"

#data set B: Start, 4 repeats
FULL_B <- "TTCCGTCGACCCCCCGGGGGGTTGGGGGGTATTGGGGGTCGACGGAAACTGTTGAGTGGGAGTAGTGTGTAGGAGGCTGTATACCCTCGAATCGGGCATGGTTACAGACGAACCCTAGTTGGGTTGAAGCGAACAGGATGGCGAACGCGGTGTCTGCGCCAGTTGTTACAGACGAACCCTAGTTGGGTTGAAGCCACGACAATCAAGTCTGGTTGCATGGCGACACGGAGTTACAGACGAACCCTAGTTGGGTTGAAGCCTGTGGCCTCCAGCGGCCGTCAGACAGTCGCATCCGAGTTACAGACGAACCCTAGTTGGGTTGAAGCAAGAAGCCGCTCGCCGTCCTCGATGACGGGCGGGCGGTTACAGACGAACCCTAGTTGGGTTGAAGCGACAAGACTCGCGACGAAGCCGAGTCGAAACGCCGCGTTACAGACGAACCCTAGTTGGGTTGAAGCCTCTTTATCCCTCCTGCCCGAATGTCTACGAATATCGTTACAGACGAACCCTAGTTGGGTTGAAGCGAACCCACTGGTGAAGAAAAAGTTGTAGAGACCCTAGTTACAGACGAACCCTAGTTGGGTTGAAGCACGACAATCAAGTCTGGTTACATGGCGACACGGATGGGTTACAGACGAACCCTAGTTGGGTTGAAGCTTCCACAACGTCGGGGAGGGCGAAATTAGCCAAGCAGTTACAGACGAACCCTAGTTGGGTTGAAGCTCCCGCTGGGGATGTCGGGAGTGCCGGGCGAGCCAGTTACAGACGAACCCTAGTTGGGTTGAAGCCCCGGCCCGTTGCCCCCCACGGCAATCGTCTGCTGTTACAGACGAACCCTAGTTGGGTTGAAGCCGTCTGTGTTATTCTGTGCGTCTGCCGCGACAACGTTACAGACGAACCCTAGTTGGGTTGAAGCATTGCCTGTCCCCGTCGTGTAATCAACTCGGAATCGTTACAGACGAACCCTAGTTGGGTTGAAGCGAGATGTGCGACCGCGGCGAAATGAGCAGTTCGTGGTTACAGACGAACCCTAGTTGGGTTGAAGCGCGACATGGGGACCGTCGAGAACGCGCTCTATGGGGAGTTACAGACGAACCCTAGTTGGGTTGAAGCCGAGGGTCCCGGTGTCGAGAGGACCGGGACGGACGGAGTTACAGACGAACCCTAGTTGGGTTGAAGCTCGGTAATCTGGGAGGCGTCAGTCTCGGCCGAGTAATCGTTACAGACGAACCCTAGTTGGGTTGAAGCCTCGCCATCGCCGCGAACTCGGTCCTCCTCGGGGTGGTTACAGACGAACCCTAGTTGGGTTGAAGCAAGCCTTGAGAGTGTCTGTTGGTATGATGAATGTTGTTACAGACGAACCCTAGTTGGGTTGAAGCAAGTAGACCGCGCTCAGTTACGACAGCTGCTCGAGTTACAGACGAACCCTAGTTGGGTTGAAGCACGATGATCTCGCCAGTCTGCAGCGTTACATTGGGTTACAAACGAATCTTTTCTCGTGAGGACTTCCGAAACTAACCTCTTCCCGGAACTAGTCGAAGGT"
REVERSE_B <- toString(reverseComplement(DNAString(FULL_B)))
REPEAT_B <- "GTTACAGACGAACCCTAGTTGGGTTGAAGC"

#data set C: Start, 4 repeats
FULL_C <- "TCATTTCCGTCGACCCCCCGGGGGGTACAGGGGAATTGAGGGTCGACGGAAACGTTGATGTGAGTTCGCTATGTACAGCCGTCTACCCCCGTGAATTCGGACGGGTTTCAGACGAACCCTTGTGGGTTTGAAGCTCGGGGTGGACGTTCTGCCGGGCATAGGCAATATATTGGTTTCAGACGAACCCTTGTGGGTTTGAAGCGAGCGGAGTCGATGGAAGCCGAGACGGAACAACTGTGTTTCAGACGAACCCTTGTGGGTTTGAAGCTCGGTGGTGTCGCGTATCGAATCACACGTCTGTGGTTTCAGACGAACCCTTGTGGGTTTGAAGCCGCTGGCGTCCGAGTAATTACGCCTGCCCCGAGTGCGGTTTCAGACGAACCCTTGTGGGTTTGAAGCAGCACCTCCGTGGGAGTTGTCGATTCTATCGAAGCGGTTTCAGACGAACCCTTGTGGGTTTGAAGCTGCTAGTCGGACCGAGGGAGGAAGCCTCTCGGATCGTTTCAGACGAACCCTTGTGGGTTTGAAGCGACAACGTTCGCCACCCCGAGACGGGTAACTTCACAGGTTTCAGACGAACCCTTGTGGGTTTGAAGCTGTGGGGAACAAGCAATAGAGGACTTTGAAGGATGGTTTCAGACGAACCCTTGTGGGTTTGAAGCCGTCCGGCGACCCGAATATGTCCGTCCCCGTGCGGGTTTCAGACGAACCCTTGTGGGTTTGAAGCGGCGACAACGACGGCAACCTCGAAGTCTGCGCCGCCTGTTTCAGACGAACCCTTGTGGGGTTGAAGCGGCACCTACGACCGCCGCTACAACGACGACCCGCGTTTCAGGCGAACCTGTGTGGGATTAAATC"
REVERSE_C <- toString(reverseComplement(DNAString(FULL_C)))
REPEAT_C <- "GGTTTCAGACGAACCCTTGTGGGTTTGAAGC"

#data set D: Start, 4 repeats
FULL_D <- "TTCGTTTCCGTCGACCCTCGGGGGGTACGGGGGAATTGAGGGTCGACGGAAACGTTGATGTGAGTTCGCTATGTACAGCCGTCTACCCCCGTGAATTCGGACGGGTTTCAGACGAACCCTTGTGGGATTGAAGCGCAGGCATCTCGACCGGCGACCTCCCGGAACACTTTGGTTTCAGACGAACCCTTGTGGGATTGAAGCTGCACAAGTCCATCAAAGAAGACGGTTACACGCAACCGTTTCAGACGAACCCTTGTGGGATTGAAGCGTCCTCGACGAAGCCGGGGTCCAATTCGGGAACCGGTTTCAGACGAACCCTTGTGGGATTGAAGCGCGACCGGGCCGAACGAAACCGCGCAGAACATCACGGGTTTCAGACGAACCCTTGTGGGATTGAAGCAGGACACGTCTGCGCTGCGTGCGCCTACTGCATGAGGTTTCAGACGAACCCTTGTGGGATTGAAGCGCGCACGCGTAGGCCCCGAGCACGAACGCGACCGGGGTTTCAGACGAACCCTTGTGGGATTGAAGCTGATTCTCGAAGGAGCCGAGCGCGGCTTCCGTGGAATGTTTCAGACGAACCCTTGTGGGATTGAAGCGTGAAGGTCAAACGCTTCTGGTAGGCAGGCCCGCGGTTTCAGACGAACCCTTGTGGGATTGAAGCAGTCTGCGGTCCCGGCGTTGTGAATGTGAACGACTTGTTTCAGACGAACCCTTGTGGGATTGAAGCGCCGATTGCGATATTTTTCAGAATCGGCGACTGTTGTTTCAGACGAACCCTTGTGGGATTGAAGCGGTGTCCAAATGTCGCCAGAGAATTCCTGTCAAGCTTGTTTCAGACGAACCCTTGTGGGATTGAAGCTGACGCTGTTCGATGTCTTTCAACTCGTTTGTCTCTCGTTTCAGACGAACCCTTGTGGGATTGAAGCGGCCTCCGCGGCGGGCGGTAAAGCCGATTCATCCTGTTTCAGACGAACCCTTGTGGGATTGAAGCGACGGTTCCGTTCTCGTCGTAGATCGTCACGTCCAGTTTCAGACGAACCCTTGTGGGATTGAAGCGCGCGACAGAGCCCGTCGAGATCTCCGAAGGCGAGAGTTTCAGACGAACCCTTGTGGGATTGAAGCTCTCGGCAACAGACCCTGCGTTAGTAATTTGGACACTGTTTCAGACGAACCCTTGTGGGATTGAAGCAGATTGTGGATGAACGTCACGCCAAGCCGACCGGGAGTTTCAGACGAACCCTTGTGGGATTGAAGCGACGAGACGATCGGGCTTCTGGAAGACCAGGAGGAAGGTTTCAGACGAACCCTTGTGGGATTGAAGCTTCAGCGACTACGACCTGATGCGCTCCGATACAGGTTTCAGACGAACCCTTGTGGGATTGAAGCATGATCCGCGTGGTCCGGTTGACGTTCGGCTGCTCTTGTTTCAGACGAACCCTTGTGGGATTGAAGCTATCACCCGTGGGCCGGTGACAACGAGGACCACGAGTTTCAGACGAACCCTTGTGGGATTGAAGCATTTTGCCGGACGCTTCGGGAGTTGATTCATAGTTGTTTCAGACGAACCCTTGTGGGATTGAAGCGTCGCTCGCGGGATCGACATCATCGCGCTCAAGGCGAGTTTCAGACGAACCCTTGTGGGATTGAAGCTAGTCCGCCCCGCACTCAGGGCACGCGTCGTCGTCCGTTTCAGACGAACCCTTGTGGGATTGAAGCTGCGCGTGCAGCTCGTTTGATGAGACGTCACCACCGGTTTCAGACGAACCCTTGTGGGATTGAAGCGAGGGCGGGGGCAACGTCGCGGAGATGCAGGACCTCGTTTCAGACGAACCCTTGTGGGATTGAAGCGAGACGATTGAGTGGGGCGGTGACTCCGCTATGGCAGGTTTCAGACGAACCCTTGTGGGATTGAAGCTATTACAATTGCCAAAAACCACGTTGCTAATAGCTCCGGTTTCAGACGAACCCTTGTGGGATTGAAGCGGGAGTGTATGCGATATCCCTTAAAACCCCCTACTTGTTTCAGACGAACCCTTGTGGGATTGAAGCACTCGTCGCACGTGTTAGTTTTCGGTTCGCTCGTTTGTTTCAGACGAACCCTTGTGGGATTGAAGCGGCGGGCTCGGCAGGAGCTCGAAGACCGCTTTTCGGGTTTCAGACGAACCCTTGTGGGATTGAAGCACCGACGATGTAACTGTCTCGGCGGCAACGCTCCCGTTTCAGACGAACCCTTGTGGGATTGAAGCACTCTCCGCGTCCGACTGCTCGACATCCCGACATCGGTTTCAGACGAACCCTTGTGGGATTGAAGCACCGCGAAATCGGCGTCAATCGTCACGTCCGGTAGGGTTTCAGACGAACCCTTGTGGGATTGAAGCCACTCTGCGAAGCCGTCGATGAAGTCCTCGAACTCCGTTTCAGACGAACCCTTGTGGGATTGAAGCTTTTCTACGAGGTCGGAAACGGACACTTCGACGTTGTTTCAGACGAACCCTTGTGGGATTGAAGCTCGAAGCCGCAGACGACGGCGCAGAGAGCGGTATCACGTTTCAGACGAACCCTTGTGGGATTGAAGCCGGCTCGACGACGATTACTTCGCGCTCGAAGCTCGTGTTTCAGACGAACCCTTGTGGGATTGAAGCGTCGCTCGCGGGATCGACATCATCGCGCTCAAGGCGATGGTTTCAGACGAACCCTTGTGGGATTGAAGT"
REVERSE_D <- toString(reverseComplement(DNAString(FULL_D)))
REPEAT_D <- "GTTTCAGACGAACCCTTGTGGGATTGAAGC"

#data set E: End, 4 repeats
FULL_E <- "CTGGGACATGCGGAGATCGTGGGTTCGGAGCCCACCCCGGGCATTCTTCGACTTCAGTACGTGATTGACTAATACGGGCATCGCAATGCCGTCTCCTGATTCTTCACTTCGATACGGGTTCGTCTGAAACATCCACCGTGCCGCCGAACGAGACCACGTTAGTCGAGCTTCAACCCCACAAGGGTTCGTCTGAAACGCGGGCAACCACCCGAGCGGAGGTGTCGACGCATGAGCGCTTCAACCCCACAAGGGTTCGTCTGAAACACTGCATGAGCGCCTGGGCAGTCCCGTCGACGGAGCTTCAACCCCACAAGGGTTCGTCTGAAACGGTCGGACCGCCACCGAGCCAGCGCTCGCGACGGCCGCTTCAACCCCACAAGGGTTCGTCTGAAACCTTTTCCTCGCGCCCCTCGCAAAGCACTCGATAGGCGGCTTCAACCCCACAAGGGTTCGTCTGAAACTCGGACGCCGACCGCCTGGACATCTACGTCTACGATCCGGCTTCAACCCCACAAGGGTTCGTCTGAAACAGTCTTCGCCGTCGCCCACCGCCATGGCGTTATCCAGCTTCAACCCCACAAGGGTTCGTCTGAAACGCGTGTCGTTTGAGGTGTCGACAACATTTGCAAAATTGCTTCAACCCCACAAGGGTTCGTCTGAAACCCGCTGCCCCCGAGTTGTCAACCGCTTCAAGTGCGGTGCTTCAACCCCACAAGGGTTCGTCTGAAACCGGTGTCGAGGTCGTTCGCTCCGGCGGGGATTTCGACGCTTCAACCCCACAAGGGTTCGTCTGAAACCGTCCGATGCGCTGTCGGCGTCGATGCGGAACTCGATGCTTCAACCCCACAAGGGTTCGTCTGAAACCACCGTGTTCAGTCTCTTGTCGCAGTCGGATGTCGCTTCAACCCCACAAGGGTTCGTCTGAAACTCGCGTCGCCATCACGGCGCTGACATTCGCCTACGCGCTTCAACCCCACAAGGGTTCGTCTGAAACGGTGTTCGAGGAATACATGGCCCCGGAGGGCGAGGACGCTTCAACCCCACAAGGGTTCGTCTGAAACAATCTTCCGTGGGATCTCGTAGTAATCGAGGATTTTGCTTCAACCCCACAAGGGTTCGTCTGAAACCGCTCGATACCGTCGTTGCTGTCAGCGACCTCCGCGCTTCAACCCCACAAGGGTTCGTCTGAAACACGATTTACGAGATTATCGAAAGAATGGACGTTGAAAGCTTCAACCCCACAAGGGTTCGTCTGAAACACCAACGACTCGGATGGAAGCTGGTGAGCGATATGCTTCAACCCCACAAGGGTTCGTCTGAAACTGGGATTCGGTCTCCACAAAGACCCCACTAACTGTGCTTCAACCCCACAAGGGTTCGTCTGAAACGGGTTCGGCATGCTCCCATAGCCGCCAGTGCCGTCGCTTCAACCCCACAAGGGTTCGTCTGAAACGACGTCGACGATGTTGTCTCGGACGAGGTCGGGCATGCTTCAACCCCACAAGGGTTCGTCTGAAACTGGGGGCGTGGAACCAATGAACCGGTCGTCGACCGGCTTCAACCCCACAAGGGTTCGTCTGAAACCGTGACTTCGAAGTCGCCGCAGACGTGCGTTCCCTCGCTTCAACCCCACAAGGGTTCGTCTGAAACTTTCTGAGATTCGAGGGCATCTTCGGACCTTTCCTCGCTTCAACCCCACAAGGGTTCGTCTGAAACCATCCCGAATTCACGGGTACAATGATGTGTATCACCCGAATTCAAAAGAGTGTTTCCGTCGACCCTCAATTCACCCGCAACCCCCGGGGGGTCGACGAAA"
REVERSE_E <- toString(reverseComplement(DNAString(FULL_E)))
REPEAT_E <- "GCTTCAACCCCACAAGGGTTCGTCTGAAAC"

#data set F: Start, 4 repeats
FULL_F <- "AGCAGTCGACCCCCCGGGGGGTAGAGGGGTATTGGAGGTCGACGGAAACACTTGAGTGGGAACCGTCCGTAGAGCACTTGTGGCGCTTAAATCGGCCATGGTTACAGACGAACCCTAGTTGGGTTGAAGCATCCGGTCCCAGTAGAGCGCGTAGACGCCGAGGGCGTTACAGACGAACCCTAGTTGGGTTGAAGCCTTTGTCGGCGAGGACGAGCGCGCCGGCTTCGAGCGTGTTACAGACGAACCCTAGTTGGGTTGAAGCGAGCTTGTCGAAGACGATTGGGTCGCTCACACTCCGTTACAGACGAACCCTAGTTGGGTTGAAGCCGACAATATTCGACCAATTGGACATATGCATCCTGTTACAGACGAACCCTAGTTGGGTTGAAGCGCCCGTTTAACTAATAATGTCATGGTTCGGAAAGCGGGTTACAGACGAACCCTAGTTGGGTTGAAGCTCGACACGACAACACGACTCGATTCGCCGTACGACGGTTACAGACGAACCCTAGTTGGGTTGAAGCTACACCGGACCCGAACTCGACGTGCACGTGAAGCCCGTTACAGACGAACCCTAGTTGGGTTGAAGCTCTACCAACTGGTTGGGTGCTGCTCGGATATCTCTCGTTACAGACGAACCCTAGTTGGGTTGAAGCGACCGGGACCCGTGATGCTCGACAAGGTTGTAGTGTTACAGACGAACCCTAGTTGGGTTGAAGCACGCTCGGTCCCGGCGAGACGGCGACGTTCGCGCCGGTTACAGACGAACCCTAGTTGGGTTGAAGCTAAACGCATCTAATGAGTTATCAGACTCGTGGACTGGTTACAGACGAACCCTAGTTGGGTTGAAGCACCTCATAGCGCGGGTCCGCCCGTCTGGGTTGGCGTCTGGTTACAGACGAACCCTAGTTGGGTTGAAGCCTTGACAACGGCTTTGACGGTCGGACAATGACCTATCCGTTACAGACGAACCCTAGTTGGGTTGAAGCTGTCGTGACTGCGCCCACGTCATCACTATCGGTGATTGTTACAGACGAACCCTAGTTGGGTTGAAGCCTCAAGAGGTGGCTCTAATGGGAATCAAGCAGTTGTTACAGACGAACCCTAGTTGGGTTGAAGCAGCGTCGACGACGTGCTCGACCGACTGGAGTCTGACCCGTTACAGACGAACCCTAGTTGGGTTGAAGCCGCGCACGGCATCCACGTCGACGAGCGGGACGCCCGGTTACAGACGAACCCTAGTTGGGTTGAAGCACGTCCTCGGTAGGATGTTTCGCCGCGAGACTCTGCAAGTTACAGACGAACCCTAGTTGGGTTGAAGCGGGAACCGTTGAGGTCCTCGTCGATGATGCGAGCGTTACAGACGAACCCTAGTTGGGTTGAAGCATCGACCGGAACACGCGAACCACGCTCGGCGAGCGTTACAGACGAACCCTAGTTGGGTTGAAGCACTCCAACGAGTTCATTAATCCCGATGGGGAGGAGGTTACAGACGAACCCTAGTTGGGTTGAAGCTTACAGACGATTCGGTCGACAGAACCCTATCGCCCGTTACAGACGCACCCTTGCGGGATTGAATCCCGATTCCAAATGTGACAGGAGGAAATTCAGACCCATCGGTAAACTTTAAATTCCGTAGTCTGAATACTCTTCCAGTACATCTATCGTATAGAGTTTCGT"
REVERSE_F <- toString(reverseComplement(DNAString(FULL_F)))
REPEAT_F <- "GTTACAGACGAACCCTAGTTGGGTTGAAGC"

#data set G: End, 4 repeats
FULL_G <- "CTACTGGGTATTTCCTAGTGGGGGCGCTCTCGCTGCTGTTCCTGTCGATAGTGGTTCTTTTGTAACGACGACTAACGTCACTGTCCGTGACCTTGGTCTCGCTTCAACCCAACTAGGGTTCGTCTGTAACCCGGCCGTGATTGGCACGACGTGCCACGTGTTCGCGACGCTTCAACCCAACTAGGGTTCGTCTGTAACCTGAGCGGCGGTCTCAAGAGGCGACTCGACGATGTCGCTTCAACCCAATTAGGGTTCGTCTGTAACTGCATGGACGATCTCGTTCACGTCCTCCGGCAGTCGGCTTCAACCCAATTAGGGTTCGTCTGTAACCAGGCGTTCTTGCGCGCCTGCCACCTTCTCGGCGGGATGCTTCAACCCAATTAGGGTTCGTCTGTAACGATATTCGTGGACGTTCGGAGCGGCGGAATAAAGAAGCTTCAACCCAATTAGGGTTCGTCTGTAACACTCCAACGAGTTCATTAATCCCGATGGGGAGCAGGCTTCAACCCAATTAGGGTTCGTCTGTAACGGGGGTGTTCCAAAATGAGCTGGCACCCGGTTCGGGCTTCAACCCAATTAGGGTTCGTCTGTAACACCGACGTAACACTCTCGGTTGATGGGGGAGACCTGCTTCAACCCAATTAGGGTTCGTCTGTAACATTCTCAGGGGGCGTGCCGAGGCGTCCCTCGGTTCCGGCTTCAACCCAATTAGGGTTCGTCTGTAACAGTTGAGACCGTCGCGCCCGGATTCCCATGTGACTGCTTCAACCCAATTAGGGTTCGTCTGTAACCATGCCCAATTCAGGGGCATACAGCCCCCTACACACTGCTCCCACTCAACAGTTTCCGTCGACCTCCAATAACCCCCGAACCCCCGGGGGGTCGACGGAA"
REVERSE_G <- toString(reverseComplement(DNAString(FULL_G)))
REPEAT_G <- "GCTTCAACCCAATTAGGGTTCGTCTGTAAC"

#data set H: Start, 4 repeats
FULL_H <- "TCCGTCGACCCCCGGGGGGGGTTGAGGGTTATTGGGGGTCGACGGAAACTGTTGAGTGTGAGCAGTGTGTAGGAGACTGTATGCCCTCAGATTGGGCATGGTTACAGACGAACCCTAGTTGGGTTGAAGCTCCCGAGGAGTCGTCCGCGTCGACGCGGATGAGTTGGTTACAGACGAACCCTAGTTGGGTTGAAGCGTGCATCGCCGAGCGGTCGTTCGGGTCCATCTTGTGTTACAGACGAACCCTAGTTGGGTTGAAGCTGATGATGGGTCCATCACCACCCACCCCCTGCGTGTTACAGACGAACCCTAGTTGGGTTGAAGCAAGAACCATTGACTGGAACGCAACGACGCATCCTCGTTACAGACGAACCCTAGTTGGGTTGAAGCTGGGGATGTCTCGGACAATCTCTGACCGTGTGATAAAGTTACAGACGAACCCTAGTTGGGTTGAAGCTGACTGGAGAGCGTGCAACGGCCTATAAGATAGACCTGTTACAGACGAACCCTAGTTGGGTTGAAGCTCAACTGTCCCCTTGTCCGCTTCGTCGCCTGAGCGGTTACAGACGAACCCTAGTTGGGTTGAAGCGTACTTCTCTCCATCTTCGAGATTAGTGACGTTCAGGTTACAGACGAACCCTAGTTGGGTTGAAGCGTACTAAAGTGTGGCACGGAAAACACGAACAGCCCTGTTACAGACGAACCCTAGTTGGGTTGAAGCCCATCGTGGGTTGATTCGATACTCGACCGCAACGGTTACAGACGAACCCTAGTTGGGTTGAAGCTCGACTTGCGAGATAATGTCGTCTTGCTCGGTGGTGTTACAGACGAACCCTAGTTGGGTTGAAGCAAATATGATACCACCGAACCCGTTTGTAGCGCCCAGTTACAGACGAACCCTAGTTGGGTTGAAGCGTTGGTGCCCGCGCGGTTGGCGTAGGCAAGGTTCGTGGTTACAGACGAACCCTAGTTGGGTTGAAGCTGAGGTCTTCCCTAACTTCTGCATTGGTGTGTCCCGGTTACAGACGAACCCTAGTTGGGTTGAAGCAAACTCCTCGGCGTCGACACCGAGCGCACGATCGACGGTTACAGACGAACCCTAGTTGGGTTGAAGCCGTTTCCTCATCGAGGACTATTAAATCTCCTGACGAATGTTACAGACGAACCCTAGTTGGGTTGAAGCTTGCATCTCTTCCAGACTCTCCTCTATGCTCCGGTTACAGACGAACCCTTTCTCGTGAGGACTTCCGAAAATAGTCGATTCTCGGAGGTACCCTGACGGGCGACATCAGTGAATCCTTCGTGGACCGATGGACATCATCACCCAAGTAACGAGGATAATCACA"
REVERSE_H <- toString(reverseComplement(DNAString(FULL_H)))
REPEAT_H <- "GTTACAGACGAACCCTAGTTGGGTTGAAGC"

#data set I: Start, 2 repeat
FULL_I <- "TTCCGTCGACCCCCCGGGGGTTTGAGGGGTATTAGGGGGTCGACGGAAACACTTGAGTGAGATGCGACAGTACTAGCTCTATAGCGGCCAGATTGGCCATGGTTACAGACGGACCCTCGTTGGGTTGAAGCGGACGATGCGCAGCGAGTAATCGCCGCCTTCCCGACGGTTACAGACGGACCCTAGTTGGGTTGAAGTGGGGTGAGAGGGGGGTTGTTCGCTTCTCAAACTTCCTTAACTAGACAGTGACACTGACTAGTAGAGTCTGCTAAACCCACAGATAATTGGGAAATGTGA"
REVERSE_I <- toString(reverseComplement(DNAString(FULL_I)))
REPEAT_I <- "GGTTACAGACGGACCCTCGTTGGGTTGAAG"
#create a single variable for the combination of all arrays, in forward and reverse, separated by spaces
FULL_AthruI <- paste(FULL_A, FULL_B, FULL_C, FULL_D, FULL_E, FULL_F, FULL_G, FULL_H, FULL_I)
REVERSE_AthruI <- paste(REVERSE_A, REVERSE_B, REVERSE_C, REVERSE_D, REVERSE_E, REVERSE_F, REVERSE_G, REVERSE_H, REVERSE_I)

#Define the path to dircetory containg the outputs from "process.sh": 
#1. otus_0.99.fna (Each sequence is the "seed" sequence for a single OTU ) , 
#2. tabfile.txt  - per-sample abundance information for each OTU

path="/scratch200/leah/Crisp2018/A/"

## MANUALLY CHANGE HERE THE REPEAT SEQUENCE TO FIT THE RIGHT ARRAY! ##
REPEAT_SEQUENCE <- REPEAT_A

REPEAT_NUMBER<-if (REPEAT_SEQUENCE=="REPEAT_I") {REPEAT_NUMBER<-2} else {REPEAT_NUMBER<-4}
leaderChoice <-if (REPEAT_SEQUENCE==REPEAT_A|REPEAT_SEQUENCE==REPEAT_E|REPEAT_SEQUENCE==REPEAT_G) {leaderChoice<-"NO"} else {leaderChoice<-"YES"}
filename=paste0(path,"otus_0.99.fna")


OTU <- readDNAStringSet(filename)		#assign to variable OTU as a bioString type
all=list()
 
num <- REPEAT_NUMBER		##how many repeats we want to find, depending on array

print(path)
print(num)
print(filename)
print(REPEAT_SEQUENCE)

if (leaderChoice == "YES") {					#allow user to select leader position
	leaderChoice <- "Start"					#rename to make sense
	pos <- 0 } else if (leaderChoice == "NO") {	#result of selection stored in variable "pos", for "position".
	leaderChoice <- "End"					#...
	pos <- (REPEAT_NUMBER - 2) } else {			#used in reference to indexed repeat locations to choose correct spacer
	stop("ERROR: Must select leader position")	#if neither is selected, throw error and terminate program
}
print(paste0("POS is ",pos))
REPEAT_ERRORS <- 2	#LOW number: strong filter
SPACER_ERRORS <- 5	#HIGH number: strong filter

ID=0.99

ori=c() #store names of otus with the original number of repeats

for(m in 1:3) {		#assume fewer than 5 acquisitions.  Loop should terminate on own, but finite loop is a backup.
	filterX <- paste0("filter", num)	#create first gen of each vector
	acqX <- paste0("acq", num)		#as yet uninitialized
	spacersX <- paste0("spacers", num)	#each new generation is newly named
	finalX <- paste0("final", num)	#for multiple acquisitions
	
	assign(filterX, character(length = length(OTU)))	#initialize empty list for filter: max if every cluster has acquisition
	filterTemp <- get(filterX)					#copy filter for temporary use in loop
  counter=1 #to store results in Ori
	for (i in 1:length(OTU)) {				
     if (m==1 & vcountPattern(REPEAT_SEQUENCE, OTU[i],		#check for repeat in each array entry
                    max.mismatch = REPEAT_ERRORS, min.mismatch = 0,		#mutation allowance
                    with.indels = FALSE, fixed = TRUE, algorithm = 'auto'		#default settings
  ) ==REPEAT_NUMBER-1) {								#if in OTU[i], count at least (num) repeats
    ori[counter]=names(OTU)[i]
    counter=counter+1
  }
		if (vcountPattern(REPEAT_SEQUENCE, OTU[i],		#check for repeat in each array entry
		    max.mismatch = REPEAT_ERRORS, min.mismatch = 0,		#mutation allowance
		    with.indels = FALSE, fixed = TRUE, algorithm = 'auto'		#default settings
		   			) >= num) {								#if in OTU[i], count at least (num) repeats
			filterTemp[i]<- paste0(">otu_", i, "\n", toString(OTU[i]))		#then add OTU number, convert to char, add to list
		}
	  
	}
  #end of 1st i loop#
	if (identical(get(filterX), filterTemp)) { 
		print(paste("Finished with a maximum of", m-1, "acquisitions."))
		break
	}								#when num exceeds max number repeats, temp will be empty: break loop
  
	assign(filterX, filterTemp)				#otherwise, copy temp filter back to current filter,
	assign(acqX, filterTemp[filterTemp != ""])	#collapse empty elements, leaving only acquisitions (for reference)

	finalTemp <- spacersTemp <- character(length = length(get(acqX)))	#create new temps of length acqX for next refinement
	assign(finalX, character(length = length(get(acqX))))			#also initialize empty finalX for identical check

	for (j in 1:length(get(acqX))) {	#for each acquisition event, extract only the new spacer(s) from the cluster						
		repeats <- vmatchPattern(REPEAT_SEQUENCE,	get(acqX)[j],		#same parameters, but now find spacer locations
				max.mismatch = REPEAT_ERRORS, min.mismatch = 0,
				with.indels = FALSE, fixed = TRUE, algorithm = 'auto')
		spacer <- substr(get(acqX)[j], endIndex(repeats)[[1]][m+pos] + 1, startIndex(repeats)[[1]][m+pos+1] - 1)	#first or final spacer, depending on pos

		if (nchar(spacer) >= 45 || nchar(spacer) <= 10 ||	#ignore long spacer original spacer, or short (from repeat misread)
		    vcountPattern(spacer, FULL_AthruI,			#compare current spacer to all full sequences, to check for duplications
		    max.mismatch = SPACER_ERRORS, min.mismatch = 0,	#again, same parameters
		    with.indels = TRUE, fixed = TRUE, algorithm = 'auto') ||
		    vcountPattern(spacer, REVERSE_AthruI,			#ALSO compare the reverse complement
		    max.mismatch = SPACER_ERRORS, min.mismatch = 0,
		    with.indels = TRUE, fixed = TRUE, algorithm = 'auto'))
		{ } else {	#if there is a match, do nothing.  If not, add to appropriate temp vectors
			spacersTemp[j] <- paste0(substr(get(acqX)[j], 1, endIndex(repeats)[[1]][m+pos]),
						" ->", spacer, "<- ",									#saves full OTU, marks spacer with arrows
						substr(get(acqX)[j], startIndex(repeats)[[1]][m+pos+1], nchar(get(acqX)[j])))
			finalTemp[j] <- paste0(substr(get(acqX)[j], 1, regexpr("\n", get(acqX)[j])), spacer)	#saves only the OTU number and the spacer itself
		}	
	}
	###end of 2nd loop ####
	if (identical(get(finalX), finalTemp)) { 			#if all spacers removed, end program
		print(paste("Finished with a maximum of", m-1, "interesting acquisitions."))
		break
	}		
	assign(spacersX, spacersTemp[spacersTemp != ""])	#collapse repeated spacers, leaving only acquired, not duplicated, sequences
	assign(finalX, finalTemp[finalTemp != ""])		#final output, including only spacer and number
  get(finalX)->z
  write(z, paste0(path,ID,"_spacers_",m,".txt"),sep="\t")
  all[[m]]=z
  
	
	num <- REPEAT_NUMBER + m	#iterate for next loop
}

write(ori,paste0(path,"ori_",ID,".txt"),sep="\t")
print(num)
all.all=unlist(all)

write(all.all,paste0(path,"all_spacers_",ID,".txt"))

