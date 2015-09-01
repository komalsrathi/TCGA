################################ download rsem.genes.normalized_results ###################################

# move all normalized counts to genecounts folder
system('mv *.rsem.genes.normalized_results ./genecounts/')
# make a file to rename the existing namings 
system('cut -f5,7 file_manifest.txt | grep rsem.genes.normalized_results > file_rename.txt')
# rename count files according to sample names
ren <- read.table('file_rename.txt')
file.rename(from=list.files(pattern=".rsem.genes.normalized_results"),to=list(as.character(ren$V1))[[1]])

################################ download rsem.genes.normalized_results ###################################
