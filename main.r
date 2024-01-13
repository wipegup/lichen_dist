## Make changes here
iterations <- 1000
samples <- 200
input_file <- 'dist/WaCoMeF97.txt'
output_directory <- 'dist'
########

path_split <- unlist(strsplit(input_file, '/'))
path_end <- path_split[length(path_split)]
time = strftime(Sys.time(), format="%Y%m%d-%H%M%S")
output_file <- paste(output_directory, "/", time,"_it", iterations, "_samp", samples,"_", path_end, sep='')

taxa_frame=read.table(input_file)
col = names(taxa_frame)[1]
taxa_list = taxa_frame[,col]
taxa_counts <- c()

for (x in 1:iterations) {
    taxa_sampled_with_replacement <- sample(taxa_list, samples, replace=TRUE)
    unique_taxa_count <- length(unique(taxa_sampled_with_replacement))
    taxa_counts <- c(taxa_counts, unique_taxa_count)
}

write.table(taxa_counts, output_file, col.names=FALSE, row.names=FALSE)