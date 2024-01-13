## Make changes here
iteration_min <- 1000
iteration_max <- 10000
iteration_step <- 1000
sample_min <- 100
sample_max <- 500
sample_step <- 50
input_file <- 'dist/WaCoMeF97.txt'
output_directory <- 'dist'
########

path_split <- unlist(strsplit(input_file, '/'))
path_end <- path_split[length(path_split)]
time = strftime(Sys.time(), format="%Y%m%d-%H%M%S")
output_file <- paste(output_directory, "/", time,"_it", iteration_min, "_", iteration_max, "_",iteration_step, "_samp", sample_min,"_", sample_max, "_", sample_step,"_", path_end, sep='')

taxa_frame=read.table(input_file)
col = names(taxa_frame)[1]
taxa_list = taxa_frame[,col]
samples <- c()
iterations <- c()
indexes <- c()
means <- c()
medians <- c()
sds <- c()

for (sample_size in seq(sample_min, sample_max, sample_step)) {
    for (iteration_num in seq(iteration_min, iteration_max, iteration_step)){
        simpsons <- c()
        shannons <- c()
        for (i in 1:iteration_num){
            taxa_sampled_with_replacement <- sample(taxa_list, sample_size, replace=TRUE)
            grouped_taxa <- data.frame(table(sapply(taxa_sampled_with_replacement, function(m) m)))
            grouped_taxa <- as.data.frame(t(grouped_taxa[,-1]))
            simpsons <- c(simpsons, vegan::diversity(grouped_taxa, index='simpson'))
            shannons <- c(shannons, vegan::diversity(grouped_taxa, index='shannon'))
            # unique_taxa_count <- length(unique(taxa_sampled_with_replacement))
            # taxa_counts <- c(taxa_counts, unique_taxa_count)
        }   
            samples <- c(samples, c(sample_size, sample_size))
            iterations <- c(iterations, c(iteration_num, iteration_num))
            indexes <- c(indexes, c('simpson', 'shannon'))
            means <- c(means, c(mean(simpsons), mean(shannons)))
            medians <- c(medians, c(median(simpsons), median(shannons)))
            sds <- c(sds, c(sd(simpsons), sd(shannons)))
    }
}
diversity_metrics <- data.frame(samples, iterations, indexes, means, medians, sds)
names(diversity_metrics) <- c('samples', 'iterations', 'index', 'mean', 'median', 'sd')
write.table(diversity_metrics, output_file, col.names=TRUE, row.names=FALSE)