 # Washington county data from APP
input_file <- 'dist/WaCoMeF97.txt'
output_directory <- 'dist'
script_name <- 'mccune_w_replacement'

## "sample" is the "little N" -- going ahead and varying from 100 to 1k by 50s
sample_min <- 100
sample_max <- 1000
sample_step <- 50

# fixed sample iteration @ 1k as suggested
iteration_num <- 1000

# Generate output file name
path_split <- unlist(strsplit(input_file, '/'))
path_end <- path_split[length(path_split)]
time = strftime(Sys.time(), format="%Y%m%d-%H%M%S")
output_file <- paste(output_directory, "/", script_name, "_", time,"_it", iteration_num,"_samp", sample_min,"_", sample_max, "_", sample_step,"_", path_end, sep='')

# Read in data
taxa_frame=read.table(input_file)
col = names(taxa_frame)[1]
taxa_list = taxa_frame[,col]
# Find sp count
big_n_max <- length(taxa_list)

# Prepare output collections
little_ns <- c()
big_ns <- c()
statistics <- c()
means <- c()
medians <- c()
sds <- c()

# Vary sample_size (little_n) from 100 to 1k by 50s
for (sample_size in seq(sample_min, sample_max, sample_step)) {
    for (i in 0:50){
        little_n <- sample_size
        big_n <- little_n + 2^i
        if (big_n > big_n_max){
            # Exit if big_n gt big_n_max
            break
        }
        # create array of size big_n sampling W/ replacement
        big_n_array <- sample(taxa_list, big_n, replace=TRUE)

        # prepare collections for storing iteration results
        simpsons <- c()
        shannons <- c()
        unique_taxa_counts <- c()

        # 1k times
        for (i in 1:iteration_num){
            # create array of size little_n sampling from big_n_array W/ replacement
            little_n_array <- sample(big_n_array, little_n, replace=TRUE)

            # Group data for simpson and shannon 
            grouped_taxa <- data.frame(table(sapply(little_n_array, function(m) m)))
            grouped_taxa <- as.data.frame(t(grouped_taxa[,-1]))

            # Append observation of iteration to collections
            unique_taxa_counts <- c(unique_taxa_counts, length(unique(little_n_array)))
            simpsons <- c(simpsons, vegan::diversity(grouped_taxa, index='simpson'))
            shannons <- c(shannons, vegan::diversity(grouped_taxa, index='shannon'))
        }

        # Insert mean, median and sd of observations
        little_ns <- c(little_ns, c(little_n, little_n, little_n))
        big_ns <- c(big_ns, c(big_n, big_n, big_n))
        statistics <- c(statistics, c('simpson', 'shannon', 'unique'))
        means <- c(means, c(mean(simpsons), mean(shannons), mean(unique_taxa_counts)))
        medians <- c(medians, c(median(simpsons), median(shannons), median(unique_taxa_counts)))
        sds <- c(sds, c(sd(simpsons), sd(shannons), sd(unique_taxa_counts)))
    }
}

# Generate output file
diversity_metrics <- data.frame(little_ns, big_ns, statistics, means, medians, sds)
names(diversity_metrics) <- c('little_n', 'big_n', 'statistic', 'mean', 'median', 'sd')
write.table(diversity_metrics, output_file, col.names=TRUE, row.names=FALSE)