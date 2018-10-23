library(Claddis)
library(mulTree)
source("MorphDistMatrix.verbose.R")
source("anc.state.R")
source("anc.state_fun.R")
source("anc.unc.R")
## matrix
Nexus_data <- ReadMorphNexus("2014-Beck-ProcB-matrix-morpho.nex")
Nexus_matrix <- Nexus_data$matrix
## tree
tree <- read.nexus("2014-Beck-ProcB-TEM.tre")

#######################################
## Cleaning the matrices and the trees
#######################################

## Removing some species
tree_tmp <- extract.clade(tree, 133)
tree_tmp <- drop.tip(tree_tmp, extract.clade(tree_tmp, 127)$tip.label)
tree_data <- drop.tip(tree_tmp, c("Erinaceus", "Ptilocercus", "Orycteropus", "Microgale"))

## Cleaning the tree and the table
tree_data <- clean.tree(tree_data, Nexus_matrix)
table <- clean.table(Nexus_matrix, tree_data)
Nexus_data$matrix <- table
tree_data <- lapply.root(tree_data)

## Removing some characters
Nexus_data$matrix <- table
to_remove <- which(apply(Nexus_data$matrix, 2, function(X) length(levels(as.factor(X)))) < 2)
to_remove <- c(to_remove, which(apply(Nexus_data$matrix, 2, function(X) length(which(is.na(X)))) > 25))
matrix_tmp <- Nexus_data$matrix[,-c(to_remove)]
matrix_tmp <- matrix_tmp[,-c(28,59,75,98,135,172,182,186,189,191,204,224,235)]
matrix_tmp <- matrix_tmp[,-c(57,239)]
Nexus_matrix <- make.nexus(matrix_tmp)
Nexus_matrix$min.vals <- as.numeric(Nexus_matrix$min.vals)
Nexus_matrix$max.vals <- as.numeric(Nexus_matrix$max.vals)

#######################################
## Estimating ancestral states
#######################################

anc_states <- anc.state(tree_data, Nexus_matrix, method='ML-claddis', verbose=TRUE)

#####################################
## Distance matrix
#####################################

## Distance matrix using tips only
matrix_tips <- Nexus_matrix
dist_tips <- MorphDistMatrix.verbose(matrix_tips, verbose=TRUE)

## Distance matrix using also nodes
matrix_nodes <- Nexus_matrix
matrix_nodes$matrix <- anc.unc(anc_states, 0.95, missing=NA)$state
dist_nodes <- MorphDistMatrix.verbose(matrix_nodes, verbose=TRUE)

#####################################
## Creating the two pco's (with and without nodes)
#####################################

ord_data_tips <- cmdscale(dist_tips$gower.dist.matrix, k=nrow(dist_tips$gower.dist.matrix) - 2, add=T)$points
ord_data_tipsNnodes <- cmdscale(dist_nodes$gower.dist.matrix, k=nrow(dist_nodes$gower.dist.matrix) - 2, add=T)$points

#####################################
## Loading the FADLAD data
#####################################

FADLAD_data <- read.csv("Beck2014_FADLAD.csv", row.names=1)

## removing taxa not present in the tree
FADLAD_data <- FADLAD_data[-which(is.na(match(rownames(FADLAD_data), tree_data$tip.label))),]

test_data <- list("tree_data"=tree_data, "ord_data_tips"=ord_data_tips, "ord_data_tips_nodes"=ord_data_tipsNnodes, "FADLAD_data"=FADLAD_data)

save(test_data, file="../test_data.Rda")


#####################################
## Ecology data
#####################################

## Loading the data
McClean_data <- read.csv("2015-McClean.csv")

## Generating the matrix
McClean_matrix <- with(McClean_data,tapply(McClean_data$Abundance,list(McClean_data$Tdepth,McClean_data$Genus),mean))
McClean_matrix <- as.matrix(cbind(McClean_matrix[,1:93]))

## Calculating the distance matrix
McClean_distance <- as.matrix(dist(McClean_matrix, method="euclidean"))

## Ordinating the distance matrix
McClean_ordination <- cmdscale(McClean_distance, eig = TRUE, k = nrow(McClean_distance)-1)$points

## Factors
treatment <- c("a","a","b","b","a","a","a","b","b","a","a","b","b","b","b","b","b","b","a","b","b","a","b","b","b","b","a","a","a","b","b","a","a","a","a","a","a","a","a","a")
depth <- c(1,2,1,2,1,1,2,1,2,1,2,1,2,1,2,1,2,2,2,1,2,1,1,2,1,2,1,1,2,1,2,1,2,1,2,1,1,1,1,1)

## Output list
McClean_data <- list("ordination"=McClean_ordination, "treatment"=treatment, "depth"=depth)



#####################################
## Model test data
#####################################

## Test data
data(BeckLee_mat99) ; data(BeckLee_ages) ; data(BeckLee_tree)
data_bootstrapped <- boot.matrix(chrono.subsets(BeckLee_mat99, BeckLee_tree, method = "continuous", rev(seq(from = 0, to = 120, by = 5)), model = "proximity"))
model_test_data <- dispRity(data_bootstrapped, c(sum, variances))

save(model_test_data, file="../model_test_data.Rda")


#####################################
## LDA test data
#####################################
library(ape)
library(geomorph)
set.seed(42)

## Loading the plethodon dataset
data(plethodon)

## Performing the GPA
procrustes <- geomorph::gpagen(plethodon$land, print.progress = FALSE)

## Sorting the species
Jord_sp <- which(plethodon$species == "Jord")
Teyah_sp <- which(plethodon$species == "Teyah")

## Naming the species in the procrustes object
Jord_sp_names <- paste0("Jord_", Jord_sp)
Teyah_sp_names <- paste0("Teyah_", Teyah_sp)
dimnames(procrustes$coords)[[3]] <- rep(NA, length(c(Jord_sp, Teyah_sp)))
dimnames(procrustes$coords)[[3]][Jord_sp] <- Jord_sp_names
dimnames(procrustes$coords)[[3]][Teyah_sp] <- Teyah_sp_names

## Procrustes matrix
procrustes_matrix <- geomorph::two.d.array(procrustes$coords)

## Species classifier
species_class <- unlist(lapply(strsplit(rownames(procrustes_matrix), split = "_"),
                               function(x) return(x[1])))

## Morpho_groups classifier
morpho_class <- c(rep("group1", 10), rep("group2", 10), rep("group3", 20))

## Random classifier
random_class <- sample(c("random1", "random2", "random3"), 40, replace = TRUE) 

## test lda objects
geomorph_df <- geomorph.data.frame(procrustes, species = as.factor(lda_test_data$species), morpho = as.factor(lda_test_data$morpho))
data <- geomorph.ordination(geomorph_df, ordinate = FALSE)
lda_test <- lda.test(data, train = 10, bootstraps = 3)

lda_test_data <- list("procrustes" = procrustes_matrix, "species" = species_class, "morpho" = morpho_class, "random" = random_class, "lda_test" = lda_test)

save(lda_test_data, file = "../lda_test_data.Rda")