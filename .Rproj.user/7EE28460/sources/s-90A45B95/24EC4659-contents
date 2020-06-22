# Load phylogenetic tree; Input 1
tree <- read.tree('./inst/control.tre')
plot(tree)
# Read the BCI data
mydatapath = './inst/'
if (!exists('bci.full7')) {
  attach(paste(mydatapath, 'bci.full7.rdata', sep = ''))
}

# Extract the species mnemonic and the geographic info; Input 2
bci_alive_data <-
  subset(bci.full7,
         status = 'A',
         select = c('sp', 'gx', 'gy', 'dbh'))
sp_code <- unique(bci_alive_data[, 1])

# Get full species names
load('./inst/bci.spptable.Rdata')

sp_name <- NULL
genus_name <- NULL
for (i in c(1:length(sp_code))) {
  index <- which(bci.spptable$sp == sp_code[i])
  if (length(index) == 0) {
    sp_name <-
      c(sp_name, 0)
    genus_name <-
      c(genus_name, 0)
  } else {
    sp_name <-
      c(sp_name, bci.spptable$Species[index])
    genus_name <-
      c(genus_name, bci.spptable$Genus[index])
  }
}

species_name <- paste0(genus_name, '_', sp_name)
sp_code_name <- cbind(sp_code, species_name)

recognized_species <- which(species_name %in% tree$tip.label)
recognized_code_name <- sp_code_name[recognized_species,]

# Calculate the abundance of each species; Input 3
abundance <- NULL
for(i in c(1:nrow(recognized_code_name))) {
  abundance_species <- length(which(bci_alive_data[, 1] == recognized_code_name[i, 1]))
  abundance <- c(abundance, abundance_species)
}

# Wrap the sp code, species names and the abundance data into a data frame.
abundance_data <- data.frame(sp = recognized_code_name[,1], species = recognized_code_name[,2], abundance = abundance)

# Subsetting the species; Input 4
# not_in_species <- subset(species_name, !(species_name %in% tree$tip.label))
species_in <- subset(species_name, (species_name %in% tree$tip.label))

# abundance_data$abundance <- 10
# Calculate AWMIPD and plot
test <- AWMIPD(tree = tree, abundance_data = abundance_data,
               geographic_dis = bci_alive_data, subset_species = species_in,
               a = 0.01)

awmipd_value <- test[[1]]
awmipd_distribution_plot <- test[[2]]
subtree <- test[[3]]


test_inv <- AWMIPD(tree = tree, abundance_data = abundance_data,
               geographic_dis = bci_alive_data, subset_species = species_in, pdmode = 'inv')
plot_inv <- test_inv[[2]]
