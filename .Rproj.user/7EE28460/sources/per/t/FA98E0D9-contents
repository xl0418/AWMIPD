#' Calculate AWMIPD values for a given phylogenetic and the corresponding abundances
#'
#' @param tree A phylogenetic tree
#' @param abundance_data A set of abundances corresponding to the tips of the tree
#' @param geographic_dis Geographic information of species distribution
#' @param subset_species Subsetting the phylogenetic by providing the subset of the species
#' @return The AWMIPD values
#' @examples
#' AWMIPD(1, 1)
#' @export

AWMIPD <-
  function(tree,
           abundance_data,
           geographic_dis = NULL,
           subset_species = NULL,
           a = 1,
           pdmode = 'exp') {
    # Check the tree class.
    if (class(tree) != "phylo") {
      return(print('Provide a phylo tree!'))
    }
    # Check whether the subset species is really the subset of the phylogeny.
    if (all(!(subset_species %in% tree$tip.label))) {
      return(print('One or more species in the subset are not in the phylogeny.'))
    }

    if (is.null(subset_species)) {
      distance_matrix <- ape::cophenetic.phylo(tree)
    } else {
      tree <- keep.tip(tree, subset_species)
      distance_matrix <- ape::cophenetic.phylo(tree)
    }

    species_name_order <- rownames(distance_matrix)
    abundance_index <-
      match(species_name_order, abundance_data$species)
    abundance <- abundance_data$abundance[abundance_index]
    sp_order <- abundance_data$sp[abundance_index]

    # phylogenetic distance
    if (pdmode == 'inv') {
      ID = 1 / distance_matrix
      diag(ID) = 1
      AD.matrix = sweep(ID, MARGIN = 2, as.matrix(abundance), `*`)
    } else if (pdmode == 'exp') {
      ID = exp(- a * distance_matrix)
      AD.matrix = sweep(ID, MARGIN = 2, as.matrix(abundance), `*`)
    }

    total.dvalues = rowSums(AD.matrix) * as.matrix(abundance)
    awmipd_real_value = total.dvalues/sum(total.dvalues)
    D.normalized = (total.dvalues - min(total.dvalues)) / (max(total.dvalues) -
                                                             min(total.dvalues))
    awmipd_df <-
      data.frame(
        species = species_name_order,
        sp = sp_order,
        awmipd = D.normalized,
        awmipd_real = awmipd_real_value,
        abundance = abundance
      )


    if (is.null(geographic_dis)) {
      return(list(awmipd = awmipd_df, subtree = tree))
    } else {
      order_in_geo_data <- match(geographic_dis$sp , sp_order)
      geographic_dis['awmipd'] <- D.normalized[order_in_geo_data]

      # low blue, high red
      plot_awmipd <-
        ggplot2::ggplot(geographic_dis, aes(gx, gy)) + geom_point(aes(color = awmipd), size = 0.5, alpha = 0.8) +
        theme(
          legend.position = 'right',
          # axis.text = element_blank(),
          # axis.ticks = element_blank(),
          axis.line.x = element_line(color = "black", size = 1),
          axis.line.y = element_line(color = "black", size = 1),
          panel.background = element_blank()
        ) +
        xlab("") + ylab("") + scale_color_gradient2(
          low = "#005CAF",
          mid = 'green',
          high = "#D0104C",
          midpoint = 0.5,
          name = "AWMIPD"
        )

      return(list(
        awmipd = awmipd_df,
        plot = plot_awmipd,
        subtree = tree
      ))
    }


  }
