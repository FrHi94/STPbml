##average data across parasite infected as well as non-infected positions (for comparison)

#data = a sc expression matrix you want to perform the analysis on 
#meta = the data you want to group your data by (e.g. clusters, celltypes, etc.)

hmat_pb <- function(data,meta) {
  
  gene_data <- data.frame(data , cluster = meta ,check.names = F)  
  average_data <- aggregate(.~cluster, gene_data, mean)
  cluster_name <- average_data[,1]
  average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
  rownames(average_data) <- cluster_name #you can also give your data individual cluster names 
  average_data <- t(average_data)
  phmat1 <- t(scale(t(average_data)))
  #create matrix for averaged expression-values between -1.5 and 1.5 (this can be adjusted) (clipping the values)
  #phmat1[phmat1> 1.5] <- 1.5
  #phmat1[phmat1 < -1.5] <- -1.5
  return(phmat1)
  
}

##heatmap to display average DEG generated in hmat_pb

library(pheatmap)
#hmat = the heatmap generated using the hmat_pb function 
#seed = the features you would like to include (charactervector of genes to display)
phmat_simple <- function(hmat, seed) {
  pheatmap(hmat[seed,], fontsize_row = 4,
           cellwidth =20,
           cellheight = 6,
           clustering_distance_cols = "correlation",
           cluster_cols = F,
           cluster_rows = T,
           #annotation_row = annotation_row,
           breaks = seq(-max(abs(hmat)), max(abs(hmat)), length.out = 1000), 
           color = colorRampPalette(c("#660156","#21908C", "white", "#FFF68F", "#FFD700"))(1001), 
           
           
           
  )
}


##Tree plot visualization of gene set enrichment 
library(gprofiler2)

#data = data to run query on 
#source = source one of GO:BP, KEGG, REAC etc. (see gprofiler documentation)
#nterms = the number of pathways to be displayed in the the plot 
#signifcant = whether only signinficant terms should be considered (default = TRUE)

golysis.tree <- function(data, source, nterms, significant = TRUE){
  dat <- list(data)
  res <- gost(query = dat, organism = "mmusculus", sources = source , significant = significant)
  res <- res$result
  res$term_go_id <- paste0(res$term_name,"::",res$term_id)
  #res <- subset(res,res$significant == T)
  g <- ggplot(head(res, nterms), aes(area= -log10(p_value), fill= term_name, label = term_go_id)) +
    geom_treemap() + 
    geom_treemap_text(colour = "white",
                      place = "centre",
                      size = 8,
                      grow = T)+
    scale_fill_manual(values = c("#7b2cbf", "#005f73", "#0a9396","#006ac7","#e85d04","#ca6702","#bb3e03","#d00000","#ae2012", "#9b2226","#b5179e", "#ffafcc", "#81b29a", "#007200", "#cbf3f0", "#40916c") ) +
    theme(legend.position = "none")
  return(list(res,g))
}

##FUNCTION correlation for distance and gene expression (cor_dist)


# dist = dataframe with pixel-cordinates (XY) as rownames and the first column containing distances in pixels and the second column containing distances in µm 
# thrs = numeric value for threshold to set as maximal distance include for correlations
# object = STUtility object (SCT normalized) if you used a different normalization you need to change line 72 accordingly

cor_dist <- function(dist, thrs, object) {
  colnames(dist) <- c("dist_px", "dist_um") 
  dist.1 <- na.omit(dist) #remove spots without distance measurement to IHSs
  dist.1 <- dist.1[dist.1[,2] <= thrs,] #set distance threshold for correlation analysis 
  dat1 <- as.matrix(t(object@assays$SCT@scale.data)) # create a transposed dataframe (or matrix) for expression data
  dat2 <- dat1[rownames(dat1) %in% rownames(dist.1),] #intersect selection 
  dist3 <- dist.1[rownames(dist.1) %in% rownames(dat2),] #intersect in both directions          
  
  
  cor.roi_gex <- sapply(1:ncol(dat2), function(x){ #sapply function to create list with values of interest and return list of dataframes
    res <- corr.test(dist3[,2], dat2[,x], 
                     method = "pearson", adjust = "bonferroni")
    return(list(data.frame(correlation_coeficent = res$r , p.value = res$p, p.adj = res$p.adj, row.names = colnames(dat2)[x])))
    
  })
  
  cor_IHS <- do.call(rbind.data.frame, cor.roi_gex) #bin list entries to new data frame
  
  cor_IHS <- cor_IHS[order(-cor_IHS$correlation_coeficent),]
  
  return(cor_IHS)
  
}

## FUNCTION expression by distance (expr_by_distance)

#object = STUTility/Seuratobject 
#seed = genes displayed in expression by distance plots
#roi_dist = distances across coordinates of regions of interest 
#plotgenes = which genes of the seed genes to show in the plot
#center = character stating what the center of the ROI (for x-axis label)
#ncol = number of columns to show in the plot 
#nrow = number of rows to show in the plot 
#thrs = distance threshold around ROI to display gene expression 
#condition = conditions to be displayed 
#show.condition = whether the graphs should be split by condition
#delta = whether to center the data around 0 when split by condition 

expr_by_dist <- function(object, seed, roi_dist, plotgenes = c(1:length(seed)), center, ncol = 5, nrow = 4, thrs = 500, condition, show.condition = T, delta = T) {
  pb.mat <- as.matrix(GetAssayData(object[["SCT"]], slot = "data"))
  pb.mat[1:5,1:5]
  
  #get all spots with at least one count for any of the genes of interest
  pb.s <- pb.mat[seed,]
  dist_parasite1 <- subset(roi_dist, rownames(roi_dist) %in% colnames(pb.s))
  condition1 <- condition[which(rownames(condition) %in% colnames(pb.s)),]
  pb.goi <- as.data.frame(t(pb.s))
  pb.goi <- subset(pb.goi, rownames(pb.goi) %in% rownames(dist_parasite1))
  pb.goi <- subset(pb.goi, rownames(pb.goi) %in% rownames(condition1))
  condition1 <- subset(condition1, rownames(condition1) %in% rownames(pb.goi))  
  
  pb.goi.dist <- cbind(pb.goi, dist_parasite1, condition1)
  colnames(pb.goi.dist)[(ncol(pb.goi.dist)-2)] <- "dist_um" #change colnames of distances so they can be used for plotting
  colnames(pb.goi.dist)[(ncol(pb.goi.dist)-3)] <- "dist_px"
  pb.goi.dist[,c(ncol(pb.goi.dist)-2)] <- as.numeric(pb.goi.dist[,c(ncol(pb.goi.dist)-2)])
  pb.goi.dist[,ncol(pb.goi.dist)-1] <- as.factor(pb.goi.dist[,ncol(pb.goi.dist)-1])
  
  pb.goi.dist.p <- subset(pb.goi.dist, pb.goi.dist$dist_um <= thrs)
  
  condition$condition <- as.factor(condition$condition)
  
  if (delta) {
    list_of_dfs <- split(pb.goi.dist.p, pb.goi.dist.p$condition)
    for (i in seq_along(list_of_dfs)) {
      curr_df <- list_of_dfs[[i]]
      curr_condition <- unique(curr_df$condition)
      if (nrow(curr_df) == 0) {
        next
      }
      for (col in colnames(curr_df)[1:c(ncol(curr_df)-4)]) {
        fit <- loess(curr_df[, col] ~ dist_um, data = curr_df)
        highest_fitted <- fit$fitted[which.min(fit$x)]
        curr_df[, col] <- curr_df[, col] - highest_fitted
      }
      list_of_dfs[[i]] <- curr_df
    }
    pb.goi.dist.p <- do.call(rbind, list_of_dfs)
  }
  
  pb.goi.dist.pl <- lapply(colnames(pb.goi.dist.p[1:c(ncol(pb.goi.dist.p)-4)]), function(x){
    ifelse(show.condition, 
           p <- ggplot(pb.goi.dist.p , aes(x= dist_um, y= pb.goi.dist.p[,x], colour = condition)) + 
             geom_smooth(aes(fill = condition),method = "loess", 
                         n = 10,
                         span = 2) + 
             ylab(ifelse(delta, paste0("Δ ", x), x)) + 
             xlab(paste0("distance to"," ",center,"[µm]")) +
             theme_classic(),
           
           p <- ggplot(pb.goi.dist.p , aes(x= dist_um, y=  pb.goi.dist.p[,x])) + 
             geom_smooth(method = "loess", 
                         n = 10, span = 2) + 
             ylab(ifelse(delta, paste0("Δ ",x), x)) + 
             xlab(paste0("distance to"," ",center,"[µm]")) +
             theme_classic()
    )
    
    return(p)
  })
  
  
  p.l <- ggarrange(plotlist = pb.goi.dist.pl, ncol = ncol, nrow = nrow)
  
  return(list(out1 = print(head(pb.goi.dist)), out2 = p.l))
  
  
}


## FUNCTION subset DG sets to contain unique ones with highest expression values}

# df = results from DGEA

keep_max_logFC <- function(df) {
  # Sort the data frame by gene and logFC in descending order
  df_sorted <- df[order(df$gene, -df$avg_log2FC),]
  
  # Initialize an empty vector to hold the row indices to keep
  indices_to_keep <- c()
  
  # Initialize a variable to hold the index of the first row with the new gene
  first_gene_index <- 1
  
  # Loop through the data frame
  for (i in 1:nrow(df_sorted)) {
    # If current gene is the same as previous gene
    if (i > 1 && df_sorted$gene[i] == df_sorted$gene[i-1]) {
      # Skip row if its index is not the highest value seen so far for that gene
      if (i != first_gene_index) {
        next
      }
      # If current gene is not the same as previous gene
    } else {
      # Add index of first row with new gene to vector of indices to keep
      indices_to_keep <- c(indices_to_keep, first_gene_index)
      # Update index of first row with new gene
      first_gene_index <- i
    }
  }
  
  # Add index of last row to vector of indices to keep
  indices_to_keep <- c(indices_to_keep, first_gene_index)
  
  # Create new data frame with only the rows whose indices are in the vector of indices to keep
  df_unique <- df_sorted[indices_to_keep,]
  
  df_unique <- df_unique[order(-df$cluster),]
  
  return(df_unique)
}


##HEATPLOT 

#go.table = dataframe with 1 column with GO-term description and one columen containing genes involved in the process
#logFC = dataframe containing DEG data including genes 

heatplot <- function(go.table, logFC){
  

  #make data frame with 1 column as the description and one column with all the genes
  nrow(go.table)
  l <- lapply(1:nrow(go.table), function(x){
    df <- data.frame(genes = do.call(cbind, strsplit(go.table[x,16], ",", fixed = T)), goterm= paste(go.table[x,11]))
  })
  library(data.table)
  df <- rbindlist(l)
  #now we can add the logFCs from our genelist
  logFC <- as.data.frame(subset(logFC, names(logFC) %in%  df$genes))
  logFC$genes <- rownames(logFC)
  colnames(logFC) <- c("logFC", "genes")
  head(logFC)
  #now combine the data 
  library(dplyr)
  df <- left_join(df, logFC, by = "genes")
  df$genes <- as.factor(df$genes)
  
  order <- rev(c(go.table$term_name))
  df$goterm <- factor(df$goterm, levels = order)
  
  
  p <- ggplot(df, aes(reorder(genes, -logFC), goterm)) +
    geom_tile(aes(fill = logFC), colour = "white") + 
    scale_fill_gradientn(colors = colorRampPalette(c("#D05D37", "#E9967A", "white",  "#76EEC6", "#20DF9F"))(40)) + 
    #  coord_flip()+
    xlab("gene")+ 
    ylab("enriched term (KEGG)") +
    theme(panel.background = element_rect(fill = "white"), axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size=0.2, color = "black"))
  #dev.off()
  
  return(p)
  
}

##Proportion by distance

#object = STUTility/Seuratobject 
#seed = celltypes displayed in expression by distance plots
#roi_dist = distances across coordinates of regions of interest 
#plotgenes = which celltypes of the seed genes to show in the plot
#center = character stating what the center of the ROI (for x-axis label)
#ncol = number of columns to show in the plot 
#nrow = number of rows to show in the plot 
#thrs = distance threshold around ROI to display proportion values
#condition = conditions to be displayed 
#show.condition = whether the graphs should be split by condition
#delta = whether to center the data around 0 when split by condition 

prop_by_dist <- function(object, type, roi_dist, plottypes = c(1:length(types)), center, ncol = 5, nrow = 4, thrs = 500, condition, show.condition = T, delta = T) {
  
  pb.goi.dist <- data.frame(distance = object[[]][colnames(object@meta.data) %in% roi_dist], object[[]][colnames(object@meta.data) %in% type], condition = object$condition)
  colnames(pb.goi.dist)[1] <- "dist_um" #change colnames of distances so they can be used for plotting
  
  pb.goi.dist.p <- subset(pb.goi.dist, pb.goi.dist$dist_um <= thrs)
  
  pb.goi.dist.p$condition <- as.factor(pb.goi.dist.p$condition)
  
  if (delta) {
    list_of_dfs <- split(pb.goi.dist.p, pb.goi.dist.p$condition)
    for (i in seq_along(list_of_dfs)) {
      curr_df <- list_of_dfs[[i]]
      curr_condition <- unique(curr_df$condition)
      if (nrow(curr_df) == 0) {
        next
      }
      for (col in colnames(curr_df)[2:c(ncol(curr_df)-1)]) {
        fit <- loess(curr_df[, col] ~ dist_um, data = curr_df)
        highest_fitted <- fit$fitted[which.min(fit$x)]
        curr_df[, col] <- curr_df[, col] - highest_fitted
      }
      list_of_dfs[[i]] <- curr_df
    }
    pb.goi.dist.p <- do.call(rbind, list_of_dfs)
  }
  
  pb.goi.dist.pl <- lapply(colnames(pb.goi.dist.p[2:c(ncol(pb.goi.dist.p)-1)]), function(x){
    ifelse(show.condition, 
           p <- ggplot(pb.goi.dist.p , aes(x= dist_um, y= pb.goi.dist.p[,x], colour = condition)) + 
             geom_smooth(aes(fill = condition),method = "loess", 
                         n = 10,
                         span = 2) + 
             ylab(ifelse(delta, paste0("Δ ", x), x)) + 
             xlab(paste0("distance to"," ",center,"[µm]")) +
             theme_classic(),
           
           p <- ggplot(pb.goi.dist.p , aes(x= dist_um, y=  pb.goi.dist.p[,x])) + 
             geom_smooth(method = "loess", 
                         n = 10, span = 2) + 
             ylab(ifelse(delta, paste0("Δ ",x), x)) + 
             xlab(paste0("distance to"," ",center,"[µm]")) +
             theme_classic()
    )
    
    return(p)
  })
  
  
  p.l <- ggarrange(plotlist = pb.goi.dist.pl, ncol = ncol, nrow = nrow)
  
  return(list(out1 = print(head(pb.goi.dist)), out2 = p.l))
  
  
}

##Correlation by distance 

#object = STUTility/Seuratobject 
#seed = genes displayed in expression by distance plots
#roi_dist = distances across coordinates of regions of interest 
#plotgenes = which genes of the seed genes to show in the plot
#center = character stating what the center of the ROI (for x-axis label)
#ncol = number of columns to show in the plot 
#nrow = number of rows to show in the plot 
#thrs = distance threshold around ROI to display gene expression 
#condition = conditions to be displayed 
#show.condition = whether the graphs should be split by condition
#delta = whether to center the data around 0 when split by condition 

cor_dist_prop <- function(object, types,dist, thrs) {
  
  dist_df <- data.frame(dist = object[[]][colnames(object@meta.data) %in% dist], object[[]][colnames(object@meta.data) %in% types])
  
  dist.1 <- na.omit(dist_df) #remove spots without distance measurement to IHSs
  dist.1 <- dist.1[dist.1[,1] <= thrs,] #set distance threshold for correlation analysis 
  
  cor.roi_gex <- sapply(2:ncol(dist.1), function(x){ #sapply function to create list with values of interest and return list of dataframes
    res <- corr.test(dist.1[,1], dist.1[,x], 
                     method = "pearson", adjust = "bonferroni")
    return(list(data.frame(correlation_coeficent = res$r , p.value = res$p, p.adj = res$p.adj, row.names = colnames(dist.1)[x])))
    
  })
  
  cor_IHS <- do.call(rbind.data.frame, cor.roi_gex) #bin list entries to new data frame
  
  cor_IHS <- cor_IHS[order(-cor_IHS$correlation_coeficent),]
  
  return(cor_IHS)
  
}

