
# helper function to make 4 blocks in a region
block_extents4 <- function(ext) {  
  blist <- list(
    ext(c(ext[1], mean(c(ext[1], ext[2])), mean(c(ext[3], ext[4])), 
          ext[4])),
    ext(c(mean(c(ext[1], ext[2])), ext[2], mean(c(ext[3], ext[4])), 
          ext[4])),
    ext(c(ext[1], mean(c(ext[1], ext[2])), ext[3], 
          mean(c(ext[3], ext[4])))),
    ext(c(mean(c(ext[1], ext[2])), ext[2], ext[3], 
          mean(c(ext[3], ext[4]))))
  )
  
  names(blist) <- paste0("block_", 1:length(blist))
  
  return(blist)
}
  

# automated mop analysis comparing one block of interest vs three of reference: 
# makes all of the combinations at a time
mop_comb <- function(block_list, vars, calculate_distance = TRUE) {  
  bcomb <- combn(1:4, 3)
  
  mop_list <- lapply(1:ncol(bcomb), function (x) {
    b3 <- union(union(vect(block_list[[bcomb[1, x]]]), 
                      vect(block_list[[bcomb[2, x]]])), 
                vect(block_list[[bcomb[3, x]]]))
    
    crs(b3) <- crs(vars)
    
    ## run MOP
    mp <- mop(m = crop(vars, b3, mask = TRUE), g = vars, type = "detailed", 
          calculate_distance = calculate_distance, where_distance = "all")
    
    c(mp, b3_comb = b3)
  })
  
  names(mop_list) <- sapply(1:ncol(bcomb), function(x) {
    paste0("Comb_", paste(bcomb[, x], collapse = ""))
  })
  
  return(mop_list)
}


# helper function to plot results comparing distances resulting from mop_comb
summarize_mop_comb <- function(mop_comb_res, block_list) {
  
  bnum <- length(block_list):1
  res <- lapply(1:length(mop_comb_res), function(x) {
    ## areas 
    ref <- mop_comb_res[[x]]$b3_comb
    int <- block_list[[bnum[x]]]
    
    ## mask to areas
    disarea_ref <- crop(mop_comb_res[[x]]$mop_distances, ref, mask = TRUE)
    disarea_int <- crop(mop_comb_res[[x]]$mop_distances, int)
    
    ## extract values
    disval_ref <- as.data.frame(disarea_ref, cell = FALSE)
    disval_int <- as.data.frame(disarea_int, cell = FALSE)
    disval_all <- as.data.frame(mop_comb_res[[x]]$mop_distances, cell = FALSE)
    
    ## organize as a table
    bxvals <- rbind(disval_ref, disval_int, disval_all)
    bxvals$Areas <- c(rep("Reference", nrow(disval_ref)), 
                      rep("Left_out", nrow(disval_int)),
                      rep("All", nrow(disval_all)))
    
    bxvals$Areas <- as.factor(bxvals$Areas)
    
    colnames(bxvals)[1] <- "Distance"
    bxvals$Distance <- as.numeric(bxvals$Distance)
    
    return(bxvals)
  })
  
  names(res) <- names(mop_comb_res)
  
  return(res)
}



# 
plot_summary <- function(summarize_mop_comb, violin = FALSE, outliers = FALSE) {
  
  # plotting settings
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  
  # change settings
  par(mfrow = c(2, 2), mar = c(4, 4, 0.5, 0.5))
  
  lnam <- names(summarize_mop_comb)
  
  # plot in loop
  for (i in 1:length(summarize_mop_comb)) {
    if (violin) {
      ggplot2::ggplot(data = summarize_mop_comb[[i]], 
                      aes(x = Areas, y = Distance)) + 
        geom_violin(fill = "gray75", trim = outliers) +
        theme(
          panel.background = element_rect(fill = "white"),     # White plotting background
          plot.background = element_rect(fill = "white"),      # White full plot background
          panel.grid = element_blank(),                        # Remove grid lines
          axis.line = element_line(color = "black")            # Black axis lines
        )
    } else {
      if (i %in% c(1, 3)) {
        yt <- "Distance"
      } else {
        yt <- ""
      }
      
      if (i >= 3) {
        xt <- "Areas"
      } else {
        xt <- ""
      }
      
      boxplot(Distance ~ Areas, data = summarize_mop_comb[[i]], 
              outline = outliers, frame = FALSE, ylab = yt, xlab = xt, las = 1)
      legend("topright", legend = lnam[i], bty = "n")
      box(bty = "l")
    }
  }
}



boxplot(Distance ~ Areas, data = bxvals, 
        outline = outliers, frame = FALSE)


