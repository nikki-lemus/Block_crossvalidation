
block_extents4 <- function(ext) {  # we need this to work with more than 4 blocks
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
  

mop_comb <- function(block_list, vars, calculate_distance = TRUE) {  # we need to change this so it can work with more than four blocsk
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



