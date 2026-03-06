
#Rebalance a hierarchy of values to ensure that parent nodes are not less than the sum of their children.
safe_hierarchy_rebalance <- function(mat, parents) {
  
  for (col in seq_len(ncol(mat))) {
    
    root_index <- which(parents == 0)
    mat[root_index, col] <- 100
    
    repeat {
      
      changed <- FALSE
      
      for (i in seq_along(parents)) {
        
        parent <- parents[i]
        if (parent == 0) next
        
        children <- which(parents == parent)
        children_sum <- sum(mat[children, col])
        
        if (children_sum > mat[parent, col] && children_sum > 0) {
          
          scale_factor <- mat[parent, col] / children_sum
          mat[children, col] <- mat[children, col] * scale_factor
          
          changed <- TRUE
        }
      }
      
      if (!changed) break
    }
    
  }
  
  mat
}