## Numerical corrections

#1. Fix zero reappearance: ensures that once a clone appears (value > 0), it cannot have zero or NA values in subsequent timepoints.
# This is important for fishplot visualization to maintain the continuity of clones over time, even if they become very small.
# steps : cleans numerical issues, detects first appearance and prevents disapperance after appearance
fix_zero_reappearance <- function(mat, eps = 1e-4) {
  
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0 #cleans numerical issues
  
  for (i in seq_len(nrow(mat))) {
    
    row_vals <- mat[i, ]
    
    first_pos <- which(row_vals > 0)[1] #detects first appearance
    
    if (!is.na(first_pos) && first_pos > 1) {
      mat[i, first_pos - 1] <- eps
    }
    
    appeared <- FALSE
    
    for (t in seq_len(ncol(mat))) {
      
      val <- mat[i, t]
      
      if (!is.na(val) && val > 0) {
        appeared <- TRUE
      }
      
      if (appeared && (is.na(val) || val == 0)) { #prevents disapperance after appearance
        mat[i, t] <- eps
      }
    }
  }
  
  mat
} # Clones do not truly disappear between measurements, they may fall below detection but still exists biologically
# Necessary because fishplot requires : continuous trajectories and no abrupt zero-positive transitions

#2. Fix parent-child conflicts: ensures that at any timepoint, the sum of child clone frequencies does not exceed the parent clone frequency.
# Ensures biological consistency
fix_parent_child_conflicts <- function(mat, parents) {
  
  n_time <- ncol(mat)
  n_clones <- nrow(mat)
  
  for (t in seq_len(n_time)) {
    
    for (i in seq_len(n_clones)) {
      
      children <- which(parents == i)
      if (length(children) == 0) next # if children exceed parent -> rescale children to fit parent frequency
      
      parent_val <- mat[i, t]
      sum_children <- sum(mat[children, t])
      
      if (sum_children > parent_val && sum_children > 0) {
        
        ratio <- parent_val / sum_children
        mat[children, t] <- mat[children, t] * ratio
        
      }
    }
  }
  mat
} # A subclone cannot exceed its ancestral clone -> This enforces evolutionary logic

#3. Fix independent clones: ensures that the sum of frequencies of independent clones does not exceed 100% at any timepoint.
# If multiple clones have no parent (parent ==0) -> rescale them proportionally
fix_independent_clones <- function(mat, parents) {
  
  root_idx <- which(parents == 0)
  if (length(root_idx) <= 1) return(mat)
  
  for (t in seq_len(ncol(mat))) {
    s <- sum(mat[root_idx, t])
    
    if (s > 100 && s > 0) {
      mat[root_idx, t] <- mat[root_idx, t] * (100 / s)
    }
  }
  
  mat
} # Total tumor fraction cannot exceed 100%

#4. Enforce all fishplot constraints: applies both parent-child and independent clone corrections iteratively until no violations remain.
# Iteratively enforce all fishplot mathematical constraints
enforce_fishplot_constraints <- function(mat, parents, tol = 1e-8) {
  
  mat <- as.matrix(mat)
  n_time <- ncol(mat)
  n_clones <- nrow(mat)
  
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0
  
  for (t in seq_len(n_time)) {
    
    repeat {
      
      changed <- FALSE
      
      # PARENT / CHILD constraint
      for (i in seq_len(n_clones)) {
        
        children <- which(parents == i)
        if (length(children) == 0) next
        
        parent_val <- mat[i, t]
        sum_children <- sum(mat[children, t])
        
        if (sum_children > parent_val + tol && sum_children > 0) {
          
          ratio <- parent_val / sum_children
          mat[children, t] <- mat[children, t] * ratio
          changed <- TRUE
        }
      }
      
      #  ROOT constraint
      root_idx <- which(parents == 0)
      if (length(root_idx) > 0) {
        
        root_sum <- sum(mat[root_idx, t])
        
        if (root_sum > 100 + tol && root_sum > 0) {
          
          mat[root_idx, t] <- mat[root_idx, t] * (100 / root_sum)
          changed <- TRUE
        }
      }
      
      if (!changed) break
    }
  }
  
  mat
}