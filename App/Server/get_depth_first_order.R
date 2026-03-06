# Helper function to get depth-first order of nodes based on parent-child relationships by recursively traversing the tree structure.
get_depth_first_order <- function(node_ids, parent_index) {
  
  children_map <- split(seq_along(parent_index), parent_index) # list where names are parent indices and values are vectors of child indices
  
  traverse <- function(idx) {
    
    current <- node_ids[idx] # start with the current node ID
    child_idx <- children_map[[as.character(idx)]] # get child indices of the current node
    
    if (is.null(child_idx)) { # if no children, return current node ID -> leaf node
      return(current)
    }
    
    # depth-first order = biological order
    c(
      current,
      unlist(lapply(child_idx, traverse)) # children are traversed recursively in the same way
    )
  }
  
  # roots = nodes without parent
  roots <- which(parent_index == 0)
  unlist(lapply(roots, traverse))
}