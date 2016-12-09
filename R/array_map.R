
  array_map = function( method, m, n ) {
    # array indexing from 1d to nd and nd to 1d  .. for higher dimensions, just follow the patterns
    # m are input coords
    # n are dimension sizes

    if (method=="2->1") { 
      return( m[,1] + (m[,2]-1)*n[1] ) 
    }
    
    if (method=="3->1") { 
      return( m[,1] + (m[,2]-1)*n[1] + (m[,3]-1)*n[1]*n[2] ) 
    }
    
    if (method=="4->1") { 
      return( m[,1] + (m[,2]-1)*n[1] + (m[,3]-1)*n[1]*n[2] + (m[,4]-1)*n[1]*n[2]*n[3] ) 
    }

    if ( method=="1->2" ) {
      j = m-1 # -1 converts to C-indexing
      x = j %%  n[1]
      j = j %/% n[1]
      y = j
      return( cbind(x,y)+1)  # +1 returns to R-indexing
    }

    if ( method=="1->3" ) {
      j = m-1 # -1 converts to C-indexing
      x = j %%  n[1]
      j = j %/% n[1]
      y = j %%  n[2]
      j = j %/% n[2]
      z = j
      return( cbind(x,y,z) + 1 ) # +1 returns to R-indexing
    }

    if ( method=="1->4" ) {
      j = m-1 # -1 converts to C-indexing
      x = j %%  n[1]
      j = j %/% n[1]
      y = j %%  n[2]
      j = j %/% n[2]
      z = j %%  n[3]
      j = j %/% n[3]
      a = j
      return( cbind(x,y,z,a) + 1 ) # +1 returns to R-indexing
    }
  }

