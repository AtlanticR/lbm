
hivemod_attach = function( method, pointer ) {
   return( 
      switch( method, 
        bigmemory.ram=attach.big.matrix(pointer), 
        bigmemory.filebacked=attach.big.matrix(pointer), 
        ff=pointer )
    )
}

