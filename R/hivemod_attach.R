
hivemod_attach = function( method, pointer ) {
  #\\ generic method to attach data pointer from bigmemory or ff 
   return( 
      switch( method, 
        bigmemory.ram=attach.big.matrix(pointer), 
        bigmemory.filebacked=attach.big.matrix(pointer), 
        ff=pointer )
    )
}

