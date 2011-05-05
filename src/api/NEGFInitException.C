#include "UptInitException.h"

const char* ETBInitException::what(void) const throw()
{
   switch(_error)
   {
   case _ERR_ALLOC_ERR:
     return "allocation error"; 
	    
   case _ERR_ATM_MISMTCH:
     return "number of atom mismatch"; 	    

   case _ERR_EGV_VOID:
     return "empty eigenvectors"; 	    

   case _ERR_EGV_MISMTCH: 
     return "eigenvector mismatch"; 	    

   case _ERR_MAT_UNALLOC: 
     return "matrix not allocated"; 	    

   case _ERR_MAT_NOTSQRE: 
     return "matrix not square"; 	    

   case _ERR_NN_LIST:  
     return "nearest neighbour list error"; 	    

   case _ERR_HAM_UNMAT:   
     return "unrecognized material"; 	    

   case _ERR_HAM_UNPAIR:  
     return "unrecognized pair";

   case _ERR_HAM_ZEROLN:
     return "unsupported matrix with whole line of zero";  

   case _ERR_REF_ATMLIST:
     return "uncorrect reference atom list: check code";

   case _ERR_REF_NAME:
     return "state or coupling symbol not found in reference table (check etb db)";  

   case _ERR_ALLOY_TBBLK:
     return "TB block parent mismatch in alloy";

   case _ERR_DB_PAIR:
     return "wrong pair (check etb database)";
     
   case _ERR_INPUT:
     return "uncorrect input or database";

   case _ERR_DB_NOINT:
     return "integer not found in database";

   case _ERR_DB_NOBOOL:
     return "logical flag not found database";

   case _ERR_FILE_OPEN:
     return "file not found";

   case _ERR_OUTPUT:
     return "output error";

   default:
     return "general internal error";
   }
}


extern "C" {
  void throw_init_exception_(const int& errcode)
  {
     throw ETBInitException(errcode);
  }
}
