with Numerics;
use  Numerics;
package Chebyshev is
   
   
   function Chebyshev_Gauss_Lobatto (N : in Nat;
				     L : in Real := -1.0;
				     R : in Real :=  1.0) return Real_Array
     with Pre => (N >= 1);

   function Derivative_Matrix (N    : in Nat;
			       L, R : in Real) return Real_Matrix
     with Pre => (N >= 1);
   
   
end Chebyshev;
