with Numerics;
use  Numerics;
package Chebyshev is
   
   
   function Chebyshev_Gauss_Lobatto (N : in Nat;
				     L : in Real := -1.0;
				     R : in Real :=  1.0) return Real_Array;

   function Derivative_Matrix (N    : in Nat;
			       L, R : in Real) return Real_Matrix;
   
   function CGL_Transform (F : in Real_Array) return Real_Array;
   
   function Interpolate (A : in Real_Array;
			 X : in Real;
			 L : in Real := -1.0;
			 R : in Real :=  1.0) return Real;
   
end Chebyshev;
