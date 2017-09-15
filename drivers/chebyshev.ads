with Numerics;
use  Numerics;
package Chebyshev is
   
   function Chebyshev_Gauss_Lobatto (N : in Nat;
				     L : in Real := 0.0;
				     R : in Real := 1.0) return Real_Vector;

   function Derivative_Matrix (N : in Nat;
			       L : in Real := 0.0;
			       R : in Real := 1.0) return Real_Matrix;
   
   procedure CGL (D :    out Real_Matrix;
		  X :    out Real_Vector;
		  N : in     Nat;
		  L : in     Real	 := 0.0;
		  R : in     Real	 := 1.0);

   
   function CGL_Transform (F : in Real_Vector) return Real_Vector;
   
   function Interpolate (A : in Real_Vector;
			 X : in Real;
			 L : in Real := 0.0;
			 R : in Real := 1.0) return Real;
   
end Chebyshev;
