with Numerics, Numerics.Sparse_Matrices, Chebyshev;
use  Numerics, Numerics.Sparse_Matrices, Chebyshev;
generic
   K    : in Nat;
package Dense_AD.Integrator is
   
   type Variable is record
      X : Vector;
      T : Real;
   end record;
   
   type Control_Type is record
      Dt  : Real := 1.0;
      Eps : Real := 1.0e-10;
      Err : Real := 1.0;
   end record;
   
   function Collocation (Lagrangian : not null access 
			   function (X : Vector) return AD_Type;
			 Var        : in     Variable;
			 Control    : in out Control_Type) return Real_Vector;
   procedure FJ (Lagrangian : not null access 
		    function (X : Vector) return AD_Type;
		  Var     : in     Variable;
		  Control : in     Control_Type;
		  Q       : in     Real_Vector;
		  F       :    out Sparse_Vector;
		  J       :    out Sparse_Matrix);
   
   
   procedure Print_Lagrangian (X : in Vector;
			       Lagrangian : not null access
				 function (X : Vector)  return AD_Type);

   
   procedure Setup;
   
private
   
   Top_Left     : constant Real_Matrix := ((1.0, 0.0), (0.0, 0.0));
   Top_Right    : constant Real_Matrix := ((0.0, 1.0), (0.0, 0.0));
   Bottom_Left  : constant Real_Matrix := ((0.0, 0.0), (1.0, 0.0));
   Bottom_Right : constant Real_Matrix := ((0.0, 0.0), (0.0, 1.0));
   
   Grid   : constant Real_Vector := Chebyshev_Gauss_Lobatto (K, 0.0, 1.0);
   Der    : constant Real_Matrix := Derivative_Matrix (K, 0.0, 1.0);
   Half_N : constant Nat := N / 2;
   NK     : constant Nat := N * K;
   Mat_A, Mat_B, Mat_C, Mat_D : Sparse_Matrix;
   
end Dense_AD.Integrator;
