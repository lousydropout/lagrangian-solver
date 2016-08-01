with Ada.Text_IO, Numerics, Numerics.Sparse_Matrices, Chebyshev;
use  Ada.Text_IO, Numerics, Numerics.Sparse_Matrices, Chebyshev;

procedure Forward_AD.Test is
   use Real_Functions;
   
   
   function Lagrangian (Pos, Vel : in Real_Array; T : in Real) return AD_Type;
   function Lagrangian (Pos, Vel : in Real_Array; T : in Real) return AD_Type is
      Q : AD_Vector := Var (Pos);
      N : constant Nat := Pos'Length;
      Result : AD_Type := Q (1); -- ** 2;
   begin

      for K in 2 .. N loop
	 Result := Result + Q (K); -- * Q (K);
      end loop;
      return Result;
   end Lagrangian;
   
   
   N : constant Nat := 4;
   A : AD_Type;
   
   Q, V : Real_Array (1 .. N);

   L : constant Real := -1.7;
   R : constant Real :=  12.0;
   X : Real_Array  := Chebyshev_Gauss_Lobatto (N, L, R);
   D : Real_Matrix := Derivative_Matrix (N, L, R);
   Mat : Sparse_Matrix := Sparse (D);
   
   B : Sparse_Vector;
   
begin
   for I in X'Range loop
      Q (I) := Sin (X (I));
      Put_Line (Real'Image (X (I)) & "       " & Real'Image (Cos (X (I))));
   end loop;
   V := Q;
   A := Lagrangian (Q, V, 0.0);
   B := Grad (A);
   Print (B);
   
   Print (D * Sparse (Q));
   Put ("2nd: ");
   Print (Mat * Sparse (Q));
   Print (D * Sparse (Q));

   
   Put_Line ("length = " & Int'Image (Length (B)));
   null;
   
   
end Forward_AD.Test;
