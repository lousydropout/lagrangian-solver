with Numerics, Numerics.Sparse_Matrices, Ada.Text_IO;
use  Numerics, Numerics.Sparse_Matrices, Ada.Text_IO;
with Ada.Numerics.Generic_Real_Arrays;

procedure Sparse_Test is
   use Real_IO, Int_IO;
   package RA is new Ada.Numerics.Generic_Real_Arrays (Real); use RA;
   
   --  U  : RA.Real_Matrix := (1 => (1 => 1.0), 
   --  			   2 => (1 => 2.0), 
   --  			   3 => (1 => 0.0),
   --  			   4 => (1 => 4.0),
   --  			   5 => (1 => 0.0),
   --  			   6 => (1 => -3.2));
   --  V  : RA.Real_Matrix := (1 => (6.0, 0.0, 7.0, 1.0, 8.0,  4.7));
   X  : Sparse_Vector := Sparse ((1.0, 2.0));
   Y  : Sparse_Vector := Sparse ((6.0, 0.0));
   A  : Sparse_Matrix;
   B  : Sparse_Matrix;
   C  : Sparse_Matrix;
   --  D : RA.Real_Matrix := U * V;
   --  E : RA.Real_Matrix := Transpose (D);
   --  F : RA.Real_Matrix (E'Range (1), E'Range (2));
begin
   
   A := X * Y; Print (A);
   B := X * X; Print (B);
   
   Print (A - B);
   Print (B - A);
   --  for I in 1 .. 10_000 loop
   --     --  C := A * B;
   --     F := D * E;
   --  end loop;
   
end Sparse_Test;
