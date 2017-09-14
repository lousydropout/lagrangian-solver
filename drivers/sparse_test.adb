with Numerics, Numerics.Sparse_Matrices;
use  Numerics, Numerics.Sparse_Matrices;

with Ada.Text_IO; use Ada.Text_IO;
procedure Sparse_Test is
   
   A  : Sparse_Matrix;
   B  : Sparse_Matrix;
   C  : Sparse_Matrix;
   X  : Sparse_Vector := Sparse ((0.0, 0.0, 2.0, 1.0));
   Y  : Sparse_Vector := Sparse ((3.0, 4.0, 0.0));
   Z  : Sparse_Vector;
begin
   
   --  B := Eye (4); Print (B);
   --  Put_Line ("---------------------------------------------------------");
   --  A := X * Y;   Print (A);
   --  Put_Line ("---------------------------------------------------------");
   --  C := A + B;   Print (C);
   --  Put_Line ("---------------------------------------------------------");
   --  B := Remove_1stN (C, 1);
   --  Print (B);
   --  Put_Line ("---------------------------------------------------------");
   --  Put_Line ("---------------------------------------------------------");
   
   Print (X);
   Put_Line ("---------------------------------------------------------");
   Print (Y);
   Put_Line ("---------------------------------------------------------");
   --  Y := Remove_1stN (X, 2);
   --  Print (Y);
   --  Put_Line ("---------------------------------------------------------");
   Z := X and Y;
   Print (Z);
   Put_Line ("---------------------------------------------------------");
   Z := Y and X;
   Print (Z);
end Sparse_Test;
