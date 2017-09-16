with Numerics, Numerics.Sparse_Matrices;
use  Numerics, Numerics.Sparse_Matrices;

with Ada.Text_IO; use Ada.Text_IO;
procedure Sparse_Test is
   use Real_IO, Int_IO;
   A  : Sparse_Matrix;
   C  : Sparse_Matrix;
   X  : Sparse_Vector := Sparse ((0.0, 2.0, 1.0));
   Y  : Sparse_Vector := Sparse ((3.0, 0.0));
   Z  : Sparse_Vector;
   M  : Real_Matrix (2 .. 3, 5 .. 6) := ((2.0, 3.0),
					 (5.0, 7.0));
   B  : Sparse_Matrix;
   
   U : Real_Vector (13 .. 14) := (1.0, 2.0);
   V : Real_Vector (1 .. 4);
begin
   
   Z := Zero (2);
   Print (Z);
   Z := One (2, 5);
   Print (Z);
   V (1 .. 2) :=  U;
   V (3 .. 4) := V (1 .. 2) - U;
   --  for Item of V loop
   --     Put (Item); New_Line;
   --  end loop;
   --  B := Eye (2); 
   --  Print (B);
   --  Put_Line ("---------------------------------------------------------");
   --  A := X * Y;   Print (A);
   --  Put_Line ("---------------------------------------------------------");
   --  C := A and B;   Print (C);
   --  Put_Line ("---------------------------------------------------------");
   --  B := Remove_1stN (C, 1);
   --  Print (B);
   --  Put_Line ("---------------------------------------------------------");
   --  Put_Line ("---------------------------------------------------------");
   
   --  Print (X);
   --  Put_Line ("---------------------------------------------------------");
   --  Print (Y);
   --  Put_Line ("---------------------------------------------------------");
   --  Z := X or Y;
   --  Print (Z);
   --  Put_Line ("---------------------------------------------------------");
   --  Z := Y or X;
   --  Print (Z);
   --  Y := Remove_1stN (Z, 2);
   --  Print (Y);
   --  Put_Line ("---------------------------------------------------------");
end Sparse_Test;
