with Numerics, Numerics.Sparse_Matrices;
use  Numerics, Numerics.Sparse_Matrices;

with Ada.Text_IO; use Ada.Text_IO;
procedure Vector_Test is
   
   X : Sparse_Vector := Sparse ((2.0,  0.0, 0.0));
   Y : Sparse_Vector := Sparse ((0.0, -5.0, 7.0));
   Z : Sparse_Vector;
begin
   Put_Line ("X");
   X.Print;
   Z := X * 11.0;
   Put_Line ("11 * X");
   Z.Print;
   Put_Line ("Y");
   Y.Print;
   Put_Line ("-Y");
   Y := -Y;
   Y.Print;
   Y := -Y;
   Put_Line ("X + Y");
   Z := X + Y;
   Z.Print;
   Put_Line ("Y + X");
   Z := Y + X;
   Z.Print;
   Put_Line ("X - Y");
   Z := X - Y;
   Z.Print;
   Put_Line ("Y - X");
   Z := Y - X;
   Z.Print;
   
end Vector_Test;
