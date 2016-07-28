with Numerics, Numerics.Sparse_Matrices;
use  Numerics, Numerics.Sparse_Matrices;

with Ada.Text_IO; use Ada.Text_IO;
procedure Vector_Test is
   
   X : Sparse_Vector := Sparse ((0.0,  0.0, 0.0));
   Y : Sparse_Vector; -- := Sparse ((0.0, -5.0, 7.0));
   Z : Sparse_Vector;
begin
   --  X.Set_Length (3);
   Print (X);
   Set (X, 1, π);
   Add (X, 2, 9.9999);
   
   Put_Line ("X");
   Print (X);
   
   --  Z := X * 11.0;
   --  Put_Line ("11 * X");
   --  Print (Z);
   
   --  Put_Line ("Y");
   --  Print (Y);
   
   --  Put_Line ("-Y");
   --  Z := -Y;
   --  Print (Z);

   --  Put_Line ("X + Y");
   --  Z := X + Y;
   --  Print (Z);
   
   --  Put_Line ("Y + X");
   --  Z := Y + X;
   --  Print (Z);
   
   --  Put_Line ("X - Y");
   --  Z := X - Y;
   --  Print (Z);

   --  Put_Line ("Y - X");
   --  Z := Y - X;
   --  Print (Z);
   
end Vector_Test;
