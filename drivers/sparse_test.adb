with Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;
use  Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;

with Ada.Text_IO; use Ada.Text_IO;
procedure Sparse_Test is
   
   A  : Sparse_Matrix;
   B  : Sparse_Matrix;
   C  : Sparse_Matrix;
   X  : Sparse_Vector := Sparse ((0.0, 0.0, 0.0, 1.0));
begin
   
   B := Eye (4); Print (B);
   Put_Line ("---------------------------------------------------------");
   A := X * X;   Print (A);
   Put_Line ("---------------------------------------------------------");
   C := A + B;   C.Print;
   Put_Line ("---------------------------------------------------------");

end Sparse_Test;
