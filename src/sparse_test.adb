with Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;
use  Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;

with Ada.Text_IO; use Ada.Text_IO;
procedure Sparse_Test is
   N : Int := 2;
   
   I1 : Int_Array  := (2,   1,   1,  2);
   J1 : Int_Array  := (1,   3,   1,  3);
   X1 : Real_Array := (1.0, 2.0, 3.0, 5.0);
   
   A  : Sparse_Matrix := Triplet_To_Matrix (I1, J1, X1);
   B  : Sparse_Matrix := Triplet_To_Matrix (J1, I1, X1);
   
   C : Sparse_Matrix;
begin
   C := B and A;
   --  A := Transpose (A) + Transpose (B);
   A.Print;
   --  B.Transposed;
   B.Print;
   C.Print;
   
   --  C.Print;
   --  C := Transpose (B * B);
   --  C.Print;
   New_Line;
   Put_Line ("Number of Elements = " & Int'Image (Number_Of_Elements (C)));

end Sparse_Test;
