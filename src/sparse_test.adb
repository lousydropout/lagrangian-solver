with Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;
use  Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;


with Ada.Text_IO; use Ada.Text_IO;
procedure Sparse_Test is
   N : Int := 2;
   
   I1 : Int_Array  := (2,   1,   1,  2);
   J1 : Int_Array  := (1,   3,   1,  3);
   X1 : Real_Array := (1.0, 2.0, 3.0, 5.0);
   
   A  : Sparse_Matrix := Triplet_To_Matrix (I1, J1, X1, 3, 3);
   B  : Sparse_Matrix := Triplet_To_Matrix (J1, I1, X1, 3, 3);
   X  : Sparse_Vector := Sparse ((0.0, 1.0, 0.0));
   
   C : Sparse_Matrix;
   Y : Real_Array := To_Array (X);
begin
   C := A or B;
   --  A := Transpose (A) + Transpose (B);
   --  A.Print;
   --  A.Add (2, 2, -2.0);
   --  A.Set (2, 2, 3.14);
   --  A.Print;
   Put_Line ("B: ");
   B.Print;
   X := Sparse ((6.789, 999.314, 0.0));
   B.Set_Diag (To => X);
   Put_Line ("B: ");
   B.Print;
   Put_Line ("X:");
   X := B.Diag;
   Print (X);
   B := Diag (X);
   Put_Line ("B: ");
   B.Print;
   New_Line; New_Line;
   --  Put_Line ("solution: ");
   --  A.Add (1, 2, -1.0);
   --  A.Print;
   --  X := A * X;
   --  --  B.Transposed;
   --  X.Print;
   --  Put_Line ("array (X): ");
   --  for I of Y loop
   --     Put_Line (I'Img);
   --  end loop;
   --  C.Print;
   
   --  C.Print;
   --  C := Transpose (B * B);
   --  C.Print;
   --  New_Line;
   --  Put_Line ("Number of Elements = " & Int'Image (Number_Of_Elements (A)));

end Sparse_Test;
