with Ada.Text_IO, Sparse_Package;
use  Ada.Text_IO, Sparse_Package;

procedure Linear_Solver_Test is

   N : Int := 2;
   
   I1 : Int_Array  := (1, 2);
   J1 : Int_Array  := (1, 2);
   X1 : Real_Array := (1.0, 2.0); -- 1.234, 2.345);
   Left : Matrix := Triplet_To_Matrix (I1, J1, X1, 2, 2);
   
   LU : LU_Type := LU_Decomposition (Left);
   Res : Real_Array := (1.0, 2.0);
begin
   Left.Print;
   New_Line;
   Res := Solve (LU, Res);
   for X of Res loop
      Real_IO.Put (X); New_Line;
   end loop;
   null;
end Linear_Solver_Test;
