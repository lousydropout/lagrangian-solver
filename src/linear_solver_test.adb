with Sparse_Package; use Sparse_Package;

procedure Linear_Solver_Test is

   I1 : Int_Array  := (1, 2);
   J1 : Int_Array  := (1, 2);
   X1 : Real_Array := (1.234, 2.345);
   Mat : Matrix := Triplet_To_Matrix (I1, J1, X1, 2, 2);
   
   LU : LU_Type := LU_Decomposition (Mat);
   Res : Real_Vector := Vectorize ((1.0, 2.0));
begin

   Mat.Print;
   New_Line;
   Res := Solve (LU, Res);
   Res := Mat * Res;
   
   for X of Res loop
      Real_IO.Put (X); New_Line;
   end loop;
   null;
end Linear_Solver_Test;
