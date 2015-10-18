with Ada.Text_IO; use Ada.Text_IO;
with Sparse_Package; use Sparse_Package; 

procedure Sparse_Test is
   N : Int := 2;
   
   I1 : Int_Array  := (3,   2,   1);
   J1 : Int_Array  := (1,   3,   2);
   X1 : Real_Array := (1.0, 2.0, 2.0);
   Left : Matrix := Triplet_To_Matrix (I1, J1, X1, 3, 3);
   
   
   I2 : Int_Array  := (1,   2, 2);
   J2 : Int_Array  := (1,   2, 1);
   X2 : Real_Array := (1.5, 3.0, 1.0);
   Right : Matrix 
     := Triplet_To_Matrix (I => I2, 
			   J => J2, 
			   X => X2, 
			   N_Row => N, 
			   N_Col => N + 1);
   
   I3 : Int_Array  := (5, 6, 6, 1, 2, 2, 3, 4, 4);
   J3 : Int_Array  := (1, 1, 2, 4, 4, 5, 7, 7, 8);
   X3 : Real_Array := (1.5, 1.0, 3.0, 3.0, 2.0, 6.0, 3.0, 2.0, 6.0);
   Right2 : Matrix := Triplet_To_Matrix (I3, J3, X3, 6, 9);
   
   X : Real_Array (1 .. N) := (others => 0.0);
   Result : Matrix;
   Vec, X0 : Real_Vector;
   
begin
   
   Left.Print;
end Sparse_Test;
