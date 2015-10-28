separate (Numerics.Sparse_Matrices)

function Zero_Vector (N : in Nat) return Sparse_Matrix is
   X : Real_Array := (1 => 0.0);
   J : Int_Array  := (1 => 1);
   Result : Sparse_Matrix := Vectorize (J, X);
begin
   Result.N_Row := N;
   return Result;
end Zero_Vector;
