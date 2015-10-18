separate (Sparse_Package)

function Zero_Vector (N : in Nat) return Matrix is
   X : Real_Array := (1 => 0.0);
   J : Int_Array  := (1 => 1);
   Result : Matrix := Vectorize (J, X);
begin
   Result.N_Row := N;
   return Result;
end Zero_Vector;
