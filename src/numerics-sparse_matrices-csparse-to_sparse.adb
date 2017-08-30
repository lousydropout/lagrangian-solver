separate (Numerics.Sparse_Matrices.CSparse)
function To_Sparse (Mat : in Sparse_Matrix) return Sparse_Ptr is
begin
   return To_CS (M     => Cint (N_Row (Mat)), 
		 N     => Cint (N_Col (Mat)), 
		 Nzmax => Cint (Mat.X.Length), 
		 I     => To_Array (Mat.I),
		 P     => To_Array (Mat.P),
		 X     => To_Array (Mat.X));
end To_Sparse;
