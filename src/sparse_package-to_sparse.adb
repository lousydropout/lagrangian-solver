separate (Sparse_Package)
function To_Sparse (Mat : in Sparse_Matrix) return Sparse_Ptr is
begin
   return To_CS (M     => Mat.N_Row, 
		 N     => Mat.N_Col, 
		 Nzmax => Int (Mat.X.Length), 
		 I     => To_Array (Mat.I),
		 P     => To_Array (Mat.P),
		 X     => To_Array (Mat.X));
end To_Sparse;
