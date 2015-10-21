separate (Sparse_Package)

procedure Matrix_To_Sparse (Mat    : in     Matrix;
			    Sparse :    out Sparse_Ptr) is
begin
   Sparse := To_CS (M     => Mat.N_Row, 
		    N     => Mat.N_Col, 
		    Nzmax => Int (Mat.X.Length), 
		    I     => To_Array (Mat.I),
		    P     => To_Array (Mat.P),
		    X     => To_Array (Mat.X));
end Matrix_To_Sparse;
