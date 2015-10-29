separate (Numerics.Sparse_Matrices)

function Triplet_To_Matrix (I      : in Int_Array;
			    J      : in Int_Array;
			    X      : in Real_Array;
			    N_Row  : in Pos := 0;
			    N_Col  : in Pos := 0;
			    Format : in Sparse_Matrix_Format := CSC) return Sparse_Matrix is
   Result : Sparse_Matrix;
begin
   Result.N_Row  := (if N_Row = 0 then Max (I) else N_Row);
   Result.N_Col  := (if N_Col = 0 then Max (J) else N_Col);
   Result.Format := Triplet;
   
   Set (X => Result.I, To => I);
   Set (X => Result.P, To => J);
   Set (X => Result.X, To => X);
   
   case Format is
      when CSC => Compress (Result);
      when Triplet => null;
      when CSR => 
	 Result.Compress;
	 Result.Convert;
   end case;
   return Result;
end Triplet_To_Matrix;
