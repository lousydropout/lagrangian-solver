separate (Numerics.Sparse_Matrices)

function Triplet_To_Matrix (I      : in Int_Array;
			    J      : in Int_Array;
			    X      : in Real_Array;
			    N_Row  : in Pos := 0;
			    N_Col  : in Pos := 0;
			    Format : in Sparse_Matrix_Format := CSC) return Sparse_Matrix is
   Result : Sparse_Matrix;
begin
   Result.N_Row  := Pos'Max (N_Row, Max (I));
   Result.N_Col  := Pos'Max (N_Col, Max (J));
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
