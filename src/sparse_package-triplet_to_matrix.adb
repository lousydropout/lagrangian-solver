separate (Sparse_Package)

function Triplet_To_Matrix (I      : Int_Array;
			    J      : Int_Array;
			    X      : Real_Array;
			    N_Row  : Pos := 0;
			    N_Col  : Pos := 0;
			    Format : Matrix_Format := CSC) return Matrix is
   Result : Matrix;
begin
   if N_Row = 0 then
      Result.N_Row  := Max (I);
   else
      Result.N_Row  := N_Row;
   end if;
   if N_Col = 0 then
      Result.N_Col  := Max (J);
   else
      Result.N_Col  := N_Col;
   end if;

   Result.Format := Triplet;
   Result.I      := Vectorize (I);
   Result.P      := Vectorize (J);
   Result.X      := Vectorize (X);
   
   case Format is
      when CSC => Compress (Result);
      when Triplet => null;
      when CSR => 
	 Compress (Result);
	 Result.Convert;
   end case;
   return Result;
end Triplet_To_Matrix;
