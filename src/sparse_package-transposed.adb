separate (Sparse_Package)

procedure Transposed (Mat : in out Matrix) is
   Tmp : Pos;
begin
   
   Mat.Convert;
   
   case Mat.Format is
      when Triplet => 
	 null;
      when CSR => 
	 Mat.Format := CSC;
      when CSC => 
	 Mat.Format := CSR;
   end case;
   
   Tmp       := Mat.N_Row;
   Mat.N_Row := Mat.N_Col;
   Mat.N_Col := Tmp;
   
end Transposed;

