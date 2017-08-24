separate (Numerics.Sparse_Matrices)

procedure Transposed (Mat : in out Sparse_Matrix) is
   Tmp : Pos;
begin
   
   Convert (Mat);
   
   case Mat.Format is
      when Triplet => null;
      when CSR => Mat.Format := CSC;
      when CSC => Mat.Format := CSR;
   end case;
   
   Tmp       := Mat.N_Row;
   Mat.N_Row := Mat.N_Col;
   Mat.N_Col := Tmp;
   
end Transposed;

