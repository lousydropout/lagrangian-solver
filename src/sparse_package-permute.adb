separate (Sparse_Package)

function Permute (Mat : in Matrix;
		  P   : in Int_Array;
		  By  : in Permute_By_Type := Column) return Matrix is
begin
   pragma Assert (Mat.Format = CSC);
   
   case By is
      when Column => 
	 return Permute_By_Col (Mat, P);
      when Row =>
	 return Transpose (Permute_By_Col (Transpose (Mat), P));
   end case;
end Permute;
