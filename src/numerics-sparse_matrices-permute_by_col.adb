separate (Numerics.Sparse_Matrices)

function Permute_By_Col (Mat : in Sparse_Matrix;
			 P   : in Int_Array) return Sparse_Matrix is
   Result : Sparse_Matrix;
   Ind    : Pos := 1;
   Tmp    : Nat := 1;
   use Ada.Text_IO;
begin
   pragma Assert (Mat.Format = CSC);
   pragma Assert (P'Length = Int (Mat.P.Length) - 1);
   Result.Format := CSC;
   Result.N_Row  := Mat.N_Row;
   Result.N_Col  := Mat.N_Col;
   
   Result.X.Set_Length (Mat.X.Length);
   Result.I.Set_Length (Mat.I.Length);
   Result.P.Reserve_Capacity (Mat.P.Length);
   pragma Assert (P'Length = Int (Mat.P.Length) - 1);
   
   for J of P loop
      Result.P.Append (Tmp);
      Tmp := Tmp + Mat.P (J + 1) - Mat.P (J);
      for K in Mat.P (J) .. Mat.P (J + 1) - 1 loop
	 Result.I (Ind) := Mat.I (K);
	 Result.X (Ind) := Mat.X (K);
	 Ind            := Ind + 1;
      end loop;
   end loop;
   Result.P.Append (Tmp);
   
   return Result;
end Permute_By_Col;
