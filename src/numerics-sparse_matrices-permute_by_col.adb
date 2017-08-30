separate (Numerics.Sparse_Matrices)

function Permute_By_Col (Mat : in Sparse_Matrix;
			 P   : in Int_Array) return Sparse_Matrix is
   Result : Sparse_Matrix;
   Tmp    : Nat := 1;
   use Ada.Text_IO;
begin
   pragma Assert (Mat.Format = CSC);
   pragma Assert (P'Length = Integer (Mat.P.Length) - 1);
   Result.Format := CSC;
   Result.N_Row  := Mat.N_Row;
   Result.N_Col  := Mat.N_Col;
   
   Result.X.Reserve_Capacity (Mat.X.Length);
   Result.I.Reserve_Capacity (Mat.I.Length);
   Result.P.Reserve_Capacity (Mat.P.Length);
   
   for J of P loop
      for K in Mat.P (J) .. Mat.P (J + 1) - 1 loop
	 Result.X.Append (Mat.X (K));
	 Result.I.Append (Mat.I (K));
      end loop;
      Result.P.Append (Tmp);
      Tmp := Tmp + Mat.P (J + 1) - Mat.P (J);
   end loop;
   Result.P.Append (Tmp);
   
   return Result;
end Permute_By_Col;
