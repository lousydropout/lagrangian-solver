separate (Numerics.Sparse_Matrices)

function Kronecker (A, B : in Sparse_Matrix) return Sparse_Matrix is
   Result : Sparse_Matrix;
   N_Col : constant Nat := B.N_Col;
   N_Row : constant Nat := B.N_Row;
   Bl    : Int_Array (1 .. N_Col);
   Al    : Int;
   Tmp   : Nat := 1;
   use Ada.Containers;
begin
   Result.Format := CSC;
   Result.N_Row := A.N_Row * B.N_Row;
   Result.N_Col := A.N_Col * B.N_Col;
   Result.I.Reserve_Capacity (A.I.Length * B.I.Length);
   Result.X.Reserve_Capacity (A.X.Length * B.X.Length);
   Result.P.Reserve_Capacity (Count_Type (Result.N_Col + 1));
   
   -- Set P
   for I in Bl'Range loop
      Bl (I) := B.P (I + 1) - B.P (I); -- #Elements in col(I) of right
   end loop;
   
   Result.P.Append (Tmp); -- Tmp = 1
   for I in 1 .. A.N_Col loop
      Al := A.P (I + 1) - A.P (I); -- #Elements in col(I) of left
      for J in Bl'Range loop       -- Assign Result.P
	 Tmp := Tmp + Pos (Al * Bl (J));
	 Result.P.Append (Tmp);
      end loop;
   end loop;
   
   -- Set I
   for PI in 1 .. A.N_Col loop
      for PJ in 1 .. B.N_Col loop
	 for I in A.P (PI) .. A.P (PI + 1) - 1 loop
	    for J in B.P (PJ) .. B.P (PJ + 1) - 1 loop
	       Result.I.Append ((A.I (I) - 1) * N_Row + B.I (J));
	       Result.X.Append (A.X (I) * B.X (J));
	    end loop;
	 end loop;
      end loop;
   end loop;

   return Result;
end Kronecker;
