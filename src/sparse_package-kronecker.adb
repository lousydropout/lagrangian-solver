separate (Sparse_Package)

function Kronecker (Left, Right : in Matrix) return Matrix is
   Result : Matrix;
   N_Col : Nat := Right.N_Col;
   N_Row : Nat := Right.N_Row;
   Al : Int_Array (1 .. Nat (Left.P.Length) - 1);
   Bl : Int_Array (1 .. Nat (Right.P.Length) - 1);
   Tmp : Nat := 1;
   use Ada.Containers;
begin
   Result.Format := CSC;
   Result.N_Row := Left.N_Row * Right.N_Row;
   Result.N_Col := Left.N_Col * Right.N_Col;
   Result.I.Set_Length (Left.I.Length * Right.I.Length);
   Result.X.Set_Length (Left.X.Length * Right.X.Length);
   Result.P.Set_Length (Count_Type (Result.N_Col) + 1);
   
   
   -- Set P
   for I in Bl'Range loop
      Bl (I) := Right.P (I + 1) - Right.P (I); -- #Elements in col(I) of right
   end loop;
   
   Result.P (1) := 1;
   for I in Al'Range loop
      Al (I) := Left.P (I + 1) - Left.P (I); -- #Elements in col(I) of left
					     -- Assign Result.P
      for J in Bl'Range loop
	 Tmp := Tmp + Pos (Al (I) * Bl (J));
	 Result.P ((I - 1) * Bl'Length + J + 1) := Tmp;
      end loop;
   end loop;
   
   -- Set I
   Tmp := 1;
   for PI in 1 .. Nat (Left.P.Length - 1) loop
      for PJ in 1 .. Nat (Right.P.Length - 1) loop
	 for I in Left.P (PI) .. Left.P (PI + 1) - 1 loop
	    for J in Right.P (PJ) .. Right.P (PJ + 1) - 1 loop
	       Result.I (Tmp) := (Left.I (I) - 1) * N_Row + Right.I (J);
	       Tmp := Tmp + 1;
	    end loop;
	 end loop;
      end loop;
   end loop;
   -- Set up X
   Tmp := 1;
   for I of Left.X loop
      for K of Right.X loop
	 Result.X (Tmp) := I * K;
	 Tmp := Tmp + 1;
      end loop;
   end loop;
   return Result;
end Kronecker;
