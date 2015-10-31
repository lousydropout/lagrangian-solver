separate (Numerics.Sparse_Matrices)

--  function Kronecker (A, B : in Sparse_Matrix) return Sparse_Matrix is
--     Result : Sparse_Matrix;
--     N_Col : Nat := B.N_Col;
--     N_Row : Nat := B.N_Row;
--     Al : Int_Array (1 .. Nat (A.P.Length) - 1);
--     Bl : Int_Array (1 .. Nat (B.P.Length) - 1);
--     Tmp : Nat := 1;
--     use Ada.Containers;
--  begin
--     Result.Format := CSC;
--     Result.N_Row := A.N_Row * B.N_Row;
--     Result.N_Col := A.N_Col * B.N_Col;
--     Result.I.Set_Length (A.I.Length * B.I.Length);
--     Result.X.Set_Length (A.X.Length * B.X.Length);
--     Result.P.Set_Length (Count_Type (Result.N_Col) + 1);
   
   
--     -- Set P
--     for I in Bl'Range loop
--        Bl (I) := B.P (I + 1) - B.P (I); -- #Elements in col(I) of right
--     end loop;
   
--     Result.P (1) := 1;
--     for I in Al'Range loop
--        Al (I) := A.P (I + 1) - A.P (I); -- #Elements in col(I) of left
--  					     -- Assign Result.P
--        for J in Bl'Range loop
--  	 Tmp := Tmp + Pos (Al (I) * Bl (J));
--  	 Result.P ((I - 1) * Bl'Length + J + 1) := Tmp;
--        end loop;
--     end loop;
   
--     -- Set I
--     Tmp := 1;
--     for PI in 1 .. Nat (A.P.Length - 1) loop
--        for PJ in 1 .. Nat (B.P.Length - 1) loop
--  	 for I in A.P (PI) .. A.P (PI + 1) - 1 loop
--  	    for J in B.P (PJ) .. B.P (PJ + 1) - 1 loop
--  	       Result.I (Tmp) := (A.I (I) - 1) * N_Row + B.I (J);
--  	       Tmp := Tmp + 1;
--  	    end loop;
--  	 end loop;
--        end loop;
--     end loop;
--     -- Set up X
--     Tmp := 1;
--     for I of A.X loop
--        for K of B.X loop
--  	 Result.X (Tmp) := I * K;
--  	 Tmp := Tmp + 1;
--        end loop;
--     end loop;
--     return Result;
--  end Kronecker;


function Kronecker (A, B : in Sparse_Matrix) return Sparse_Matrix is
   C : Sparse_Matrix;
   N_Col : constant Nat := B.N_Col;
   N_Row : constant Nat := B.N_Row;
   Al : Int;
   Bl : Int_Array (1 .. B.N_Col);
   Tmp : Nat := 1;
   use Ada.Containers;
begin
   C.Format := CSC;
   C.N_Row := A.N_Row * B.N_Row;
   C.N_Col := A.N_Col * B.N_Col;
   C.I.Reserve_Capacity (A.I.Length * B.I.Length);
   C.X.Reserve_Capacity (A.X.Length * B.X.Length);
   C.P.Reserve_Capacity (Count_Type (C.N_Col + 1));
   
   -- Set X
   for X of A.X loop
      for Y of B.X loop C.X.Append (X * Y); end loop;
   end loop;
   -- Set I
   for I of A.I loop
      for J of B.I loop C.I.Append ((A.I (I) - 1) * N_Row + B.I (J)); end loop;
   end loop;
   
   -- Set P
   for I in Bl'Range loop
      Bl (I) := B.P (I + 1) - B.P (I); -- #Elements in col(I) of right
   end loop;
   
   C.P.Append (1);
   for I in 1 .. A.N_Col loop
      Al := A.P (I + 1) - A.P (I); -- #Elements in col(I) of left
      for J in Bl'Range loop           -- Assign C.P
	 Tmp := Tmp + Al * Bl (J);
	 C.P.Append (Tmp);
      end loop;
   end loop;
   
   return C;
end Kronecker;
