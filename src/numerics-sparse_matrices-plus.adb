separate (Numerics.Sparse_Matrices)

function Plus (Left  : in Sparse_Matrix;
	      Right : in Sparse_Matrix) return Sparse_Matrix is
   N_Col  : Nat renames Left.N_Col;
   N_Row  : Nat renames Left.N_Row;
   Result : Sparse_Matrix;
   Sum    : Pos;
   Value  : Real_Array (1 .. N_Row);
   Exist  : Int_Array  (1 .. N_Row);
   P      : Int_Array  (1 .. N_Col + 1) := (others => 0);
   
   L1, L2, R1, R2 : Nat;
   use Ada.Containers;
begin
   
   Result.Format := CSC;
   Result.N_Col  := N_Col;
   Result.N_Row  := N_Row;
   Result.P.Set_Length (Count_Type (N_Col + 1));
   
   Sum := 0;
   for I in 1 .. Nat (N_Col) loop
      Exist := (others => 0);
      
      if I < Nat (Left.P.Length) then
	 L1 := Left.P (I); 
	 L2 := Left.P (I + 1) - 1;
	 if L1 <= L2 then
	    for J in L1 .. L2 loop
	       Exist (Left.I (J)) := 1;
	    end loop;
	 end if;
      end if;
      
      if I < Nat (Right.P.Length) then
	 R1 := Right.P (I); 
	 R2 := Right.P (I + 1) - 1;
	 if R1 <= R2 then 
	    for J in R1 .. R2 loop
	       Exist (Right.I (J)) := 1;
	    end loop;
	 end if;
      end if;
      
      -- Sum over Exist
      for E of Exist loop
	 P (I) := P (I) + E;
      end loop;
      Sum := Sum + P (I);
   end loop;
   P := Cumulative_Sum (P);
   
   Set (Result.P, P);
   Result.X.Set_Length (Count_Type (Sum));
   Result.I.Set_Length (Count_Type (Sum));
   
   
   Sum := 0;      
   for I in 1 .. Nat (N_Col) loop
      Value := (others => 0.0);
      Exist := (others => 0);
      if I < Nat (Left.P.Length) then
	 L1 := Left.P (I); 
	 L2 := Left.P (I + 1) - 1;
	 if L1 <= L2 then
	    for J in L1 .. L2 loop
	       Exist (Left.I (J)) := 1;
	       Value (Left.I (J)) := Left.X (J);
	    end loop;
	 end if;
      end if;
      
      if I < Nat (Right.P.Length) then
	 R1 := Right.P (I); 
	 R2 := Right.P (I + 1) - 1;
	 if R1 <= R2 then 
	    for J in R1 .. R2 loop
	       Exist (Right.I (J)) := 1;
	       Value (Right.I (J)) := Value (Right.I (J)) + Right.X (J);
	    end loop;
	 end if;
      end if;
      
      for K in Value'Range loop
	 if Exist (K) = 1 then
	    Sum := Sum + 1;
	    Result.I (Sum) := K;
	    Result.X (Sum) := Value (K);
	 end if;
      end loop;
   end loop;
   
   return Result;
end Plus;
