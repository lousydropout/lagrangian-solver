separate (Numerics.Sparse_Matrices)

function Plus (Left  : in Sparse_Matrix;
	       Right : in Sparse_Matrix) return Sparse_Matrix is
   type Bool_Array is array (Nat range <>) of Boolean with
     Pack => True;
   N_Col  : Nat renames Left.N_Col;
   N_Row  : Nat renames Left.N_Row;
   Result : Sparse_Matrix;
   Sum    : Pos;
   Value  : Real_Array (1 .. N_Row);
   Exist  : Bool_Array (1 .. N_Row);
   P      : Int_Array  (1 .. N_Col + 1) := (others => 0);

   use Ada.Containers;
begin
   
   Result.Format := CSC;
   Result.N_Col  := N_Col;
   Result.N_Row  := N_Row;
   
   Sum := 0;
   for I in 1 .. Nat (N_Col) loop
      Exist := (others => False);
      
      for J in Left.P (I) .. Left.P (I + 1) - 1 loop
	 Exist (Left.I (J)) := True;
      end loop;
      
      for J in Right.P (I) .. Right.P (I + 1) - 1 loop
	 Exist (Right.I (J)) := True;
      end loop;
      
      -- Sum over Exist
      for E of Exist loop
	 P (I) := P (I) + 1;
      end loop;
      Sum := Sum + P (I);
   end loop;
   P := Cumulative_Sum (P);
   
   Set (Result.P, P);
   Result.X.Reserve_Capacity (Count_Type (Sum));
   Result.I.Reserve_Capacity (Count_Type (Sum));
   
   
   for I in 1 .. Nat (N_Col) loop
      Value := (others => 0.0);
      Exist := (others => False);

      for J in Left.P (I) .. Left.P (I + 1) - 1 loop
	 Exist (Left.I (J)) := True;
	 Value (Left.I (J)) := Left.X (J);
      end loop;
   
      for J in Right.P (I) .. Right.P (I + 1) - 1 loop
	 Exist (Right.I (J)) := True;
	 Value (Right.I (J)) := Value (Right.I (J)) + Right.X (J);
      end loop;
      
      for K in Value'Range loop
	 if Exist (K) then
	    Result.I.Append (K);
	    Result.X.Append (Value (K));
	 end if;
      end loop;
   end loop;
   
   return Result;
end Plus;
