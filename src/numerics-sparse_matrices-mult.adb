separate (Numerics.Sparse_Matrices)

function Mult (Left, Right : in Sparse_Matrix) return Sparse_Matrix is
   Result : Sparse_Matrix;
   Tmp    : constant Sparse_Matrix := Convert (Left);
   
   Left_I  : constant Int_Array  := To_Array (Tmp.I);
   Right_J : constant Int_Array  := To_Array (Right.I);
   Left_X  : constant Real_Array := To_Array (Tmp.X);
   Right_Y : constant Real_Array := To_Array (Right.X);
   
   P : Int_Array (1 .. Nat (Right.P.Length)) := (others => 0);
   
   Sum    : Pos;
   L1, L2 : Pos;
   R1, R2 : Pos;
   use Ada.Containers;
begin
   pragma Assert (Left.Format = CSC);
   pragma Assert (Right.Format = CSC);
   
   Result.Format := CSC;
   Result.N_Row  := Left.N_Row;
   Result.N_Col  := Right.N_Col;
   
   -----------------------------------------------------------------------
   Result.X.Reserve_Capacity (Left.X.Length + Right.X.Length);
   Result.I.Reserve_Capacity (Left.X.Length + Right.X.Length);
   Result.P.Set_Length (Right.P.Length);
   
   for I in 1 .. Nat (Right.P.Length - 1) loop
      R1 := Right.P (I);   
      R2 := Right.P (I + 1) - 1;
      
      Sum := 0;
      for J in 1 .. Nat (Tmp.P.Length - 1) loop
	 L1 := Tmp.P (J);     
	 L2 := Tmp.P (J + 1) - 1;
	 
	 
	 if Left_I (L1 .. L2) * Right_J (R1 .. R2) then
	    Sum := Sum + 1;
	    Result.I.Append (J);
	    Result.X.Append (Dot_Product (Left_I (L1 .. L2), Right_J (R1 .. R2),
					  Left_X (L1 .. L2), Right_Y (R1 .. R2)));
	 end if;
      end loop;
      P (I) := Sum;
   end loop;
   P := Cumulative_Sum (P);
   for I in P'Range loop 
      Result.P (I) := P (I); 
   end loop;
   Result.I.Reserve_Capacity (Result.I.Length);
   Result.X.Reserve_Capacity (Result.I.Length);
   
   return Result;
end Mult;
