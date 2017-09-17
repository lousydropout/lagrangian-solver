separate (Numerics.Sparse_Matrices)

--  function Minus (Left  : in Sparse_Matrix;
--  		Right : in Sparse_Matrix) return Sparse_Matrix is
--     Result : Sparse_Matrix := Right;
--  begin
--     for X of Result.X loop
--        X := -X;
--     end loop;
--     return Left + Result;
--  end Minus;

function Minus (Left  : in Sparse_Matrix;
	       Right : in Sparse_Matrix) return Sparse_Matrix is
   use Ada.Containers, Ada.Text_IO;
   A : Sparse_Matrix renames Left;
   B : Sparse_Matrix renames Right;
   C : Sparse_Matrix;
   N_Col  : Nat renames A.N_Col;
   N_Row  : constant Count_Type := Count_Type (A.N_Row);
   Res    : constant Count_Type := Count_Type'Max (A.X.Length, B.X.Length);
   Sum    : Pos := 1;
   P      : Pos;
   N, M   : Pos;
   L, R   : Pos;
   AI, BI : Nat;
begin
   if A.X.Length = 0 then return B;
   elsif B.X.Length = 0 then return A;
   end if;
   C.Format := CSC; C.N_Col := N_Col; C.N_Row := Pos (N_Row);
   
   C.P.Reserve_Capacity (Count_Type (N_Col + 1));
   C.X.Reserve_Capacity (Res);
   C.I.Reserve_Capacity (Res);
   
   C.P.Append (1); -- 1st element is always 1
   for I in 1 .. Nat (N_Col) loop
      if C.X.Capacity < C.X.Length + N_Row then
      	 C.X.Reserve_Capacity (C.X.Capacity + N_Row);
      	 C.I.Reserve_Capacity (C.X.Capacity + N_Row);
      end if;

      N := A.P (I);  L := A.P (I + 1) - 1;
      M := B.P (I);  R := B.P (I + 1) - 1;
      P := 0;
      while N <= L and M <= R loop
	 P := P + 1;
	 AI := A.I (N); BI := B.I (M);
	 if AI = BI then
	    C.X.Append (A.X (N) - B.X (M));
	    C.I.Append (AI);
	    N := N + 1; M := M + 1;
	 elsif AI < BI then
	    C.X.Append (A.X (N));
	    C.I.Append (AI);
	    N := N + 1;
	 elsif AI > BI then
	    C.X.Append (-B.X (M));
	    C.I.Append (BI);
	    M := M + 1;
	 end if;
      end loop;
      while M > R and then N <= L loop
	 P := P + 1;
	 C.X.Append (A.X (N));
	 C.I.Append (A.I (N));
	 N := N + 1;
      end loop;
      while N > L and then M <= R loop
	 P := P + 1;
	 C.X.Append (-B.X (M));
	 C.I.Append  (B.I (M));
	 M := M + 1;
      end loop;
      Sum := Sum + P;
      C.P.Append (Sum);
   end loop;
   Sum := Sum - 1;
   
   C.X.Reserve_Capacity (Count_Type (Sum));
   C.I.Reserve_Capacity (Count_Type (Sum));
   
   return C;
end Minus;
