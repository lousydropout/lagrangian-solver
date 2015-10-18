separate (Sparse_Package)

function Eye (N : in Nat) return Matrix is
   Result : Matrix;
   use Ada.Containers;
begin
   Result.Format := CSC;
   Result.I.Set_Length (Count_Type (N));
   Result.X.Set_Length (Count_Type (N));
   Result.P.Set_Length (Count_Type (N + 1));
   
   for I in 1 .. N loop
      Result.I (I) := I;
      Result.P (I) := I;
      Result.X (I) := 1.0;
   end loop;
   Result.P (N + 1) := N + 1;
   return Result;
end Eye;
