separate (Numerics.Sparse_Matrices)

function Eye (N : in Nat) return Sparse_Matrix is
   Result : Sparse_Matrix;
   use Ada.Containers;
begin
   Result.Format := CSC;
   Result.X := RV_Package.To_Vector (1.0, Count_Type (N));
   Result.I.Reserve_Capacity (Count_Type (N));
   Result.P.Reserve_Capacity (Count_Type (N + 1));
   
   for I in 1 .. N loop
      Result.I.Append (I);
      Result.P.Append (I);
   end loop;
   Result.P.Append (N + 1);
   return Result;
end Eye;
