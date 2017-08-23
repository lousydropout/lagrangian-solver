separate (Numerics.Sparse_Matrices)

function Zero (N : in Nat) return Sparse_Matrix is
   Result : Sparse_Matrix;
   use Ada.Containers;
begin
   Result.Format := CSC;
   Result.N_Row  := N;
   Result.N_Col  := N;

   Result.X.Reserve_Capacity (Count_Type (1));
   Result.I.Reserve_Capacity (Count_Type (1));
   Result.P.Reserve_Capacity (Count_Type (N + 1));
   
   Result.I.Append (1);
   Result.X.Append (0.0);
   Result.P.Append (1);
   for I in 2 .. N + 1 loop
      Result.P.Append (2);
   end loop;
   return Result;
end Zero;
