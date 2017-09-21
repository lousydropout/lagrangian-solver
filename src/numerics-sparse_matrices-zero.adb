separate (Numerics.Sparse_Matrices)

function Zero (N : in Pos) return Sparse_Matrix is
   Result : Sparse_Matrix;
   use Ada.Containers;
begin
   Result.Format := CSC;
   Result.N_Row  := N;
   Result.N_Col  := N;
   
   Result.X.Reserve_Capacity (0);
   Result.I.Reserve_Capacity (0);
   Result.P.Reserve_Capacity (Count_Type (N + 1));
   
   for I in 1 .. N + 1 loop
      Result.P.Append (1);
   end loop;
   return Result;
end Zero;
