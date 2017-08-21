separate (Numerics.Sparse_Matrices)

function Omega (N : in Nat;
		M : in Pos := 0) return Sparse_Matrix is
   Result : Sparse_Matrix;
   use Ada.Containers;
begin
   
   -----------------------------------------------
   -- omega (N, M) = [ 0_N   -1_N      ]
   --                [ 1_N    0_N      ]
   --                [             0_M ]
   -- and is such that for
   --     X = [q, p, lambda]
   -- where 
   --     q, p \in R^N 
   -- and
   --     lambda \in R^M,
   -- the Hamiltonian equations are
   --     omega (N, M) (d/dt) X = grad H
   ----------------------------------------------
   Result.Format := CSC;
   Result.N_Row  := 2 * N + M;
   Result.N_Col  := 2 * N + M;

   Result.P.Reserve_Capacity (Count_Type (2 * N + M + 1));
   Result.I.Reserve_Capacity (Count_Type (2 * N));
   Result.X.Reserve_Capacity (Count_Type (2 * N));
   
   for I in 1 .. N loop
      Result.P.Append (I);
      Result.I.Append (N + I);
      Result.X.Append (1.0);
   end loop;
   
   for I in 1 .. N loop
      Result.P.Append (N + I);
      Result.I.Append (I);
      Result.X.Append (-1.0);
   end loop;
   
   for I in 2 * N + 1 .. 2 * N + M + 1 loop
      Result.P.Append (I);
   end loop;
   
   return Result;
end Omega;
