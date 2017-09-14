separate (Numerics.Sparse_Matrices)

function Plus2 (Left  : in Sparse_Matrix;
		Right : in Sparse_Matrix) return Sparse_Matrix is
   use Ada.Text_IO, Ada.Containers;
   X, Y : RVector;
   Ai, Aj, Bi, Bj : IVector;
   N_Row, N_Col : Pos;
   N : constant Count_Type := Left.X.Length + Right.X.Length;
   C : Sparse_Matrix;
begin
   pragma Assert (Left.N_Row = Right.N_Row);
   pragma Assert (Left.N_Col = Right.N_Col);
   
   To_Triplet (Left, Ai, Aj, X, N_Row, N_Col);
   To_Triplet (Right, Bi, Bj, Y, N_Row, N_Col);
   
   X.Reserve_Capacity (N);
   Ai.Reserve_Capacity (N);
   Aj.Reserve_Capacity (N);
   
   for I in 1 .. Pos (Y.Length) loop
      Ai.Append (Bi (I));
      Aj.Append (Bj (I));
      X.Append  (Y (I));
   end loop;
   
   Triplet_To_Matrix (C, Ai, Aj, X, N_Row, N_Col);
   return C;
end Plus2;
