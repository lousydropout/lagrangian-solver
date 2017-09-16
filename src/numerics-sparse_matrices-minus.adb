separate (Numerics.Sparse_Matrices)

function Minus (Left  : in Sparse_Matrix;
		Right : in Sparse_Matrix) return Sparse_Matrix is
   Result : Sparse_Matrix := Right;
begin
   for X of Result.X loop
      X := -X;
   end loop;
   return Left + Result;
end Minus;
