separate (Sparse_Package)

function Minus (Left  : in Matrix;
	      Right : in Matrix) return Matrix is
   Result : Matrix := Right;
begin
   for X of Result.X loop
      X := -X;
   end loop;
   return Left + Result;
end Minus;
