separate (Sparse_Package)

function Transpose (Mat : in Matrix) return Matrix is
   Result : Matrix := Mat;
begin
   Result.Transposed;
   return Result;
end Transpose;
