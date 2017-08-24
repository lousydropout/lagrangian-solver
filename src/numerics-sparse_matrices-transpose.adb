separate (Numerics.Sparse_Matrices)

function Transpose (Mat : in Sparse_Matrix) return Sparse_Matrix is
   Result : Sparse_Matrix := Mat;
begin
   Transposed (Result);
   return Result;
end Transpose;
