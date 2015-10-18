separate (Sparse_Package)

function Norm_RV (X : in Real_Vector) return Real is
   use Real_Functions;
begin
   return Sqrt (Norm2 (X));
end Norm_RV;
