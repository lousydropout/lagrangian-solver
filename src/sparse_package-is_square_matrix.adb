separate (Sparse_Package)

function Is_Square_Matrix (A : in Sparse_Matrix) return Boolean is
begin
   if A.N_Row = A.N_Col then
      return True;
   else
      return False;
   end if;
end Is_Square_Matrix;
