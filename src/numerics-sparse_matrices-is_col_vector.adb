separate (Numerics.Sparse_Matrices)

function Is_Col_Vector (A : in Sparse_Matrix) return Boolean is
begin
   if A.N_Col = 1 then 
      return True;
   else 
      return False; 
   end if;
end Is_Col_Vector;