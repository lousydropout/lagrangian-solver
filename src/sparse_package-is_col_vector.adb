separate (Sparse_Package)

function Is_Col_Vector (A : in Matrix) return Boolean is
begin
   if A.N_Col = 1 then 
      return True;
   else 
      return False; 
   end if;
end Is_Col_Vector;
