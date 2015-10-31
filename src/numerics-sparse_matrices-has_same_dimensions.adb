separate (Numerics.Sparse_Matrices)

function Has_Same_Dimensions (Left, Right : in Sparse_Matrix) return Boolean is
begin
   if Left.N_Row /= Right.N_Row or else Left.N_Col /= Right.N_Col then
      return False;
   end if;
   return True;
end Has_Same_Dimensions;
