separate (Numerics.Sparse_Matrices)

function Has_Same_Dimensions (Left, Right : in Sparse_Matrix) return Boolean is
   Result : Boolean := True;
begin
   if Left.N_Row /= Right.N_Row or else Left.N_Col /= Right.N_Col then
      Result := False;
   end if;
   return Result;
end Has_Same_Dimensions;
