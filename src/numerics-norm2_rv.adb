separate (Numerics)

function Norm2_RV (X : in Real_Vector) return Real is
   Result : Real := 0.0;
begin
   for Item of X loop
      Result := Result + Item ** 2;
   end loop;
   return Result;
end Norm2_RV;
