separate (Numerics)

function Dot_Product_RV (X, Y : in Real_Vector) return Real is
   Result : Real := 0.0;
begin
   pragma Assert (Nat (X.Length) = Nat (Y.Length));
   for I in 1 .. Nat (X.Length) loop
      Result := Result + X (I) * Y (I);
   end loop;
   return Result;
end Dot_Product_RV;
