separate (Numerics)

function Minus_RV_RV (Left, Right : in Real_Vector) return Real_Vector is
   Result : Real_Vector := RV_Package.To_Vector (Left.Length);
begin
   
   for I in 1 .. Nat (Left.Length) loop
      Result (I) := Left (I) - Right (I);
   end loop;
   
   return Result;
end Minus_RV_RV;      
