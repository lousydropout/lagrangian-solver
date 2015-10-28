separate (Numerics)

function Add_RV_RV (Left, Right : in Real_Vector) return Real_Vector is 
   Result : Real_Vector;
begin
   pragma Assert (Nat (Left.Length) = Nat (Right.Length));
   Result.Set_Length (Left.Length);
   for I in 1 .. Nat (Left.Length) loop
      Result (I) := Left (I) + Right (I);
   end loop;
   return Result;
end Add_RV_RV;
