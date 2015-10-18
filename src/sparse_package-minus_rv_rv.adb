separate (Sparse_Package)

function Minus_RV_RV (Left, Right : in Real_Vector) return Real_Vector is
   use RV_Package;
   Result : Real_Vector := To_Vector (Left.Length);
begin
   pragma Assert (Nat (Left.Length) = Nat (Right.Length));
   
   for I in 1 .. Nat (Left.Length) loop
      Result (I) := Left (I) - Right (I);
   end loop;
   return Result;
end Minus_RV_RV;      
