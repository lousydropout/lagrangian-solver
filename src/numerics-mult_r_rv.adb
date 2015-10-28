separate (Numerics)


function Mult_R_RV (Left  : in Real;
		    Right : in Real_Vector) return Real_Vector is
   Result : Real_Vector;
begin
   Result.Set_Length (Right.Length);
   for I in 1 .. Nat (Right.Length) loop
      Result (I) := Left * Right (I);
   end loop;
   return Result;
end Mult_R_RV;

