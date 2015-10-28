separate (Numerics)

function Abs_Max_RA (Item : in Real_Array) return Real is
   Result : Real := 0.0;
begin
   for N of Item loop
      Result := Real'Max (Result, abs (N));
   end loop;
   return Result;
end Abs_Max_RA;

