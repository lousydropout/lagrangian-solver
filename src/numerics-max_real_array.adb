separate (Numerics)

function Max_Real_Array (Item : in Real_Vector) return Real is
   Result : Real := Item (Item'First);
begin
   for N of Item loop
      Result := Real'Max (Result, N);
   end loop;
   return Result;
end Max_Real_Array;
