separate (Numerics)

function Max_Int_Array (Item : in Int_Array) return Integer is
   Result : Integer := Item (Item'First);
begin
   for N of Item loop
      Result := Integer'Max (Result, N);
   end loop;
   return Result;
end Max_Int_Array;
