separate (Sparse_Package)

function Max_Real_Array (Item : in Real_Array) return Real is
   Result : Real := Item (Item'First);
begin
   for N of Item loop
      Result := Real'Max (Result, N);
   end loop;
   return Result;
end Max_Real_Array;
