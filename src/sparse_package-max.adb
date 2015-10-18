separate (Sparse_Package)

function Max (Item : in Int_Array) return Int is
   Result : Int := Item (Item'First);
begin
   for N of Item loop
      Result := Int'Max (Result, N);
   end loop;
   return Result;
end Max;


function Max (Item : in Real_Array) return Real is
   Result : Real := Item (Item'First);
begin
   for N of Item loop
      Result := Real'Max (Result, N);
   end loop;
   return Result;
end Max;
