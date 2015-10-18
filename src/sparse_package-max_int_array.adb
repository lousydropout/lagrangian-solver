separate (Sparse_Package)

function Max_Int_Array (Item : in Int_Array) return Int is
   Result : Int := Item (Item'First);
begin
   for N of Item loop
      Result := Int'Max (Result, N);
   end loop;
   return Result;
end Max_Int_Array;
