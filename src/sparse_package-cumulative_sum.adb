separate (Sparse_Package)

function Cumulative_Sum (Item : in Int_Array) return Int_Array is
   Result : Int_Array (Item'Range);
   Tmp    : Int := 1;
begin
   for I in Item'Range loop
      Result (I) := Tmp;
      Tmp := Tmp + Item (I);
   end loop;
   return Result;
end Cumulative_Sum;
