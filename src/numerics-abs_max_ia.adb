separate (Numerics)

function Abs_Max_IA (Item : in Int_Array) return Int is
   Result : Int := 0;
begin
   for N of Item loop
      Result := Int'Max (Result, abs (N));
   end loop;
   return Result;
end Abs_Max_IA;
