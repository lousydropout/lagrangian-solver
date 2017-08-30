separate (Numerics)

function Abs_Max_IA (Item : in Int_Array) return Integer is
   Result : Integer := 0;
begin
   for N of Item loop
      Result := Integer'Max (Result, abs (N));
   end loop;
   return Result;
end Abs_Max_IA;
