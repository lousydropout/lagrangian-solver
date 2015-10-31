separate (Numerics.Sparse_Matrices)

--  function Cumulative_Sum (Item : in Int_Array) return Int_Array is
--     Result : Int_Array (Item'Range);
--     Tmp    : Int := 1;
--  begin
--     for I in Item'Range loop
--        Result (I) := Tmp;
--        Tmp := Tmp + Item (I);
--     end loop;
--     return Result;
--  end Cumulative_Sum;

procedure Cumulative_Sum (Item : in out Int_Array) is
   use Ada.Text_IO;
   N : Pos := 1;
   M : Pos;
begin
   for I in Item'Range loop
      M := Item (I);
      Item (I) := N;
      N := N + M;
   end loop;
end Cumulative_Sum;   
