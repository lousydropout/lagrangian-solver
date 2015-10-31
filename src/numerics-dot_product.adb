separate (Numerics)

function Dot_Product (Left_I, Right_J : in Int_Array;
		      Left_X, Right_Y : in Real_Array) return Real is
   Result : Real     := 0.0;
   I      : Nat  := Left_I'First;
   J      : Nat  := Right_J'First;
   L      : constant Nat := Left_I'Last;
   R      : constant Nat := Right_J'Last;
begin
   
   while I <= L and J <= R loop
      if Left_I (I) = Right_J (J) then 
	 
	 Result := Result + Left_X (I) * Right_Y (J);
	 I := I + 1; J := J + 1;
	 
      elsif I <= L then
	 
	 while J <= R and then Right_J (J) < Left_I (I) loop
	    J := J + 1;
	 end loop;
	 
      elsif J <= R then
	 
	 while I <= L and then Left_I (I) < Right_J (J) loop
	    I := I + 1;
	 end loop;
      end if;
   end loop;

   return Result;
end Dot_Product;
