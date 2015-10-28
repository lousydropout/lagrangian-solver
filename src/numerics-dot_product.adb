separate (Numerics)

function Dot_Product (Left_I, Right_J : in Int_Array;
		      Left_X, Right_Y : in Real_Array) return Real is
   Result : Real     := 0.0;
   I      : Int  := Left_I'First;
   J      : Int  := Right_J'First;
begin
   while I <= Left_I'Last and J <= Right_J'Last loop
      if Left_I (I) = Right_J (J) then 
	 Result := Result + Left_X (I) * Right_Y (J);
	 I := I + 1;
      end if;
      
      if I <= Left_I'Last then
	 while J <= Right_J'Last and then Right_J (J) < Left_I (I) loop
	    J := J + 1;
	 end loop;
      end if;
      
      if J <= Right_J'Last then
	 while I <= Left_I'Last and then Left_I (I) < Right_J (J) loop
	    I := I + 1;
	 end loop;
      end if;
   end loop;
   return Result;
end Dot_Product;
