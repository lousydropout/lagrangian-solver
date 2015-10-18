separate (Sparse_Package)

function Mult_Int_Array (Left, Right : in Int_Array) return Boolean is
   I : Int := Left'First;
   J : Int := Right'First;
begin
   while I <= Left'Last and J <= Right'Last loop
      if Left (I) = Right (J) then  return True; end if;
      
      if I <= Left'Last then
	 while J <= Right'Last and then Right (J) < Left (I) loop
	    J := J + 1;
	 end loop;
      end if;
      
      if J <= Right'Last then
	 while I <= Left'Last and then Left (I) < Right (J) loop
	    I := I + 1;
	 end loop;
      end if;
   end loop;
   return False;
end Mult_Int_Array;
