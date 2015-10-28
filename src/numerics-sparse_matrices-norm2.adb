separate (Numerics.Sparse_Matrices)

function Norm2 (Item : in Sparse_Matrix) return Real is
   Sum, Result : Real := 0.0;
   X : Real_Vector renames Item.X;
   P : Int_Vector renames Item.P;
   --  use Real_Functions;
begin
   
   if Item.N_Col = 1 or else Item.N_Row = 1 then
      -- 2-norm for vectors
      for Item of X loop
	 Result := Result + Item ** 2;
      end loop;
      return (Result);
   else
      -- 1-norm for matrices
      for J in 1 .. Item.N_Col loop
	 Sum := 0.0;
	 for I in Item.P (J) .. Item.P (J + 1) - 1 loop
	    Sum := abs (X (I));
	 end loop;
	 Result := Real'Max (Result, Sum);
      end loop;
      return Result ** 2;
   end if;
end Norm2;
