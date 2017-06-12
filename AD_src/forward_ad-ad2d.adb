package body Forward_AD.AD2D is
   
   function "+" (A, B : in AD2D) return AD2D is
   begin
      return (X => A.X + B.X,
	      Y => A.Y + B.Y);
   end "+";
   
   function "-" (A, B : in AD2D) return AD2D is
   begin
      return (X => A.X - B.X,
	      Y => A.Y - B.Y);
   end "-";

   function "*" (A : in AD2D;
		 B : in AD2D) return AD_Type is
      Result : AD_Type := A.X * B.X + A.Y * B.Y;
   begin
      return Result;
   end "*";
   
   function "-" (A : in AD2D) return AD2D is
   begin
      return (X => -A.X, Y => -A.Y);
   end "-";

   
   
   
   function To_AD2D_Vector (Pos_Vector : in Pos2D_Vector) return AD2D_Vector is
      Rvec   : AD_Vector := Var (To_Array (Pos_Vector));
      Result : AD2D_Vector (Pos_Vector'First .. Pos_Vector'Last);
      Idx    : Int := Pos_Vector'First;
   begin
      for K in Rvec'Range loop
	 if K mod 2 = 1 then  		--  K is odd
	    Result (Idx).X := Rvec (K);
	 else				--  K is even
	    Result (Idx).Y := Rvec (K);
	    Idx := Idx + 1;
	 end if;
      end loop;
      return Result;
   end To_AD2D_Vector;

   
end Forward_AD.AD2D;
