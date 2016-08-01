with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

procedure Forward_AD.Verlet is
   use Real_Functions, Real_IO, Int_IO;
   
   function Potential (X : in Real_Array) return AD_Type;
   function Potential (X : in Real_Array) return AD_Type is
      Q : AD_Vector := Var (X);
      U : AD_Type   := Zero (X'Length);
   begin
      for I in Q'Range loop
	 U := U + Q (I) ** 2;
      end loop;
      return U;
   end Potential;
   
   
   X : Real_Array := (1.0, 2.0, 3.0);
   U : AD_Type := Potential (X);
   --  F : Sparse_Vector := U.Grad;
   F : Real_Array := Grad (U);

begin
   
   null;
   --  Print (U);
   for Item of F loop
      Put (Item, Exp => 0, Aft => 3); New_Line;
   end loop;
   
end Forward_AD.Verlet;
