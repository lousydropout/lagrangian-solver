with Numerics, Ada.Text_IO, Forward_AD.AD2D;
use  Numerics, Ada.Text_IO, Forward_AD.AD2D;

procedure Forward_AD.Bubble is
   use Real_Functions, Real_IO, Int_IO;
   
   type Potential_Function is not null access 
     function (X : in Pos2D_Vector) return AD_Type;
   
   function Potential (Q : in Pos2D_Vector) return AD_Type;
   function Potential (Q : in Pos2D_Vector) return AD_Type is
      X : AD2D_Vector := To_AD2D_Vector (Q);
      U : AD_Type   := Zero (X'Length);
   begin
      for I in X'Range loop
	 for J in I + 1 .. X'Last loop
	    U := U + Square (X (I) - X (J));
	 end loop;
      end loop;
      return U;
   end Potential;
   
   
   procedure Verlet (Potential : in Potential_Function;
		     X	       : in Pos2D_Vector;
		     File      : in File_Type);
   procedure Verlet (Potential : in Potential_Function;
		     X	       : in Pos2D_Vector;
		     File      : in File_Type) is
      U : AD_Type    :=  Potential (X);
      F : Real_Array := -Grad (U);
   begin

      for Item of F loop
	 Put (File => File, Item => Item, Exp => 0, Aft => 3); 
	 New_Line (File);
      end loop;
   end Verlet;

   
   N : constant Nat := 5;
   M : constant Nat := 8;
   
   R2D  : Pos2D_Vector (1 .. 2);
   
   
   --  X    : Real_Array := (1.0, 2.0, 3.0);
   File : File_Type;
   
begin
   
   R2D (1) := (0.0, 0.0);
   R2D (2) := (0.0, 1.0);
   
   null;
   
   
   Create (File => File, Name => "out.xyz");
   
   
   Verlet (Potential => Potential'Access, 
	   X         => R2D, 
	   File      => File);
   
   
   Close (File => File);
end Forward_AD.Bubble;
   
