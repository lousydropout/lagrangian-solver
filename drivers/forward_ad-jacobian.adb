with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

procedure Forward_AD.Jacobian is
   use Real_IO, Int_IO;
   
   function Func (X : in Real_Array) return AD_Vector is
      Y : constant AD_Vector := Var (X);
      Zero : constant Real_Array (X'Range) := (others => 0.0);
      Result : AD_Vector := Var (Zero);
   begin
      Result (1) := Y (1) * Y (2);
      Result (2) := Y (1) * Y (3);
      Result (3) := Y (2) * Y (3);
      --  Result (4) := Y (2) * Y (3);
      return Result;
   end Func;
   
   procedure Print_AD_Type (X : in AD_Type) is
      G : Real_Array := Grad (X);
   begin
      Put (Val (X), Aft => 4, Exp => 0);
      for I in G'Range loop
	 Put (",  "); Put (G (I), Aft => 4, Exp => 0);
      end loop;
      New_Line;
   end Print_AD_Type;
   
   X : Real_Array := (1.0, 2.0, 3.0);
   Y : AD_Vector (1 .. 3);
   J : Sparse_Matrix;
begin
   Y := Func (X);
   J := Jacobian_Matrix (Y);
   
   for I in 1 .. Nat (3) loop
      Print_AD_Type (Y (I));
   end loop;
   New_Line; New_Line;
   
   Print (J);
   null;
end Forward_AD.Jacobian;

