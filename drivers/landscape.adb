with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

procedure Landscape is
   use Real_IO, Int_IO, Real_Functions;
   α : Real := 100.0;
   
   function PE (Q : in Real_Vector) return Real is
      PE_G, PE_M : Real;
      T     : Real renames Q (1);
      S     : Real renames Q (2);
      Tmp   : Real := Cos (0.5 * (T + S));
      Cs    : constant Real := Cos (S);
      Ct    : constant Real := Cos (T);
      C2tps : constant Real := Cos (2.0 * T +       S);
      Ctp2s : constant Real := Cos (T       + 2.0 * S);
      TpS   : constant Real := T + S;
   begin
      if TpS > 0.99 * π then
	 Tmp := 0.5 / (Tmp + 1.0e-10);
      elsif TpS < -0.99 * π then
	 Tmp := 0.5 / (Tmp - 1.0e-10);
      end if;
      --  Tmp  := (0.5 * Sign (Tmp)) / (abs (Tmp) + 1.0e-10);
      PE_G := Ct + C2tps;
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2
	   +  Cos (2.0 * S) - 3.0 * Cs ** 2
  	   + (Tmp ** 3) * Cos (2.0 * (T + S))
	   - 3.0 * (Tmp ** 5) * ((Ct + C2tps) * (Cs + Ctp2s));
      return ((α / 6.0) * PE_M + PE_G);
   end PE;
   
   
   function Coordinate_Transform (X : in Real_Vector) return Real_Vector is
      Y   : Real_Vector (1 .. 2);
   begin
      Y (1) := 0.5 * (X (1) - X (2));
      Y (2) := 0.5 * (X (1) + X (2));
      return Y;
   end Coordinate_Transform;
   
   X, Y : Real_Vector (1 .. 2);
begin
   null;
   
   X := (1.0, 0.0); Y := Coordinate_Transform (X);
   Put ("(1, 0) ----> (");
   Put (Y (1), Aft => 2, Fore => 1, Exp => 0); Put (", ");
   Put (Y (2), Aft => 2, Fore => 1, Exp => 0); Put_Line (")");
   X := (0.0, 1.0); Y := Coordinate_Transform (X);
   Put ("(0, 1) ----> (");
   Put (Y (1), Aft => 2, Fore => 1, Exp => 0); Put (", ");
   Put (Y (2), Aft => 2, Fore => 1, Exp => 0); Put_Line (")");
   X := (1.0, 1.0); Y := Coordinate_Transform (X);
   Put ("(1, 1) ----> (");
   Put (Y (1), Aft => 2, Fore => 1, Exp => 0); Put (", ");
   Put (Y (2), Aft => 2, Fore => 1, Exp => 0); Put_Line (")");
end Landscape;
