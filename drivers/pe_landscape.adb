
with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

procedure PE_Landscape is
   use Real_IO, Int_IO, Real_Functions;
   α : constant Real := 1.0e-2;
   R : constant Real := π * 0.7028;
   
   
   function Phi (R : in Real) return Real is
   begin
      return 0.5 * (1.0 + Tanh (50.0 * (R - 0.5)));
   end Phi;
   
   function PE (Q : in Real_Vector) return Real is
      PE_G, PE_M : Real;
      T     : Real renames Q (1);
      S     : Real renames Q (2);
      R     : constant Real := 2.0 * abs (Cos (0.5 * (T + S)));
      Cs    : constant Real := Cos (S);
      Ct    : constant Real := Cos (T);
      C2tps : constant Real := Cos (2.0 * T +       S);
      Ctp2s : constant Real := Cos (T       + 2.0 * S);
      Vo, Vi, Tmp, Ro : Real;
   begin
      Ro := 0.9;
      if R < Ro then return -250.0; end if;
      Tmp := 1.0 / R;
      PE_G := Ct + 2.0 * C2tps;
      Vi := 0.01 * Exp (-12.0 * (R - Ro)) / (R - Ro) ** 2;
      Vo := (Tmp ** 3) * Cos (2.0 * (T + S))
	- 3.0 * (Tmp ** 5) * ((Ct + C2tps) * (Cs + Ctp2s));
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2
	   +  Cos (2.0 * S) - 3.0 * Cs ** 2 
	   + Vo * Phi (R) + Vi * (1.0 - Phi (R));
      return (1.0 * ((α / 6.0) * PE_M - PE_G));
   end PE;
   
   
   function PE2 (Q : in Real_Vector) return Real is
      PE_G, PE_M : Real;
      T     : Real renames Q (1);
      S     : Real renames Q (2);
      R     : constant Real := 2.0 * abs (Cos (0.5 * (T + S)));
      Cs    : constant Real := Cos (S);
      Ct    : constant Real := Cos (T);
      C2tps : constant Real := Cos (2.0 * T +       S);
      Ctp2s : constant Real := Cos (T       + 2.0 * S);
      Vo : Real;
   begin
      PE_G := Ct + 2.0 * C2tps;
      Vo := (Cos (2.0 * (T + S)) * R ** 2 - 3.0 * ((Ct + C2tps) * (Cs + Ctp2s)))
	/ R ** 5;
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2 + Cos (2.0 * S) - 3.0 * Cs ** 2 + Vo;
      return ((α / 6.0) * PE_M - PE_G);
   end PE2;
   
   
   function Coordinate_Transform (X : in Real_Vector) return Real_Vector is
      Y   : Real_Vector (1 .. 2);
   begin
      Y (1) := 0.5 *  (3.0 * X (1) + X (2));
      Y (2) := 0.5 * (-3.0 * X (1) + X (2));
      return Y;
   end Coordinate_Transform;
   
   
   
   N  : constant Nat := 1_000;
   Dx : constant Real := 2.0 / Real (N);
   X, Y : Real_Vector (1 .. 2);
   R13, Tmp : Real;
   Max_Diff : Real := 0.0;
begin
   
   for I in 1 .. N loop
      X (1) := -1.0 + (Real (I) - 0.5) * Dx;
      for J in 1 .. N loop
	 X (2) := -1.0 + (Real (J) - 0.5) * Dx;
   	 Y := R * Coordinate_Transform (X);
	 Tmp := abs ((PE (Y) - PE2 (Y)) / PE2 (Y));
	 R13 := 2.0 * Cos (0.5 * (Y (1) + Y (2)));
	 if R13 > 1.0 then
	    if Tmp > Max_Diff then
	       Max_Diff := Tmp;
	       Put (Y (1)); Put (", "); Put (Y (2)); New_Line;
	       Put ("    difference = "); Put (Max_Diff); New_Line; New_Line;
	    end if;
	 end if;
      end loop;
   end loop;
end PE_Landscape;
