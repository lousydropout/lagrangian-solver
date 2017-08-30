with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

procedure Dynamics is
   use Real_IO, Int_IO, Real_Functions;
   
   Dt_Exception : exception;
   
   function Func (X : in Real_Array) return Real_Array is
      Result : Real_Array (X'Range);
      Q : Real_Array renames X (1 .. 2);
      P : constant Real_Array (1 .. 2) := X (3 .. 4);
      Q_Dot : Real_Array (1 .. 2);
      Cts : constant Real := Cos (Q (1) + Q (2));
      Det : constant Real := 1.1 * (6.1 + 4.0 * Cts) - (2.0 + Cts) ** 2;
      C : constant Real := Cos (0.5 * (Q (1) + Q (2)));
      S : constant Real := Sin (0.5 * (Q (1) + Q (2)));
      A, B, D : Real;
      α : constant Real := 1.0;
      PE_G, PE_M, KE : Real_Array (1 .. 2);
   begin
      if abs (Det) < 1.0e-5 then
	 Put ("Det = "); Put (Det); New_Line;
	 Put ("t = "); Put (X (1)); New_Line;
	 Put ("s = "); Put (X (2)); New_Line;
      end if;
      pragma Assert (abs (Det) > 1.0e-5, "Determinant in Func is zero");
      Q_Dot (1) := (1.1 *               P (1) - (2.0 + Cts) * P (2)) / Det;
      Q_Dot (2) := ((6.1 + 4.0 * Cts) * P (2) - (2.0 + Cts) * P (1)) / Det;
      
      Result (1 .. 2) := Q_Dot;
      
      A := -0.25 * Sin (2.0 * (Q (1) + Q (2)))
	+ (3.0 * S) / (16.0 * C ** 4) * Cos (2.0 * (Q (1) + Q (2)))
	- (3.0 * S) / (64.0 * C ** 6)
	* (Cos (Q (1)) + Cos (2.0 * Q (1) + Q (2)))
	* (Cos (Q (2)) + Cos (2.0 * Q (2) + Q (1)));
      D := 3.0 / (32.0 * C ** 5);
      B := D * (Cos (Q (2)) + Cos (2.0 * Q (2) + Q (1)));
      D := D * (Cos (Q (1)) + Cos (2.0 * Q (1) + Q (2)));
      
      KE (1) := -Q (1) * (Q (1) + Q (2)) * Sin (Q (1) + Q (2));
      KE (2) := KE (1);
      
      PE_G (1) := -Sin (Q (1)) - 2.0 * Sin (2.0 * Q (1) + Q (2));
      PE_G (2) := -Sin (2.0 * Q (1) + Q (2));
      
      PE_M (1) := Sin (2.0 * Q (1)) 
	- A + B * (Sin (Q (1)) + 2.0 * Sin (2.0 * Q (1) + Q (2)))
	+ D * Sin (Q (1) + 2.0 * Q (2));
      PE_M (2) := Sin (2.0 * Q (2)) 
	- A + B * Sin (2.0 * Q (1) + Q (2))
	+ D * (Sin (Q (2)) + 2.0 * Sin (2.0 * Q (2) + Q (1)));
      
      Result (3 .. 4) := KE - PE_G - (α / 6.0) * PE_M;
      return Result;
   end Func;
   
   function Bogack_Shampine (X	 : in     Real_Array;
			     Dt	 : in     Real;
			     Err :    out Real) return Real_Array is
      K1, K2, K3, K4, Y, Z : Real_Array (X'Range);
   begin
      K1 := Func (X);
      K2 := Func (X + 0.50 * Dt * K1);
      K3 := Func (X + 0.75 * Dt * K2);
      Y  := X + (Dt / 9.0) * (2.0 * K1 + 3.0 * K2 + 4.0 * K3);
      K4 := Func (Y);
      Z  := X + (Dt / 24.0) * (7.0 * K1 + 6.0 * K2 + 8.0 * K3 + 3.0 * K4);
      
      Err := Norm (Z - Y);
      return Z;
   end Bogack_Shampine;
   
   
   procedure Update (X	: in out Real_Array;
		     T	: in out Real;
		     Dt	: in out Real) is
      Eps : constant Real := 1.0e-8;
      Err : Real := 1.0;
      Y   : Real_Array (X'Range);
      Tstep : Real := Dt;
   begin
      while Err > Eps loop
      	 Y := Bogack_Shampine (X, Tstep, Err);
      	 if (Err <= Eps) then
      	    X  := Y;
      	    T  := T + Tstep;
      	 end if;
	 if Tstep < 1.0e-12 then raise Dt_Exception; end if;
      	 Tstep := 0.8 * Tstep * (Eps / (Err + 1.0e-20)) ** 0.3;
      end loop;
      Dt := Tstep;
   end Update;

   
   procedure Print (File : in File_Type;
		    X	 : in Real_Array;
		    T	 : in Real) is
      X1, Y1, X2, Y2 : Real;
   begin
      X1 := -Sin (X (1));
      Y1 :=  Cos (X (1));
      X2 := X1 - Sin (2.0 * X (1) + X (2));
      Y2 := Y1 + Cos (2.0 * X (1) + X (2));
      -- time
      Put (File, T); Put (File, ",  ");
      -- position of ball 2
      Put (File, X1); Put (File, ",  "); 
      Put (File, Y1); Put (File, ",  ");
      -- position of ball 3
      Put (File, X2); Put (File, ",  "); 
      Put (File, Y2); New_Line (File);
   end Print;
   
   procedure Print_XYZ (File : in File_Type;
			X    : in Real_Array;
			T    : in Real) is
      X1, Y1, X2, Y2 : Real;
      R : constant Real := 10.0;
   begin
      X1 := -R * Sin (X (1));
      Y1 :=  R * Cos (X (1));
      X2 := X1 - R * Sin (2.0 * X (1) + X (2));
      Y2 := Y1 + R * Cos (2.0 * X (1) + X (2));
      -- print header
      Put_Line (File, "3");
      Put_Line (File, "Properties=pos:R:2");
      -- position of ball 1
      Put_Line (File, "0.0     0.0     5.0");
      -- position of ball 2
      Put (File => File, Item => X1);
      Put (File => File, Item => "     ");
      Put (File => File, Item => Y1);
      Put (File => File, Item => "     ");
      Put (File => File, Item => "5.0");
      New_Line (File => File);
      -- position of ball 3
      Put (File => File, Item => X2);
      Put (File => File, Item => "     ");
      Put (File => File, Item => Y2);
      Put (File => File, Item => "     ");
      Put (File => File, Item => "5.0");
      New_Line (File => File);
   end Print_XYZ;
   
   function Momenta (X : in Real_Array) return Real_Array is
      use Real_Functions;
      C : constant Real := Cos (X (1) + X (2));
      P : Real_Array := X;
   begin
      P (3) := (6.1 + 4.0 * C) * X (3) + (2.0 + C) * X (4);
      P (4) := (2.0 + C)       * X (3) + 1.1       * X (4);
      return P;
   end Momenta;

   
   -- Initial Conditions ----
   X : Real_Array := (0.0, 0.0, 0.1, 0.0);
   T : Real := 0.0;
   Dt : Real := 1.0e-1;
   File : File_Type;
   XYZ : File_Type;
   Time : Real := T;
   Iter : Nat := 1;
begin
   X := Momenta (X);
   
   
   Create (File, Name => "output.csv");
   Create (XYZ, Name => "out.xyz");
   
   
   Put_Line (File, "t, x1, y1, x2, y2");
   Print (File, X, T);
   Put (T); New_Line;
   
   while T < 1.0e2 loop
      Time := T + 1.0e-2;
      while T < Time loop
	 if Dt > Time - T then Dt := Time - T; end if;
	 Update (X, T, Dt);
      end loop;
      Print (File, X, T);
      Print_XYZ (XYZ, X, T);
      Put (T); Put ("    ");
      Put (Dt); New_Line;
   end loop;

   
end Dynamics;
