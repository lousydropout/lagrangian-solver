with Numerics, Ada.Text_IO, Auto_Differentiation.Integrator;
use  Numerics, Ada.Text_IO, Auto_Differentiation.Integrator;

procedure Auto_Differentiation.Dynamics is
   use Real_IO, Int_IO;
   --  Set Up Parameters -----------------
   N       : constant Nat := 2;
   Control : Control_Type (N => N);
   -------------------------------
   function KE (Q, P : in AD_Vector; N : in Nat) return AD_Type is
      C        : constant AD_Type := Cos (Q (1) + Q (2));
      Const6_1 : constant AD_Type := Const (6.1, 2 * N);
      Const1_1 : constant AD_Type := Const (1.1, 2 * N);
      C_Plus_2 : constant AD_Type := Const (2.0, 2 * N) + C;
      C_Plus_6 : constant AD_Type := Const6_1 + 4.0 * C;
      Det : constant AD_Type := 1.1 * C_Plus_6 - C_Plus_2 ** 2;
      T_Dot, S_Dot : AD_Type;
   begin
      T_Dot := 1.1 * P (1) - C_Plus_2 * P (2);
      S_Dot := -C_Plus_2 * P (1) + C_Plus_6 * P (2);
      T_Dot := T_Dot / Det;
      S_Dot := S_Dot / Det;
      
      return (0.5 * C_Plus_6 * T_Dot ** 2 + 0.55 * S_Dot ** 2
		+ C_Plus_2 * T_Dot * S_Dot);
   end KE;
   
   function PE (Q : in AD_Vector; N : in Nat) return AD_Type is
      PE_G, PE_M : AD_Type;
      T : AD_Type renames Q (1);
      S : AD_Type renames Q (2);
      Tmp : AD_Type := Cos (0.5 * (T + S));
      Alpha : constant Real := 10.0;
   begin
      pragma Assert (Val (Tmp) /= 0.0);
      Tmp := Const (0.5, 2 * N) / Tmp;
      PE_G := Cos (T) + Cos (2.0 * T + S);
      PE_M
	:= Cos (2.0 * T) - 3.0 * Cos (T) ** 2
	+  Cos (2.0 * S) - 3.0 * Cos (S) ** 2
	+ (Tmp ** 3) * Cos (2.0 * (T + S))
	- 3.0 * (Tmp ** 5) * ((Cos (T) + Cos (2.0 * T + S)) 
			      * (Cos (S) + Cos (2.0 * S + T)));
      return ((Alpha / 6.0) * PE_M + PE_G);
   end PE;
   

   function Hamiltonian (X : in Real_Array; N : in Nat) return AD_Type;
   --- print data ------
   procedure Print_Data (Var : in Variable) is
      use Real_Functions;
      T : Real renames Var.X (1);
      S : Real renames Var.X (2);
      X, Y : Real_Array (1 .. 2);
   begin
      X (1) := -Sin (T);
      Y (1) :=  Cos (T);
      X (2) := X (1) - Sin (2.0 * T + S);
      Y (2) := Y (1) + Cos (2.0 * T + S);
      
      ---------------------------------------
      Put (Var.T, Aft => 3, Exp => 0); -- print time
      for I in 1 .. Nat (2) loop
	 Put (",  "); Put (X (I), Aft => 3, Exp => 0);
	 Put (",  "); Put (Y (I), Aft => 3, Exp => 0); -- print positions
      end loop;
      Put (",  "); 
      -- print total energy
      Put (Val (Hamiltonian (Var.X, N)), Aft => 10, Exp => 0); New_Line;
   end Print_Data;
   
   -------------------------------
   --- Set up Hamiltonian -----
   function Hamiltonian (X : in Real_Array; N : in Nat) return AD_Type is
      Q : AD_Vector := Var  (X (1     ..     N), 2 * N,     1);
      P : AD_Vector := Var  (X (N + 1 .. 2 * N), 2 * N, N + 1);
   begin
      return (KE (Q, P, N) + PE (Q, N));
   end Hamiltonian;
   -------------------------------
   
   function Momenta (Var : in Variable) return Real_Array is
      use Real_Functions;
      C : constant Real := Cos (Var.X (1) + Var.X (2));
      P : Real_Array := Var.X;
   begin
      P (3) := (6.1 + 4.0 * C) * Var.X (3) + (2.0 + C) * Var.X (4);
      P (4) := (2.0 + C)       * Var.X (3) + 1.1       * Var.X (4);
      return P;
   end Momenta;
   
   
   -- Initial Conditions ----
   Var : Variable :=  (N2 => 2 * N,
		       X  => (0.0, 0.0, 0.05, 0.0),
		       T  => 0.0);
   X : Real_Array renames Var.X;
   T : Real       renames Var.T;
   -------------------------------
   Time : Real;
   XYZ : File_Type;
   DT : constant Real := 1.0e-2;
begin
   Level := Gradient;
   Var.X := Momenta (Var);
   
   Create (XYZ, Name => "out.xyz");
   
   Put_Line ("T, x1, y1, x2, y2, E");
   Print_Data (Var);
   
   while T < 1.0e1 loop
      Time := T + Dt;
      while T < Time loop
	 if Control.Dt > Time - T then Control.Dt := Time - T; end if;
	 Update (Hamiltonian'Access, Var, Control);
      end loop;
      Print_Data (Var);
      Print_XYZ (XYZ, Var);
   end loop;
   

   
end Auto_Differentiation.Dynamics;
