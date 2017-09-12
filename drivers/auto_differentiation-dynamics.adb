with Numerics, Ada.Text_IO, Auto_Differentiation.Integrator;
use  Numerics, Ada.Text_IO, Auto_Differentiation.Integrator;

procedure Auto_Differentiation.Dynamics is
   use Real_IO, Int_IO;
   --  Set Up Parameters -----------------
   N       : constant Nat := 2;
   Control : Control_Type (N => N);
   α       : constant Real := 100.0;
   -------------------------------
   function KE (Q, P : in AD_Vector; N : in Nat) return AD_Type is
      C    : constant AD_Type := Cos (Q (1) + Q (2));
      Cp2  : constant AD_Type := C + 2.0;
      C4p6 : constant AD_Type := 4.0 * C + 6.1;
      Det  : constant AD_Type := 1.1 * C4p6 - Cp2 ** 2;
      T_Dot, S_Dot : AD_Type;
   begin
      T_Dot :=  (1.1 * P (1) -  Cp2 * P (2)) / Det;
      S_Dot := (-Cp2 * P (1) + C4p6 * P (2)) / Det;
      return (0.5 * C4p6 * T_Dot ** 2 + 0.55 * S_Dot ** 2 + Cp2 * T_Dot * S_Dot);
   end KE;
   
   function PE (Q : in AD_Vector; N : in Nat) return AD_Type is
      PE_G, PE_M : AD_Type;
      T     : AD_Type renames Q (1);
      S     : AD_Type renames Q (2);
      Tmp   : AD_Type := Cos (0.5 * (T + S));
      Cs    : constant AD_Type := Cos (S);
      Ct    : constant AD_Type := Cos (T);
      C2tps : constant AD_Type := Cos (2.0 * T +       S);
      Ctp2s : constant AD_Type := Cos (T       + 2.0 * S);
   begin
      Tmp  := 0.5 / Tmp;
      PE_G := Ct + C2tps;
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2
	   +  Cos (2.0 * S) - 3.0 * Cs ** 2
  	   + (Tmp ** 3) * Cos (2.0 * (T + S))
	   - 3.0 * (Tmp ** 5) * ((Ct + C2tps) * (Cs + Ctp2s));
      return ((α / 6.0) * PE_M + PE_G);
   end PE;
   
   --- Set up Hamiltonian -----
   function Hamiltonian (X : in Real_Array; N : in Nat) return AD_Type is
      Q : AD_Vector := Var  (X (1     ..     N), 2 * N,     1);
      P : AD_Vector := Var  (X (N + 1 .. 2 * N), 2 * N, N + 1);
   begin
      pragma Assert (abs (X (1) + X (2)) < π, "singularity error: t + s = pi");
      return (KE (Q, P, N) + PE (Q, N));
   end Hamiltonian;
   -------------------------------
   
   function Momenta (Var : in Variable) return Real_Array is
      use Real_Functions;
      X : Real_Array renames Var.X;
      C : constant Real := Cos (X (1) + X (2));
   begin
      return (X (1), 
	      X (2),
	      (6.1 + 4.0 * C) * X (3) + (2.0 + C) * X (4),
	      (2.0 +       C) * X (3) +  1.1      * X (4));
   end Momenta;
   
   
   -- Initial Conditions ----
   Var : Variable :=  (N2 => 2 * N,
		       X  => (0.0, 0.0, -2.0, 10.0),
		       T  => 0.0);
   X : Real_Array renames Var.X;
   T : Real       renames Var.T;
   -------------------------------
   Time : Real;
   XYZ  : File_Type;
   DT   : constant Real := 1.0e-2;
   
begin
   Control.Eps := 1.0e-8;
   Level := Gradient;
   
   X := Momenta (Var);
   
   Create (XYZ, Name => "out.xyz");

   Put_Line ("T, x1, y1, x2, y2, E");
   Print_Data (Var, Hamiltonian'Access);
   Print_XYZ (XYZ, Var);
   
   while T < 1.0e2 loop
      Time := T + Dt;
      while T < Time loop
	 if Control.Dt > Time - T then Control.Dt := Time - T; end if;
	 Update (Hamiltonian'Access, Var, Control);
      end loop;
      Print_Data (Var, Hamiltonian'Access);
      Print_XYZ (XYZ, Var);
   end loop;
   
end Auto_Differentiation.Dynamics;
