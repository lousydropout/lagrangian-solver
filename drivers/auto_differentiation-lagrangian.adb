with Numerics, Ada.Text_IO, Auto_Differentiation.Integrator, Chebyshev;
use  Numerics, Ada.Text_IO, Auto_Differentiation.Integrator, Chebyshev;

procedure Auto_Differentiation.Lagrangian is

   --  Set Up Parameters -----------------
   Control : Control_Type
     := (N => 2, Dt => 0.005, Eps => 1.0e-10, Err => 1.0, K => 4);
   N       : Nat renames Control.N;
   α       : constant Real := 100.0;
   -------------------------------
   function KE (Q, Q_Dot : in AD_Vector; N : in Nat) return AD_Type is
      C    : constant AD_Type := Cos (Q (1) + Q (2));
      Cp2  : constant AD_Type := C + 2.0;
      C4p6 : constant AD_Type := 4.0 * C + 6.1;
      T_Dot : AD_Type renames Q_Dot (1);
      S_Dot : AD_Type renames Q_Dot (2);
   begin
      return (0.5 * C4p6 * T_Dot ** 2 + 0.55 * S_Dot ** 2 + Cp2 * T_Dot * S_Dot);
   end KE;
   -------------------------------
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
   -------------------------------
   function Lagrangian (X : in Real_Vector; N : in Nat) return AD_Type is
      Q     : AD_Vector := Var  (X (1     ..     N), 2 * N,     1);
      Q_Dot : AD_Vector := Var  (X (N + 1 .. 2 * N), 2 * N, N + 1);
   begin
      pragma Assert (abs (X (1) + X (2)) < π, "singularity error: t + s = pi");
      return  KE (Q, Q_Dot, N) - PE (Q, N);
   end Lagrangian;
   -------------------------------
   -- Initial Conditions ----
   Var : Variable :=  (N2 => 2 * N,
		       X  => (0.0, 0.0, -2.0, 10.0),
		       T  => 0.0);
   X : Real_Vector renames Var.X;
   -------------------------------
   XYZ   : File_Type;
   Y     : Real_Vector (1 .. 2 * N * Control.K);
   Dt    : constant Real := 0.01;
   Time  : Real := Var.T;
   PTime : Real;
   A     : array (1 .. 2 * N) of Real_Vector (1 .. Control.K);
begin
   
   Create (XYZ, Name => "out.xyz");

   Put_Line ("T, x1, y1, x2, y2");
   Print_Data_L (XYZ, Var);
   Print_XYZ (XYZ, Var);
   
   Setup (N, Control.K);
   
   while Var.T < 5.2 loop
      Y := Collocation (Lagrangian'Access, Var, Control);
      
      for I in 1 .. 2 * N loop
	 for K in 1 .. Control.K loop
	    A (I) (K) := Y (2 * N * (K - 1) + I);
	 end loop;
      end loop;
      
      for I in 1 .. 2 * N loop
	 A (I) := CGL_Transform (A (I));
      end loop;
      
      PTime := Var.T + Control.Dt;
      while Time < PTime loop
	 for I in 1 .. 2 * N loop
	    Var.X (I) := Interpolate (A (I), Time, Var.T, Var.T + Control.Dt);
	 end loop;
	 Var.T     := Time;
	 Print_XYZ (XYZ, Var);
	 Print_Data_L (XYZ, Var);
	 Time := Time + Dt;
      end loop;
      Var.T := PTime;
      for I in 1 .. 2 * N loop
	 Var.X (I) := Y (Y'Last - 2 * N + I);
      end loop;

   end loop;

   Print_Data_L (XYZ, Var);
   Close (XYZ);
end Auto_Differentiation.Lagrangian;
