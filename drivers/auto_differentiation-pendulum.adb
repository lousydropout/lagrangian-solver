with Numerics, Ada.Text_IO, Auto_Differentiation.Integrator, Chebyshev;
use  Numerics, Ada.Text_IO, Auto_Differentiation.Integrator, Chebyshev;

procedure Auto_Differentiation.Pendulum is
   use Real_IO, Real_Functions;
   --  Set Up Parameters -----------------
   Control : Control_Type
     := (N => 1, Dt => 2.0, Eps => 1.0e-10, Err => 1.0, K => 11);
   N : Nat renames Control.N;
   -------------------------------
   
   function Lagrangian (X : in Real_Vector; N : in Nat) return AD_Type is
      θ : AD_Type := Var (X => X (1), I => 1, N => 2);
      ω : AD_Type := Var (X => X (2), I => 2, N => 2);
      M : Real    := 2.0;
      L : Real    := 1.0;
      G : Real    := 10.0;
   begin
      return (0.5 * M * L**2) * (ω ** 2) - M * G * (1.0 - Sin (θ));
   end Lagrangian;

   -------------------------------
   -- Initial Conditions ----
   Var : Variable :=  (N2 => 2 * N,
		       X  => (1.0e-10, 0.0),
		       T  => 0.0);
   θ : Real renames Var.X (1);
   ω : Real renames Var.X (2);
   T : Real renames Var.T;
   -------------------------------
   Y    : Real_Vector (1 .. 2 * N * Control.K);
   A, B : Real_Vector (1 .. N * Control.K);
   File : File_Type;
   Dt   : constant Real := 0.05;
   Time : Real := T - Dt;
begin
   
   Create (File, Name => "pendulum.csv");
   Put_Line (File, "time, x, u, y, v");
   
   while T < 5.0 loop
      Y := Collocation (Lagrangian'Access, Var, Control);
      for I in 1 .. Control.K loop
	 A (I) := Y (2 * I - 1);
	 B (I) := Y (2 * I);
      end loop;
      
      A := CGL_Transform (A);
      B := CGL_Transform (B);
      
      while Time <= T + Control.Dt loop
	 Time := Time + Dt;
	 θ    := Interpolate (A => A, X => Time, L => T, R => T + Control.Dt);
	 ω    := Interpolate (A => B, X => Time, L => T, R => T + Control.Dt);
	 Put (File, Time);          Put (File, ",  ");
	 Put (File,       Cos (θ)); Put (File, ",  ");
	 Put (File, -ω  * Sin (θ)); Put (File, ",  ");
	 Put (File, 1.0 - Sin (θ)); Put (File, ",  ");
	 Put (File, -ω  * Cos (θ)); New_Line (File);
      end loop;
      
      θ := Y (Y'Last - 1);
      ω := Y (Y'Last);
      T := T + Control.Dt;
   end loop;
   
   Close (File);
   
end Auto_Differentiation.Pendulum;
