with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Pendulum_DAE is
   use Int_IO, Real_IO, Real_Functions;
   -----------------------------------------------
   N   : constant Nat := 2;
   Nc  : constant Nat := 1;
   K   : constant Nat := 7;
   Num : constant Nat := 2 * N + Nc;
   -----------------------------------------------
   package AD_Package is new Dense_AD (Num); 
   package Integrator is new AD_Package.Integrator (K, Nc);
   use AD_Package, Integrator;
   -----------------------------------------------
   Control : Control_Type := New_Control_Type (Tol => 1.0e-10);
   -----------------------------------------------
   function Lagrangian (T : in Real;
			Q : in Vector) return AD_Type is
      X : AD_Type := Var (Q (1), 1);
      Y : AD_Type := Var (Q (2), 2);
      U : AD_Type := Var (Q (3), 3);
      V : AD_Type := Var (Q (4), 4);
      λ : AD_Type := Var (Q (5), 5);
   begin
      return 0.5 * (U * U + V * V) - 10.0 * Y - λ * (100.0 - X * X - Y * Y);
   end Lagrangian;
   -----------------------------------------------
   -- Initial Conditions ----
   Var  : Variable :=  (X => (10.0, 0.0, 0.0, 0.0, 1.0e-20),
			T => 0.0);
   State : Variable := Var;
   -------------------------------

   Y		: Real_Vector (1 .. Num * K);
   A		: Array_Of_Vectors;
   File         : File_Type;
   Dt		: constant Real := 0.01;
   
begin

   Create (File, Name => "pendulum.csv");
   Put_Line (File, "time, x, u, y, v");
   
   while Var.T < 10.0 loop
      Y := Update (Lagrangian'Access, Var, Control, Dense);
      
      A := Chebyshev_Transform (Y);
      while State.T <= Var.T + Control.Dt loop
   	 State.T := State.T + Dt;
	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
	 -- Print data --------------------------------------------------
   	 Put (File, State.T);       Put (File, ",  ");
   	 Put (File, State.X (1)); Put (File, ",  ");
   	 Put (File, State.X (2)); Put (File, ",  ");
   	 Put (File, State.X (3)); Put (File, ",  ");
   	 Put (File, State.X (4)); New_Line (File);
	 ----------------------------------------------------------------
      end loop;
      
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;

   Close (File);
   
end Pendulum_DAE;
