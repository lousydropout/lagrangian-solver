with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Pendulum is
   use Int_IO, Real_IO, Real_Functions;
   -----------------------------------------------
   N   : constant Nat := 1;
   K   : constant Nat := 7;
   Num : constant Nat := 2 * N;
   -----------------------------------------------
   package AD_Package is new Dense_AD (Num); 
   package Integrator is new AD_Package.Integrator (K);
   use AD_Package, Integrator;
   -----------------------------------------------
   Control : Control_Type := New_Control_Type (Tol => 1.0e-10);
   -----------------------------------------------
   function Lagrangian (T : in Real;
			X : in Vector) return AD_Type is
      Q : AD_Vector := Var (X);
      θ : AD_Type renames Q (1);
      ω : AD_Type renames Q (2);
      M : Real    := 2.0;
      L : Real    := 1.0;
      G : Real    := 10.0;
   begin
      return (0.5 * M * L**2) * (ω ** 2) - M * G * (1.0 - Sin (θ));
   end Lagrangian;
   -----------------------------------------------
   -- Initial Conditions ----
   Var  : Variable :=  (X => (1.0e-10, 0.0),
			T => 0.0);
   State : Variable := Var;
   θ     : Real renames State.X (1);
   ω     : Real renames State.X (2);
   -------------------------------

   Y		: Real_Vector (1 .. Num * K);
   A		: Array_Of_Vectors;
   File, Phase  : File_Type;
   Dt		: constant Real := 0.01;
   
begin

   Create (File, Name => "pendulum.csv");
   Put_Line (File, "time, x, u, y, v");
   Create (Phase, Name => "phase.csv");
   Put_Line (Phase, "time, q, q_dot, p, H");
   
   while Var.T < 10.0 loop
      Y := Update (Lagrangian'Access, Var, Control, Dense);
      
      A := Chebyshev_Transform (Y);
      while State.T <= Var.T + Control.Dt loop
   	 State.T := State.T + Dt;
	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
	 -- Print data --------------------------------------------------
	 Print_Lagrangian (Phase, State, Lagrangian'Access, 2, 3, 2);
   	 Put (File, State.T);       Put (File, ",  ");
   	 Put (File,       Cos (θ)); Put (File, ",  ");
   	 Put (File, -ω  * Sin (θ)); Put (File, ",  ");
   	 Put (File, 1.0 - Sin (θ)); Put (File, ",  ");
   	 Put (File, -ω  * Cos (θ)); New_Line (File);
	 ----------------------------------------------------------------
      end loop;
      
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;

   Close (File);
   Close (Phase);
   
end Pendulum;
