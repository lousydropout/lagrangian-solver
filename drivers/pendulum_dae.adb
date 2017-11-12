with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Pendulum_DAE is
   use Int_IO, Real_IO, Real_Functions;
   N  : constant Nat  := 2;
   Nc : constant Nat  := 1;
   K  : constant Nat  := 13;
   package AD_Package is new Dense_AD (2 * N + Nc); use AD_Package;
   -----------------------------------------------
   Param_G : constant Real := 10.0;
   Param_L : constant Real := 10.0;
   -----------------------------------------------
   function Lagrangian (Q : in AD_Vector) return AD_Type is
      X : AD_Type renames Q (1);
      Y : AD_Type renames Q (2);
      U : AD_Type renames Q (3);
      V : AD_Type renames Q (4);
      λ : AD_Type renames Q (5);
      L : Real renames Param_L;
      G : Real renames Param_G;
   begin
      return 0.5 * (U ** 2 + V ** 2) - G * Y + λ * (X ** 2 + Y ** 2 - L ** 2);
   end Lagrangian;
   -----------------------------------------------
   procedure Gradient (T : in     Real;
		       Z : in     Real_Vector;
		       G :    out Real_Vector) is
      Q   : AD_Vector := Var (Z);
      L   : AD_Type;
      Old : Evaluation_Level := Get_Evaluation_Level;
   begin
      Set_Evaluation_Level (Gradient);
      L := Lagrangian (Q);
      G := Grad (L);
      Set_Evaluation_Level (Old);
   end Gradient;
   procedure Deriv (T : in     Real;
		    Z : in     Real_Vector;
		    G :    out Real_Vector;
		    H :    out Real_Matrix) is
      Q   : AD_Vector := Var (Z);
      L   : AD_Type;
      Old : Evaluation_Level := Get_Evaluation_Level;
   begin
      Set_Evaluation_Level (Hessian);
      L := Lagrangian (Q);
      G := Grad (L);
      H := Hessian (L);
      Set_Evaluation_Level (Old);
   end Deriv;
   -----------------------------------------------
   package Integrator_Pkg is new Integrator (K => K,
					     Num => 2 * N + Nc,
					     N_Constraints => Nc,
					     Gradient => Gradient,
					     Deriv => Deriv);
   use Integrator_Pkg;
   Control : Control_Type := New_Control_Type (Tol => 1.0e-3);
   -----------------------------------------------
   function Set_IC (Z : in Real_Vector) return Real_Vector is
      X : Real renames Z (1);
      Y : Real renames Z (2);
      U : Real renames Z (3);
      V : Real renames Z (4);
      L : Real renames Param_L;
      G : Real renames Param_G;
      Result : Real_Vector := Z;
   begin
      Result (2) := Sqrt (L ** 2 - X ** 2);
      Result (5) := 0.5 * (U ** 2 + V ** 2 - G * Y) / L ** 2;
      return Result;
   end Set_IC;
   -- Initial Conditions ----
   Var, State : Variable := (0.0, (0.4 * Param_L, 0.1, 0.0, 0.0, 0.0));
   X : Real_Vector renames Var.X;
   T : Real renames Var.T;
   -----------------------------------------------
   File : File_Type;
   T_Output : Real;
   Dt : constant Real := 0.1;
   T_Final : constant Real := 0.1;
   -----------------------------------------------
   Y    : Real_Vector (1 .. NK);
   A    : Array_Of_Vectors;
   Iter : Pos := 0;
begin
   X := Set_IC (X);
   for Item of X loop
      Put (Item); New_Line;
   end loop;
   State := Var;
   
   Create (File, Name => "pendulum.csv");
   Put_Line (File, "t, x, y, u, v, lambda");
   
   while T < T_Final loop
      Put (Var.T); Put (Iter); New_Line; Iter := Iter + 1;
      Y := Update (Var, Control, Sparse);
      A := Chebyshev_Transform (Y);
      
      T_Output := T + Dt;
      while State.T + Dt < T + Control.Dt loop
	 State.T := State.T + Dt;
   	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
	 
	 Put (File, State.T); 
	 for Item of State.X loop
	    Put (File, ", "); Put (File, Item);
	 end loop;
	 New_Line (File);
      end loop;
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;
   
   Close (File);
   
end Pendulum_DAE;
