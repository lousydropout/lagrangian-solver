with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Henon_Heiles is
   use Int_IO, Real_IO, Real_Functions;
   N : constant Nat  := 2;
   K : constant Nat  := 5;
   package AD_Package is new Dense_AD (2 * N); 
   package Integrator is new AD_Package.Integrator (K);
   use AD_Package, Integrator;
   -----------------------------------------------
   Control : Control_Type := New_Control_Type;
   function KE (Q : in AD_Vector) return AD_Type is
   begin
      return 0.5 * (Q (3) ** 2 + Q (4) ** 2);
   end KE;
   function PE (Q : in AD_Vector) return AD_Type is
   begin
      return 0.5 * (Q (1) ** 2 + Q (2) ** 2) 
	+ Q (1) ** 2 * Q (2) - Q (2) ** 3 / 3.0;
   end PE;
   function Lagrangian (T : in Real;
			X : in Vector) return AD_Type is
      Q : constant AD_Vector := Var (X);
   begin
      return KE (Q) - PE (Q);
   end Lagrangian;
   function Hamiltonian (T : in Real;
			 X : in Vector) return AD_Type is
      Q : constant AD_Vector := Var (X);
   begin
      return KE (Q) + PE (Q);
   end Hamiltonian;
   -----------------------------------------------
   function Get_IC (X : in Vector;
		    E : in Real) return Vector is
      use Real_IO;
      Y : Vector := X;
      G : Vector;
      H : AD_Type;
      F, Dw : Real := 1.0;
      W : Real renames Y (3); -- Y(3) is ω_t
   begin
      -- use Newton's method to solve for E - H = 0
      W := 1.0;
      while abs (F) > 1.0e-10 loop
	 H  := Hamiltonian (0.0, Y);
	 F := E - Val (H);
	 G  := Grad (H);
	 Dw := (E - Val (H)) / G (3); -- G(3) is \partial H / \partial ω_t
	 W  := W + Dw;
      end loop;
      H := Hamiltonian (0.0, Y);
      F := E - Val (H);
      return Y;
   end Get_IC;
   
   
   
   -- Initial Conditions ----
   Var, State, Guess : Variable := (0.0, (0.0, -0.1, 1.0, 0.0));
   X : Vector renames Var.X;
   T : Real renames Var.T;
   -------------------------------
   
   Y    : Real_Vector (1 .. 2 * N * K);
   A    : Array_Of_Vectors;
   Dt   : constant Real := 1.0e-3;
   Cname : String := "out.csv";
   Fcsv : File_Type;
   H0 : constant Real := 1.0 / 12.0;
   T_Final : constant Real := 1_000.0;
   Old : Real;
   AD : AD_Type;
begin
   X     := Get_IC (X, H0);
   AD    := Hamiltonian (0.0, X);
   State := Var;
   
   Put (Val (AD)); New_Line;
   Create (Fcsv, Name => Cname);
   Put_Line (Fcsv, "t, q1, q2, q_dot1, q_dot2, p1, p2, E");
   Print_Lagrangian (Fcsv, Var, Lagrangian'Access);

   Old := 0.0;
   while T < T_Final loop
      Put (Var.T); New_Line;
      Y := Update (Lagrangian'Access, Var, Control, Sparse);
      ----------------------------------------------------------------------  
      A := Chebyshev_Transform (Y);
      while State.T <= T + Control.Dt loop
      	 State.T := State.T + Dt;
      	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
      	 if Old * State.X (1) < 0.0 and State.X (3) > 0.0 then
      	    Print_Lagrangian (Fcsv, State, Lagrangian'Access);
      	 end if;
	 Old := State.X (1);
      end loop;
      ----------------------------------------------------------------------  
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;
   Close (Fcsv);
end Henon_Heiles;
