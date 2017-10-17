with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Henon_Heiles is
   use Int_IO, Real_IO, Real_Functions;
   N : constant Nat  := 2;
   K : constant Nat  := 13;
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
   
   function Func (X : in Vector) return Real is
   begin
      return X (1);
   end Func;
   
   function Sgn (X : in Real) return Real is
   begin
      if X >= 0.0 then return 1.0;
      else return -1.0;
      end if;
   end Sgn;
   
   function Find_State_At_Level (Level : in Real;
				 A : in Array_Of_Vectors;
				 T  : in Real;
				 Dt : in Real;
				 Lower : in out Real;
				 Upper : in out Real;
				 Func : not null access function (X : Vector)
				   return Real) return Variable is
      Guess : Variable;
      Est   : Real := 1.0e10;
      Sign  : constant Real := Sgn (Func (Interpolate (A, Upper, T, T + Dt)) -
				      Func (Interpolate (A, Lower, T, T + Dt)));
   begin
      while abs (Est - Level) > 1.0e-10 loop
	 Guess.T  := 0.5 * (Lower + Upper);
	 Guess.X := Interpolate (A, Guess.T, T, T + Dt);
	 Est      := Func (Guess.X);
	 if Est * sign > 0.0 then Upper := Guess.T;
	 else Lower := Guess.T; end if;
      end loop;
      return Guess;
   end Find_State_At_Level;
   
   -- Initial Conditions ----
   Var, State, Guess : Variable := (0.0, (0.0, -0.1, 1.0, 0.0));
   X : Vector renames Var.X;
   T : Real renames Var.T;
   -------------------------------
   
   Y    : Real_Vector (1 .. 2 * N * K);
   A, Q : Array_Of_Vectors;
   Fcsv : File_Type;
   H0 : constant Real := 1.0 / 12.0;
   T_Final : constant Real := 1_000.0;
   Lower, Upper : Real;
   AD : AD_Type;

begin
   Control.Max_Dt := 100.0;
   X     := Get_IC (X, H0);
   AD    := Hamiltonian (0.0, X);
   State := Var;
   
   Put (Val (AD)); New_Line;
   Create (Fcsv, Name => "out.csv");
   Put_Line (Fcsv, "t, q1, q2, q_dot1, q_dot2, p1, p2, E");
   Print_Lagrangian (Fcsv, Var, Lagrangian'Access);
   
   Put (Var.T); New_Line;
   while T < T_Final loop
      Y := Update (Lagrangian'Access, Var, Control, Sparse);
      A := Chebyshev_Transform (Y);
      ----------------------------------------------------------------------  
      Q := Split (Y);
      for I in 2 .. K loop
      	 if Q (1) (I - 1) * Q (1) (I) < 0.0 then -- If there's a zero, bisect
	    Lower   := Var.T + Control.Dt * Grid (I - 1);
	    Upper   := Var.T + Control.Dt * Grid (I);
	    Guess   := Find_State_At_Level
	      (0.0, A, Var.T, Control.Dt, Lower, Upper, Func'Access);
	    ----------------------------------------------------------------
	    if Guess.X (3) > 0.0 then
	       Print_Lagrangian (Fcsv, Guess, Lagrangian'Access);
	    end if;
      	 end if;
      end loop;
      ----------------------------------------------------------------------  
      Put (Var.T); New_Line;
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;
   Close (Fcsv);
end Henon_Heiles;
