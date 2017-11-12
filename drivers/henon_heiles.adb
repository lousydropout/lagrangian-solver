with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

with Numerics.Sparse_Matrices; use Numerics.Sparse_Matrices;
procedure Henon_Heiles is
   use Int_IO, Real_IO, Real_Functions;
   N : constant Nat  := 2;
   K : constant Nat  := 13;
   package AD_Package is new Dense_AD (2 * N); use AD_Package;
   -----------------------------------------------
   function KE (Q : in Real_Vector) return Real is
   begin
      return 0.5 * (Q (3) ** 2 + Q (4) ** 2);
   end KE;
   function PE (Q : in Real_Vector) return Real is
   begin
      return 0.5 * (Q (1) ** 2 + Q (2) ** 2) 
	+ Q (1) ** 2 * Q (2) - Q (2) ** 3 / 3.0;
   end PE;
   function Hamiltonian (T : in Real;
			 X : in Real_Vector) return Real is
   begin
      return KE (X) + PE (X);
   end Hamiltonian;
   function Lagrangian (X : in Real_Vector) return AD_Type is
      Q : AD_Vector := Var (X);
      T : AD_Type := 0.5 * (Q (3) ** 2 + Q (4) ** 2);
      V : AD_TYpe := 0.5 * (Q (1) ** 2 + Q (2) ** 2) 
	+ Q (1) ** 2 * Q (2) - Q (2) ** 3 / 3.0;
   begin
      return T - V;
   end Lagrangian;
   -----------------------------------------------
   function Gradient (T : in Real;
		      Z : in Real_Vector) return Real_Vector is
      X : Real renames Z (1);
      Y : Real renames Z (2);
      U : Real renames Z (3);
      V : Real renames Z (4);
   begin
      return (-X * (1.0 + 2.0 * Y),
	      Y * (Y - 1.0) - X ** 2,
	      U,
	      V);
   end Gradient;
   function Hessian (T : in Real;
		     Z : in Real_Vector) return Real_Matrix is
      X : Real renames Z (1);
      Y : Real renames Z (2);
      U : Real renames Z (3);
      V : Real renames Z (4);
      Mat : Real_Matrix (1 .. 4, 1 .. 4) := (others => (others => 0.0));
   begin
      Mat (1, 1) := -(1.0 + 2.0 * Y);
      Mat (1, 2) := -2.0 * X;
      Mat (2, 1) := Mat (1, 2);
      Mat (2, 2) := 2.0 * Y - 1.0;
      Mat (3, 3) := 1.0;
      Mat (4, 4) := 1.0;
      return Mat;
   end Hessian;
   procedure Gradient (T : in     Real;
		       Z : in     Real_Vector;
		       G :    out Real_Vector) is
      L : AD_Type;
      Old : Evaluation_Level := Get_Evaluation_Level;
   begin
      Set_Evaluation_Level (Gradient);
      L := Lagrangian (Z);
      G := Grad (L);
      Set_Evaluation_Level (Old);
      --  G := Gradient (T, Z);
   end Gradient;
   procedure Deriv (T : in     Real;
		    Z : in     Real_Vector;
		    G :    out Real_Vector;
		    H :    out Real_Matrix) is
      L : AD_Type;
      Old : Evaluation_Level := Get_Evaluation_Level;
   begin
      L := Lagrangian (Z);
      G := Grad (L);
      H := Hessian (L);
      Set_Evaluation_Level (Old);
      --  G := Gradient (T, Z);
      --  H := Hessian (T, Z);
   end Deriv;
   -----------------------------------------------
   package Integrator_Pkg is new Integrator (K => K,
					     Num => 2 * N,
					     Gradient => Gradient,
					     Deriv => Deriv);
   use Integrator_Pkg;
   Control : Control_Type := New_Control_Type (Tol => 1.0e-7);
   -----------------------------------------------
   function Get_IC (X : in Real_Vector;
		    E : in Real) return Real_Vector is
      function Hamiltonian (T : in Real;
			    X : in Vector) return AD_Type is
	 Q : AD_Vector := Var (X);
      begin
	 return 0.5 * (Q (3) ** 2 + Q (4) ** 2)
	   + 0.5 * (Q (1) ** 2 + Q (2) ** 2) 
	   + Q (1) ** 2 * Q (2) - Q (2) ** 3 / 3.0;
      end Hamiltonian;
      use Real_IO;
      Y : Real_Vector := X;
      G : Real_Vector (1 .. 2 * N);
      H : AD_Type;
      F, Dw : Real := 1.0;
      W : Real renames Y (3); -- Y(3) is ω_t
   begin
      -- use Newton's method to solve for E - H = 0
      W := 1.0;
      while abs (F) > 1.0e-14 loop
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
   
   function Func (X : in Real_Vector) return Real is
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
				 Func : not null access function
				   (X : Real_Vector) return Real)
				return Variable is
				   
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
   X : Real_Vector renames Var.X;
   T : Real renames Var.T;
   -------------------------------
   
   Y    : Real_Vector (1 .. 2 * N * K);
   A, Q : Array_Of_Vectors;
   Fcsv : File_Type;
   H0 : constant Real := 1.0 / 12.0;
   T_Final : constant Real := 1_000.0;
   Lower, Upper : Real;
   Iter : Pos := 0;
begin
   Control.Max_Dt := 100.0;
   X     := Get_IC (X, H0);
   State := Var;
   
   Create (Fcsv, Name => "out.csv");
   Put_Line (Fcsv, "t, q1, q2, q_dot1, q_dot2, p1, p2, E");
   Print_Lagrangian (Fcsv, Var, Hamiltonian'Access);
   
   while T < T_Final loop
      Put (Var.T); Put (Iter); New_Line;
      Iter := Iter + 1;
      Y := Update (Var, Control, Sparse);
      A := Chebyshev_Transform (Y);
      
      Q := Split (Y);
      for I in 2 .. K loop
      	 if Q (1) (I - 1) * Q (1) (I) < 0.0 then -- If there's a zero, bisect
      	    Lower   := Var.T + Control.Dt * Grid (I - 1);
      	    Upper   := Var.T + Control.Dt * Grid (I);
      	    Guess   := Find_State_At_Level
      	      (0.0, A, Var.T, Control.Dt, Lower, Upper, Func'Access);
      	    ----------------------------------------------------------------
      	    if Guess.X (3) > 0.0 then
      	       Print_Lagrangian (Fcsv, Guess, Hamiltonian'Access);
      	    end if;
      	 end if;
      end loop;
      
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;
   Close (Fcsv);
end Henon_Heiles;
