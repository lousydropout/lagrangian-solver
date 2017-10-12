with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Steel_Balls is
   use Int_IO, Real_IO, Real_Functions;
   -----------------------------------------------
   N   : constant Nat  := 2;
   K   : constant Nat  := 11;
   α   : Real;
   -----------------------------------------------
   package E_Solver   is new Dense_AD (1); 
   package AD_Package is new Dense_AD (2 * N); 
   package Integrator is new AD_Package.Integrator (K);
   use AD_Package, Integrator;
   -----------------------------------------------
   Control : Control_Type := New_Control_Type;
   -----------------------------------------------
   function Phi (R : in AD_Type) return AD_Type is
   begin
      return 0.5 * (1.0 + Tanh (50.0 * (R - 0.5)));
   end Phi;
   -----------------------------------------------
   function KE (Q : in AD_Vector) return AD_Type is
      C    : constant AD_Type := Cos (Q (1) + Q (2));
      Cp2  : constant AD_Type := C + 2.0;
      C4p6 : constant AD_Type := 4.0 * C + 6.1;
      T_Dot : AD_Type renames Q (3);
      S_Dot : AD_Type renames Q (4);
   begin
      return (0.5 * C4p6 * T_Dot ** 2 + 0.55 * S_Dot ** 2 + Cp2 * T_Dot * S_Dot);
   end KE;
   -------------------------------
   function PE (Q : in AD_Vector) return AD_Type is
      use Real_Functions;
      PE_G, PE_M : AD_Type;
      T     : AD_Type renames Q (1);
      S     : AD_Type renames Q (2);
      R     : constant AD_Type := 2.0 * Cos (0.5 * (T + S));
      Cs    : constant AD_Type := Cos (S);
      Ct    : constant AD_Type := Cos (T);
      C2tps : constant AD_Type := Cos (2.0 * T +       S);
      Ctp2s : constant AD_Type := Cos (T       + 2.0 * S);
      TpS   : constant Real := Val (T) + Val (S);
      Vo, Vi, Ro, Tmp : AD_Type;
   begin
      Ro := Const (0.9);
      if Val (R) < Val (Ro) then return Const (1.0e2); end if;
      Tmp  := 1.0 / R;
      PE_G := Ct + 2.0 * C2tps;
      Vi := 0.01 * Exp (-12.0 * (R - Ro)) / (R - Ro) ** 2;
      Vo := (Tmp ** 3) * Cos (2.0 * (T + S))
	- 3.0 * (Tmp ** 5) * ((Ct + C2tps) * (Cs + Ctp2s));
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2
	   +  Cos (2.0 * S) - 3.0 * Cs ** 2
	   +  Vo * Phi (R) + Vi * (1.0 - Phi (R));
      return ((α / 6.0) * PE_M + PE_G);
   end PE;
   -------------------------------
   function Lagrangian (T : in Real;
			X : in Vector) return AD_Type is
      Q : AD_Vector := Var (X);
   begin
      return  KE (Q) - PE (Q);
   end Lagrangian;
   -----------------------------------------------
   function Hamiltonian (T : in Real;
			 X : in Vector) return AD_Type is
      Q : AD_Vector := Var (X);
   begin
      return  KE (Q) + PE (Q);
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
   
   function R13 (X : in Vector) return Real_Vector is
      use Real_Functions;
      T   : Real renames X (1);
      S   : Real renames X (2);
      R13 : Real_Vector (1 .. 2);
   begin
      R13 (1) := -Sin (T) - Sin (2.0 * T + S);
      R13 (2) :=  Cos (T) + Cos (2.0 * T + S);
      return R13;
   end R13;
   
   function R13 (X : in Vector) return Real is
      V : Real_Vector := R13 (X);
   begin
      return Norm (R13 (X));
   end R13;
   
   function V13 (X : in Vector) return Real_Vector is
      use Real_Functions;
      T   : Real renames X (1);
      S   : Real renames X (2);
      Ω_t : Real renames X (3);
      Ω_s : Real renames X (4);
      Ct  : constant Real := Cos (T);
      St  : constant Real := Sin (T);
      C2tps : constant Real := Cos (2.0 * T + S);
      S2tps : constant Real := Sin (2.0 * T + S);
      V13 : Real_Vector (1 .. 2);
   begin
      V13 (1) := -(Ct + 2.0 * C2tps) * Ω_T - C2tps * Ω_S;
      V13 (2) := -(St + 2.0 * S2tps) * Ω_T - S2tps * Ω_S;
      return V13;
   end V13;
   
   function V13_New (X : in Vector) return Real_Vector is
      N : Real_Vector := R13 (X);
      V : Real_Vector := V13 (X);
   begin
      N := (1.0 / Norm (N)) * N;
      return V - (2.0 * Dot (V, N)) * N;
   end V13_New;
   
   function New_Vel (X : in Vector) return Vector is
      T     : Real renames X (1);
      S     : Real renames X (2);
      Ct    : constant Real := Cos (T);
      St    : constant Real := Sin (T);
      Stps  : constant Real := Sin (T + S);
      C2tps : constant Real := Cos (2.0 * T + S);
      S2tps : constant Real := Sin (2.0 * T + S);
      Vel   : Real_Vector := V13_New (X);
      Nvel  : Real_Vector (1 .. 2);
   begin
      Nvel (1) := S2tps * Vel (1) - C2tps * Vel (2);
      Nvel (2) := -(St + S2tps) * Vel (1) + (Ct + C2tps) * Vel (2);
      Nvel     := (-1.0 / Stps) * Nvel;
      return (X (1), X (2), Nvel (1), Nvel (2));
   end New_Vel;
   
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
      Est   : Real;
      Rh    : constant Real := Func (Interpolate (A, Upper, T, T + Dt));
      Rl    : constant Real := Func (Interpolate (A, Lower, T, T + Dt));
      Sign  : constant Real := Sgn (Rh - Rl);
      Iter : Nat := 1;
   begin
      pragma Assert ((Rh - 1.0) * (Rl - 1.0) < 0.0);
      Guess.T := 0.5 * (Lower + Upper);
      Guess.X := Interpolate (A, Guess.T, T, T + Dt);
      Est     := Func (Guess.X);
      while abs (Est - Level - 1.1e-10) > 1.0e-10 loop
	 if (Est - Level - 1.1e-10) * sign > 0.0 then Upper := Guess.T;
	 else Lower := Guess.T; end if;
	 Guess.T := 0.5 * (Lower + Upper);
	 Guess.X := Interpolate (A, Guess.T, T, T + Dt);
	 Est     := Func (Guess.X);
	 Put (Iter); New_Line;
	 Iter := Iter + 1;
      end loop;
      return Guess;
   end Find_State_At_Level;

   -- Initial Conditions ----
   Guess, Var, State : Variable;
   X : Vector renames Var.X;
   T : Real renames Var.T;
   -------------------------------

   Y    : Real_Vector (1 .. 2 * N * K);
   A    : Array_Of_Vectors;
   File : File_Type;
   Fcsv : File_Type;
   Dt   : Real;
   Name : String := "out.xyz";
   Cname : String := "out.csv";
   Total_Energy, T_Final  : Real;
   Line : String (1 .. 50);
   Last : Natural;
   Okay : Boolean;
   --  Upper, Lower : Real;
   
begin
   Control.Max_Dt := 1.0e2;
   Control.Tol := 1.0e-12;
   ------------------------------------------------------------
   -- Read initial conditions
   T     := 0.0;
   X (1) := 0.0;
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- s
   X (2) := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- ω_s
   X (4) := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- total energy
   Total_Energy := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- alpha
   α  := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- dt
   Dt := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- T_final
   T_Final := Real'Value (Line (1 .. Last));
   ------------------------------------------------------------
   X := Get_IC (X, Total_Energy);
   Total_Energy := Val (Hamiltonian (0.0, X));
   Put ("Total Energy = "); Put (Total_Energy); New_Line;
   ------------------------------------------------------------
   State   := Var;
   ------------------------------------------------------------
   ------------------------------------------------------------
   Print_XYZ (File, State, Name, Create);
   Print_CSV (Fcsv, State, CName, Lagrangian'Access, Create);
   
   while T < T_Final loop
      Y := Update (Lagrangian'Access, Var, Control, Sparse);
      Put (Var.T); New_Line;
      
      Okay := True;
      A := Chebyshev_Transform (Y);
      while State.T + Dt <= T + Control.Dt and Okay loop
      	 State.T := State.T + Dt;
	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
	 ---------------------------------------------------------
	 if R13 (State.X) < 1.0 then
	    Put ("*****   ");
	    Put (State.T);
	    Put_Line ("*****   collision");
	 --     Lower   := Real'Max (Var.T, State.T - Dt);
	 --     Upper   := State.T;
	 --     Guess   := Find_State_At_Level (1.0, A, Var.T, Control.Dt, 
	 --  				    Lower, Upper, R13'Access);
	 --     Print_CSV (Fcsv, Guess, CName, Lagrangian'Access);
	 --     Var.X   := New_Vel (Guess.X);
	 --     Var.T   := Guess.T;
	 --     Print_CSV (Fcsv, Var, CName, Lagrangian'Access);
	 --     State.T := Var.T; -- State.T - Dt;
	 --     State.X := Var.X;
	 --     if R13 (State.X) < 1.0 then Put_Line ("less than 1");
	 --     else Put_Line ("Okay");
	 --     end if;
	 --     Okay    := False;
	 end if;
	 ---------------------------------------------------------
	 Print_XYZ (File, State, Name);
	 Print_CSV (Fcsv, State, CName, Lagrangian'Access);
	 
      end loop;
      -- Update variable Var
      if Okay then Update (Var => Var, Y => Y, Dt => Control.Dt); end if;
   end loop;

end Steel_Balls;
