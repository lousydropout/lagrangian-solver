with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Steel_Balls is
   use Int_IO, Real_IO, Real_Functions;
   -----------------------------------------------
   N   : constant Nat  := 2;
   K   : constant Nat  := 8;
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
      return 0.5 * (1.0 + Tanh (30.0 * (R - 0.9)));
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
      Ro := Const (0.8);
      if Val (R) < Val (Ro) then return Const (1.0e2); end if;
      Tmp  := 1.0 / R;
      PE_G := Ct + 2.0 * C2tps;
      Vi := 0.2 * Exp (-10.0 * (R - Ro)) / (R - Ro) ** 2;
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
   
   
   
   -- Initial Conditions ----
   Var, State : Variable;
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
      
      A := Chebyshev_Transform (Y);
      while State.T <= T + Control.Dt loop
      	 State.T := State.T + Dt;
	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
	 Print_XYZ (File, State, Name);
	 Print_CSV (Fcsv, State, CName, Lagrangian'Access);
      end loop;
      
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;

end Steel_Balls;
