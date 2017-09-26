with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Steel_Balls is
   use Int_IO, Real_IO, Real_Functions;
   -----------------------------------------------
   N   : constant Nat  := 2;
   K   : constant Nat  := 7;
   α   : Real;
   -----------------------------------------------
   package AD_Package is new Dense_AD (2 * N); 
   package Integrator is new AD_Package.Integrator (K);
   use AD_Package, Integrator;
   -----------------------------------------------
   Control : Control_Type := New_Control_Type;
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
      PE_G, PE_M : AD_Type;
      T     : AD_Type renames Q (1);
      S     : AD_Type renames Q (2);
      Tmp   : AD_Type := Cos (0.5 * (T + S));
      Cs    : constant AD_Type := Cos (S);
      Ct    : constant AD_Type := Cos (T);
      C2tps : constant AD_Type := Cos (2.0 * T +       S);
      Ctp2s : constant AD_Type := Cos (T       + 2.0 * S);
      TpS   : constant Real := Val (T) + Val (S);
   begin
      if TpS > 0.99 * π then
	 Tmp := 0.5 / (Tmp + 1.0e-10);
      elsif TpS < -0.99 * π then
	 Tmp := 0.5 / (Tmp - 1.0e-10);
      end if;
      --  Tmp  := (0.5 * Sign (Tmp)) / (abs (Tmp) + 1.0e-10);
      PE_G := Ct + C2tps;
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2
	   +  Cos (2.0 * S) - 3.0 * Cs ** 2
  	   + (Tmp ** 3) * Cos (2.0 * (T + S))
	   - 3.0 * (Tmp ** 5) * ((Ct + C2tps) * (Cs + Ctp2s));
      return ((α / 6.0) * PE_M + PE_G);
   end PE;
   -------------------------------
   function Lagrangian (T : in Real;
			X : in Vector) return AD_Type is
      Q : AD_Vector := Var  (X);
   begin
      return  KE (Q) - PE (Q);
   end Lagrangian;
   -----------------------------------------------
      
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
   T_Final  : Real := 10.0;
   Line : String (1 .. 50);
   Last : Natural;
   
begin
   
   -- Read initial conditions
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- T_init
   T     := Real'Value (Line (1 .. Last));
   Get_Line (Line, Last); -- t
   X (1) := Real'Value (Line (1 .. Last));
   Get_Line (Line, Last); -- s
   X (2) := Real'Value (Line (1 .. Last));
   Get_Line (Line, Last); -- t_dot
   X (3) := Real'Value (Line (1 .. Last));
   Get_Line (Line, Last); -- s_dot
   X (4) := Real'Value (Line (1 .. Last));
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- T_final
   T_Final := Real'Value (Line (1 .. Last));
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- dt
   Dt := Real'Value (Line (1 .. Last));
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- alpha
   α  := Real'Value (Line (1 .. Last));
   State   := Var;
   ------------------------------------------------------------
   Put ("Total energy of the system = "); 
   Put (Hamiltonian (T, X, Lagrangian'Access));
   New_Line;
   ------------------------------------------------------------
   Print_XYZ (File, Var, Name, Create);
   Print_CSV (Fcsv, Var, CName, Lagrangian'Access, Create);
   
   while T < T_Final loop
      Y := Update (Lagrangian'Access, Var, Control, Sparse);
      Put (Var.T); New_Line;
      
      A := Chebyshev_Transform (Y);
      while State.T <= T + Control.Dt loop
      	 State.T := State.T + Dt;
	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
	 Print_XYZ (File, Var, Name);
	 Print_CSV (Fcsv, Var, CName, Lagrangian'Access);
      end loop;
      
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;

end Steel_Balls;
