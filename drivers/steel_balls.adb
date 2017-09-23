with Numerics, Ada.Text_IO, Chebyshev, Dense_AD, Dense_AD.Integrator;
use  Numerics, Ada.Text_IO, Chebyshev;

procedure Steel_Balls is
   use Int_IO, Real_IO, Real_Functions;
   -----------------------------------------------
   N   : constant Nat  := 2;
   K   : constant Nat  := 7;
   α   : constant Real := 100.0;
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
   begin
      Tmp  := 0.5 / Tmp;
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
      Q     : AD_Vector := Var  (X);
   begin
      --  pragma Assert (abs (X (1) + X (2)) < π, "singularity error: t + s = pi");
      return  KE (Q) - PE (Q);
   end Lagrangian;
   -----------------------------------------------
   -- Initial Conditions ----
   Var : Variable :=  (X  => (0.0, 1.5, 10.0, 1000.0),
		       T  => 0.0);
   State : Variable := Var;
   X : Vector renames Var.X;
   T : Real renames Var.T;
   -------------------------------

   Y    : Real_Vector (1 .. 2 * N * K);
   A    : Array_Of_Vectors;
   File : File_Type;
   Dt   : constant Real := 1.0e-4;
   Name : String := "out.xyz";
   
begin

   Print_XYZ (File, Var, Name, Create);
   
   while T < 1.0e-2 loop
      Y := Update (Lagrangian'Access, Var, Control);
      Put (Var.T); New_Line;
      
      A := Chebyshev_Transform (Y);
      while State.T <= T + Control.Dt loop
      	 State.T := State.T + Dt;
	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
	 Print_XYZ (File, Var, Name);
      end loop;
      
      
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;

end Steel_Balls;
