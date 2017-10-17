package body Sb_Package is
   
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
   
   function X1 (T : in Real) return Real is
   begin
      return -Real_Functions.Sin (T);
   end X1;
   
   function X1 (X : in Vector) return Real is
   begin
      return -Real_Functions.Sin (X (1));
   end X1;
   
   function Y1 (T : in Real) return Real is
   begin
      return Real_Functions.Cos (T);
   end Y1;
   
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
   
   function Sgn (X : in Real) return Real is
   begin
      if X >= 0.0 then return 1.0;
      else return -1.0;
      end if;
   end Sgn;
   
   function Func (X : in Vector) return Real is
   begin
      return X (1);
   end Func;

   
   
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
      Rh    : constant Real := Func (Interpolate (A, Upper, T, T + Dt));
      Rl    : constant Real := Func (Interpolate (A, Lower, T, T + Dt));
      Sign  : constant Real := Sgn (Rh - Rl);
      Iter  : Nat := 1;
   begin
      while abs (Est - Level) > 1.0e-13 loop
	 Guess.T := 0.5 * (Lower + Upper);
	 Guess.X := Interpolate (A, Guess.T, T, T + Dt);
	 Est     := Func (Guess.X);
	 if (Est - Level) * Sign > 0.0 then Upper := Guess.T;
	 else Lower := Guess.T; end if;
	 Put (Iter); Put (Est - Level); New_Line;
	 Iter := Iter + 1;
      end loop;
      return Guess;
   end Find_State_At_Level;
end Sb_Package;
