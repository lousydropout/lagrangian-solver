with Numerics, Chebyshev, Numerics.Sparse_Matrices.CSparse;
use  Numerics, Chebyshev;

package body Auto_Differentiation.Integrator is
   
   function Is_Setup return Boolean is
      Okay : Boolean := True;
   begin
      if N_Col (MatA) = 0 or else N_Row (MatA) = 0 or else
	 N_Col (MatB) = 0 or else N_Row (MatB) = 0 or else
	 N_Col (MatC) = 0 or else N_Row (MatC) = 0 or else
	 N_Col (MatD) = 0 or else N_Row (MatD) = 0 
      then 
	 Okay := False;
      end if;
      return Okay;
   end Is_Setup;
   
   procedure Setup (N : in Nat;
		    K : in Nat) is
      Top_Left     : Sparse_Matrix := Sparse (((1.0, 0.0), (0.0, 0.0)));
      Top_Right    : Sparse_Matrix := Sparse (((0.0, 1.0), (0.0, 0.0)));
      Bottom_Left  : Sparse_Matrix := Sparse (((0.0, 0.0), (1.0, 0.0)));
      Bottom_Right : Sparse_Matrix := Sparse (((0.0, 0.0), (0.0, 1.0)));
      EyeN  : constant Sparse_Matrix := Eye (N);
      EyeK  : constant Sparse_Matrix := Eye (K);
      Eye2N : constant Sparse_Matrix := Eye (2 * N);
      D     : constant Sparse_Matrix := Sparse (Derivative_Matrix (K, 0.0, 1.0));
   begin
      Top_Left     := Top_Left     and EyeN;
      Top_Right    := Top_Right    and EyeN;
      Bottom_Left  := Bottom_Left  and EyeN;
      Bottom_Right := Bottom_Right and EyeN;
      
      MatC := D and Eye2N; -- temporary, overwritten below
      MatA := MatC * (EyeK and Top_Left);
      MatB := MatC * (EyeK and Bottom_Right);
      MatC := EyeK and Top_Right;
      MatD := EyeK and Bottom_Left;
   end Setup;
   
   
   
   function Collocation (Lagrangian : not null access 
			   function (X : Real_Vector; N : Nat) return AD_Type;
			 Var        : in     Variable;
			 Control    : in out Control_Type) return Real_Vector is
      use Ada.Text_IO, Real_IO;
      K   : Nat  renames Control.K;
      N   : Nat  renames Control.N;
      Old : constant Evaluation_Level :=  Level;
      Tmp : constant Integer := 2 * N;
      
      Q   : Real_Vector (1 .. 2 * N * K);
      DQ  : Sparse_Vector;
      F   : Sparse_Vector;
      J   : Sparse_Matrix;
      Res : Real := 1.0;
   begin
      Level := Hessian;
      ------------------------------------------------
      for I in 1 .. K loop
	 Q ((I - 1) * 2 * N + 1 .. I * 2 * N) := Var.X;
      end loop;
      ------------------------------------------------
      while Res > 1.0e-10 loop
	 FJ (Lagrangian, Var, Control, Q, F, J);
	 DQ := Numerics.Sparse_Matrices.CSparse.Solve (J, F);
	 Q (Tmp + 1 .. Tmp * K) := Q (Tmp + 1 .. Tmp * K) - To_Array (DQ);
	 Res := Norm (F);
	 --  Put (Res); New_Line;
      end loop;
      --  New_Line;
      ------------------------------------------------
      Control.Err := Norm (DQ) / Norm (Q);
      Level := Old;

      return Q;
   end Collocation;
   
   procedure FJ (Lagrangian : not null access 
		    function (X : Real_Vector; N : Nat) return AD_Type;
		  Var     : in     Variable;
		  Control : in     Control_Type;
		  Q       : in     Real_Vector;
		  F       :    out Sparse_Vector;
		  J       :    out Sparse_Matrix) is
      N    : Nat  renames Control.N;
      L    : AD_Type;
      X    : Real_Vector (1 .. 2 * N);
      Tmp  : Integer;
      U    : Sparse_Vector;
      V    : Sparse_Matrix;
      A    : constant Sparse_Matrix := MatA - Control.Dt * MatC;
      B    : constant Sparse_Matrix := MatB - Control.Dt * MatD;
   begin
      pragma Assert (Q'Length = 2 * N * Control.K);
      -------------------------------------
      X := Q (1 .. 2 * N);
      L := Lagrangian (X, N);
      U := Grad (L); 
      V := Hessian (L);
      for K in 2 .. Control.K loop
	 Tmp := 2 * N * (K - 1);
	 X := Q (Tmp + 1 .. Tmp + 2 * N);
	 L := Lagrangian (X, N);
	 U := U or Grad (L);
	 V := V or Hessian (L);
      end loop;
      -------------------------------------
      F := A * Q - B * U; F := Remove_1stN (F, 2 * N);
      J := A     - B * V; J := Remove_1stN (J, 2 * N);
   end FJ;
   
   
   
   function Bogack_Shampine (Hamiltonian : not null access 
			       function (X : Real_Vector; N : Nat) return AD_Type;
			     Var	 : in     Variable;
			     Control     : in out Control_Type)
			    return Real_Vector is
      X   : Real_Vector renames Var.X;
      N   : Nat  renames Control.N;
      Dt  : Real renames Control.Dt;
      Err : Real renames Control.Err;
      J   : constant Sparse_Matrix    := -Omega (N);
      Old : constant Evaluation_Level :=  Level;
      K1, K2, K3, K4, Y, Z : Real_Vector (X'Range);
   begin
      pragma Assert (2 * N = Var.N2);
      -- Turn off the calculation of Hessians (not used for explicit schemes):
      Level := Gradient; 
      
      K1 := To_Array (J * Grad (Hamiltonian (X, N)));
      K2 := To_Array (J * Grad (Hamiltonian (X + (0.50 * Dt) * K1, N)));
      K3 := To_Array (J * Grad (Hamiltonian (X + (0.75 * Dt) * K2, N)));
      Y  := X + (Dt / 9.0) * (2.0 * K1 + 3.0 * K2 + 4.0 * K3);
      K4 := To_Array (J * Grad (Hamiltonian (Y, N)));
      Z  := X + (Dt / 24.0) * (7.0 * K1 + 6.0 * K2 + 8.0 * K3 + 3.0 * K4);

      Err := Norm (Z - Y);

      Level := Old; -- return to previous evaluation level
      return (Z);
   end Bogack_Shampine;
   
   
   procedure Update (Hamiltonian : not null access 
		    	       function (X : Real_Vector; N : Nat) return AD_Type;
		     Var         : in out Variable;
		     Control     : in out Control_Type) is
      use Real_Functions;
      X   : Real_Vector renames Var.X;
      T   : Real renames Var.T;
      N   : Nat  renames Control.N;
      Dt  : Real renames Control.Dt;
      Err : Real renames Control.Err;
      Eps : Real renames Control.Eps;
      Y   : Real_Vector (X'Range);
   begin
      pragma Assert (2 * N = Var.N2);
      
      Err := 1.0;
      while Err > Eps loop
	 Y := Bogack_Shampine (Hamiltonian, Var, Control);
	 if (Err <= Eps) then
	    X  := Y;
	    T  := T + Dt;
	 end if;
	 Dt := 0.8 * Dt * (Eps / (Err + 1.0e-20)) ** 0.3;
      end loop;
   end Update;
   
   
   procedure Print_XYZ (File : in File_Type;
			Var  : in Variable) is
      use Real_Functions, Real_IO;
      X :  Real_Vector renames Var.X;
      X1, Y1, X2, Y2 : Real;
      R : constant Real := 10.0;
   begin
      X1 := -R * Sin (X (1));
      Y1 :=  R * Cos (X (1));
      X2 := X1 - R * Sin (2.0 * X (1) + X (2));
      Y2 := Y1 + R * Cos (2.0 * X (1) + X (2));
      -- print header
      Put_Line (File, "3");
      Put (File, "Properties=pos:R:2   Time=");
      Put (File, Var.T, Aft => 5, Exp => 0);
      New_Line (File);
      -- position of ball 1
      Put_Line (File, "0.0     0.0     5.0");
      -- position of ball 2
      Put (File => File, Item => X1);
      Put (File => File, Item => "     ");
      Put (File => File, Item => Y1);
      Put (File => File, Item => "     ");
      Put (File => File, Item => "5.0");
      New_Line (File => File);
      -- position of ball 3
      Put (File => File, Item => X2);
      Put (File => File, Item => "     ");
      Put (File => File, Item => Y2);
      Put (File => File, Item => "     ");
      Put (File => File, Item => "5.0");
      New_Line (File => File);
   end Print_XYZ;
   
   --- print data ------
   procedure Print_Data (Var : in Variable;
			 Hamiltonian : not null access 
			   function (X : Real_Vector; N : Nat) return AD_Type) is
      use Real_Functions, Real_IO;
      T : Real renames Var.X (1);
      S : Real renames Var.X (2);
      X, Y : Real_Vector (1 .. 2);
   begin
      X (1) := -Sin (T);
      Y (1) :=  Cos (T);
      X (2) := X (1) - Sin (2.0 * T + S);
      Y (2) := Y (1) + Cos (2.0 * T + S);
      
      ---------------------------------------
      Put (Var.T, Aft => 6, Exp => 0); -- print time
      for I in 1 .. 2 loop
	 Put (",  "); Put (X (I), Aft => 4, Exp => 0);
	 Put (",  "); Put (Y (I), Aft => 4, Exp => 0); -- print positions
      end loop;
      Put (",  "); 
      -- print total energy
      Put (Val (Hamiltonian (Var.X, 2)), Aft => 10, Exp => 0); New_Line;
   end Print_Data;
   
   procedure Print_Data_L (File	: in File_Type;
			   Var	: in Variable) is
      use Real_Functions, Real_IO;
      T : Real renames Var.X (1);
      S : Real renames Var.X (2);
      X, Y : Real_Vector (1 .. 2);
   begin
      X (1) := -Sin (T);
      Y (1) :=  Cos (T);
      X (2) := X (1) - Sin (2.0 * T + S);
      Y (2) := Y (1) + Cos (2.0 * T + S);
      
      ---------------------------------------
      Put (Var.T, Aft => 6, Exp => 0); -- print time
      for I in 1 .. 2 loop
	 -- print positions
	 Put (",  "); Put (X (I), Aft => 4, Exp => 0);
	 Put (",  "); Put (Y (I), Aft => 4, Exp => 0);
      end loop;
      New_Line;
   end Print_Data_L;
   
   
end Auto_Differentiation.Integrator;
